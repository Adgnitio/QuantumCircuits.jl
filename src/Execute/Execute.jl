# This code is part of QuantumCircuits.
#
# (C) Copyright Rafał Pracht 2022.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

module Execute

using PyCall

using QuantumCircuits.QCircuits.Math
using QuantumCircuits.QCircuits.QBase
using QuantumCircuits.QCircuits.QBase: QuantumDevice, tomatrix
using QuantumCircuits.QCircuits.Circuit
using QuantumCircuits.QCircuits.Circuit: toQiskit
using QuantumCircuits.QCircuits.Qiskit
using QuantumCircuits.QCircuits.Qiskit: qiskit
using QuantumCircuits.QCircuits.Gates: ParamT
using QuantumCircuits.QCircuits.Registers

export execute, QuantumSimulator, QiskitQuantum, @ket_str, @bra_str,
       loss_expected_zero_state, qderivative, qexecute

"Quantum Simulator"
struct QuantumSimulator <: QuantumDevice end

"Qiskit Quantum Device"
struct QiskitQuantum <: QuantumDevice
    backend::PyObject
end
QiskitQuantum() = QiskitQuantum(qiskit.providers.aer.QasmSimulator())

######################################################

# Initial stata
const _0 = [1+0im; 0]
const _1 = [0im; 1]

const qubit_map = Dict('0' => _0, '1' => _1)

#const _00 = kron(_0, _0)

function str2state(s)
    @assert all([c in ('0', '1') for c in s]) "Only 0 and 1 are allowed."

    # Reverse the bit order
    s = reverse(s)
    state = qubit_map[s[1]]

    if length(s) == 1
        return state
    end

    for d in s[2:end]
        state = kron(qubit_map[d], state)
    end

    return state
end

"Macro proces ket binary vector to quantum state"
macro ket_str(s)
    return str2state(s)
end

"Macro proces ket binary vector to quantum state"
macro bra_str(s)
    return dagger(str2state(s))
end

"Method convert state to probability"
state2probability(s) = abs.(s) .^ 2

"Function execute the quantum circuit on simulator"
function execute(device::QuantumSimulator, qc::QCircuit, params=nothing; prev_matrix=nothing)
    # The initial, all zeros vector
    s = vcat([1], zeros(2^qc.qubits-1))

    m = tomatrix(qc, params)
    if prev_matrix != nothing
        s = prev_matrix * s
    end

    p = state2probability(m * s)
    return qc.measures_matrix * p
end


function extractProbability(counts::Dict, qubits::Integer, shots::Integer)
    p = zeros(2^qubits)

    for (k, n) in counts
        p[parse(Int, k; base=2) + 1] = n/shots
    end

    return p
end

function execute(device::QiskitQuantum, qc::QiskitCircuit; shots=8192, debug=true)
    job = qiskit.execute(qc.qc, device.backend, shots=shots)
    debug && println("Job ID: $(job.job_id()) run single circuit.")
    #debug && qiskit.tools.monitor.job_monitor(job)

    res = job.result().get_counts()

    return extractProbability(res, length(qc.cbits), shots)
end

function execute(device::QuantumSimulator, circuits::Vector{QCircuit})
    results = Vector{Float64}[]
    for qc in circuits
        push!(results, execute(device, qc))
    end

    return results
end

function execute(device::QiskitQuantum, circuits::Vector{QiskitCircuit}; shots=8192, debug=true)
    jobs = [qc.qc for qc in circuits]
    job = qiskit.execute(jobs, device.backend, shots=shots)
    debug && println("Job ID: $(job.job_id()) run $(length(jobs)) circuits.")
    debug && qiskit.tools.monitor.job_monitor(job)

    # Get results
    job_results = job.result().get_counts()

    results = Vector{Float64}[]
    for (res, qc) in zip(job_results, circuits)
        push!(results, extractProbability(res, length(qc.cbits), shots))
    end

    return results
end

######################################################

function generate_mesuere_circuit(circuit::QCircuit, i)
    qr = QuantumRegister(length(circuit.vqubits), "q")
    cr = ClassicalRegister(length(circuit.vcbits))
    qc = QCircuit(qr, cr)

    measureQubits = [getid(q) for (q, b) in circuit.measures]

    for (i, v) in enumerate(reverse(string(i, base=2, pad=length(circuit.vcbits))))
        if v == '1'
            qc.x(measureQubits[i])
        end
    end

    for (q, c) in circuit.measures
        qc.measure(qr[getid(q)], cr[getid(c)])
    end
    return qc
end

function generate_mesuere_circuits(circuit::QCircuit)
    circuits = QiskitCircuit[]
    for i in 0:(2^length(circuit.vcbits)-1)
        push!(circuits, toQiskit(generate_mesuere_circuit(circuit, i)))
    end

    return circuits
end

function correctMeasures(correctMatrix, state)
    nstate = correctMatrix * state
    nstate = min.(max.(nstate, 0.0), 1.0)

    return nstate / sum(nstate)
end

######################################################

"Loss method to check if final state is zero."
function loss_expected_zero_state(state)
    return -log(real(state[1])+1e-32) + sum([real(v)^2 for v in state[2:end]])
end

######################################################

"Set the parameters to the cicquit and convert to qiskit"
function setAndConvert(qc, params::Vector{T}) where T <: ParamT
    setparameters!(qc, params)
    return toQiskit(qc)
end

"Function return the shift"
function getShift(n, idx)
    shift = zeros(ParamT, n)
    shift[idx] = π/2

    return shift
end

"Caclulate the jacobian on device."
function qjacobian(backend, qc, params::Vector{T}, runN=1, corrMes=true) where T <: ParamT
    @assert (length(params) * runN * 2 + runN) < 300 "To mamny run."
    circuits = QiskitCircuit[]

    if corrMes
        mc = generate_mesuere_circuits(qc)
        mcN = length(mc)
        append!(circuits, mc)
    end

    # function call
    rqc = setAndConvert(qc, params)
    for i in 1:runN
        push!(circuits, rqc)
    end

    # shift rule for all params
    n = length(params)
    for i in 1:n
        shift = getShift(n, i)
        slqc = setAndConvert(qc, params + shift)
        srqc = setAndConvert(qc, params - shift)
        for i in 1:runN
            push!(circuits, slqc)
        end
        for i in 1:runN
            push!(circuits, srqc)
        end
    end

    results = execute(backend, circuits)

    # correct the measures
    if corrMes
        corrMes = hcat(results[1:mcN]...)
        correctMeas = inv(corrMes)
        results = [correctMeasures(correctMeas, r) for r in results[(mcN+1):end]]
        #results = results[(mcN+1):end]
    end

    f_results = results[1:runN]
    df_results_tmp = results[(runN+1):end]


    if runN > 1
        df_results = Vector{Float64}[]
        df = zeros(ParamT, length(df_results_tmp[1]))
        for i in 1:length(df_results_tmp)
            df = df + df_results_tmp[i]
            if i % runN == 0
                push!(df_results, df / runN)
                df = zeros(ParamT, length(df_results_tmp[1]))
            end
        end
        push!(df_results, df / runN)
    else
        df_results = df_results_tmp
    end

    # proces result
    mat = (df_results[1] - df_results[2]) / 2
    for i in 2:n
        idx = 2 * i
        m = (df_results[idx - 1] - df_results[idx]) / 2
        mat = hcat(mat, m)
    end

    return sum(f_results)/runN, mat
end

"Caclulate loss value and it derivateive"
function qderivative(backend, qc, loss, params::Vector{T}, runN=1) where T <: ParamT
    y, jak = qjacobian(backend, qc, params, runN)
    der = transpose(jak) * loss'(y)
    return loss(y), der
end

"Caclulate the the value."
function qexecute(backend, qc, loss, params::Vector{T}, runN=1, corrMes=true) where T <: ParamT

    # set the parameters
    setparameters!(qc, params)
    # convert to qiskit
    qqc = toQiskit(qc)

    if corrMes
        circuits = generate_mesuere_circuits(qc)
        mcN = length(circuits)

        # add the orginal one
        push!(circuits, qqc)

        results = execute(backend, circuits)

        # correct the measures
        corrMes = hcat(results[1:mcN]...)
        correctMeas = inv(corrMes)

        ϕ = correctMeasures(correctMeas, results[end])
    else
        ϕ = execute(backend, qqc)
    end

    return loss(ϕ)
end

end  # module Execute
