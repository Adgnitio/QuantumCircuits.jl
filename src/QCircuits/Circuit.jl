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

module Circuit

using QuantumCircuits.QCircuits.QBase
using QuantumCircuits.QCircuits.Gates
using QuantumCircuits.QCircuits.Gates: ParamT, appendparams!, appendRandParams!
using QuantumCircuits.QCircuits.ComplexGates
using QuantumCircuits.QCircuits.ComplexGates: U4_params
using QuantumCircuits.QCircuits.OtherGates
using QuantumCircuits.QCircuits.Instructions
using QuantumCircuits.QCircuits.Registers
using QuantumCircuits.QCircuits.Qiskit
using QuantumCircuits.QCircuits.Graph
using QuantumCircuits.QCircuits.Math

using LinearAlgebra: I

import QuantumCircuits.QCircuits.QBase: add!, tomatrix, setparameters!, simplify,
       standardGateError, decompose, measure!, bindparameters!
import Base.show

export QCircuit, getparameters, getRandParameters, toString, setClassicalRegister!

"Nothing function"
const nop = () -> nothing

"Quantum circuit."
mutable struct QCircuit <: QuantumCircuit
    qubits::Integer
    qRegisters::Vector{QuantumRegister}
    cRegisters::Vector{ClassicalRegister}
    dcg::DirectedGraph
    has_code::Bool
    code::Vector{QuantumGate}
    vqubits::Vector{Qubit}
    vcbits::Vector{Cbit}
    measures::Vector{Pair{Qubit, Cbit}}
    measures_matrix::Matrix{Float64}

    # gates
    x::Function
    sx::Function
    y::Function
    z::Function
    h::Function
    cx::Function
    s::Function
    sdg::Function
    t::Function
    tdg::Function
    u::Function
    u3::Function
    rx::Function
    ry::Function
    rz::Function
    rzx::Function
    u4::Function
    barrier::Function
    measure::Function

    function QCircuit(qRegs::Vector{QuantumRegister}, cRegs::Vector{ClassicalRegister})
        n = sum([length(i) for i in qRegs])

        qc = new(n, qRegs, cRegs, DirectedGraph(n), false, QuantumGate[], Qubit[], ClassicalRegister[], Pair{Qubit, Cbit}[], eye(2^n), nop)

        # assign the qregister to circuit
        for r in qRegs
            append!(qc.vqubits, r.bits)
        end
        for (i, q) in enumerate(qc.vqubits)
            setid(q, i - 1)
        end

        # and classical one
        for r in cRegs
            append!(qc.vcbits, r.bits)
        end
        for (i, q) in enumerate(qc.vcbits)
            setid(q, i - 1)
        end

        # add gates
        qc.x = q -> add!(qc, X, q)
        qc.sx = q -> add!(qc, Sx, q)
        qc.y = q -> add!(qc, Y, q)
        qc.z = q -> add!(qc, Z, q)
        qc.h = q -> add!(qc, H, q)
        qc.cx = (c, t) -> add!(qc, CX, c, t)
        qc.s = q -> add!(qc, S, q)
        qc.sdg = q -> add!(qc, Sd, q)
        qc.t = q -> add!(qc, T, q)
        qc.tdg = q -> add!(qc, Td, q)
        qc.u = (q, θ, ϕ, λ) -> add!(qc, U, q, θ, ϕ, λ)
        qc.u3 = (q, θ=ParameterT(rand()*2π), ϕ=ParameterT(rand()*2π), λ=ParameterT(rand()*2π)) -> add!(qc, U3, q, θ, ϕ, λ)
        qc.rx = (q, θ=ParameterT(rand()*2π)) -> add!(qc, Rx, q, θ)
        qc.ry = (q, θ=ParameterT(rand()*2π)) -> add!(qc, Ry, q, θ)
        qc.rz = (q, θ=ParameterT(rand()*2π)) -> add!(qc, Rz, q, θ)
        qc.rzx = (q1, q2, θ=ParameterT(rand()*2π)) -> add!(qc, Rzx, q1, q2, θ)
        qc.u4 = (q1, q2, params=[ParameterT(rand() * 2π) for i in 1:U4_params]) -> add!(qc, U4, q1, q2, params)
        qc.barrier = () -> add!(qc, Barrier(qc.vqubits))
        qc.measure = (q, c) -> measure!(qc, q, c)

        return qc
    end
end
QCircuit(n::Integer) = QCircuit([QuantumRegister(n)], [ClassicalRegister(n)])
QCircuit(qReg::QuantumRegister, cReg::ClassicalRegister) = QCircuit([qReg], [cReg])
QCircuit(reg::QuantumRegister) = QCircuit([reg], ClassicalRegister[])
QCircuit(regs::Vector{QuantumRegister}) = QCircuit(regs, [ClassicalRegister(sum([length(i) for i in regs]))])
QCircuit(regs::Vector{QuantumRegister}, cReg::ClassicalRegister) = QCircuit(regs, [cReg])
function QCircuit(qc::QCircuit)
    qregs = [QuantumRegister(length(r), r.name) for r in qc.qRegisters]
    cregs = [ClassicalRegister(length(r), r.name) for r in qc.cRegisters]

    QCircuit(qregs, cregs)
end

function getCode(c::QCircuit)
    if c.has_code
        return c.code
    else
        c.code = to_vector(c.dcg)
        c.has_code = true

        return c.code
    end
end

Base.:(==)(c1::QCircuit, c2::QCircuit) = c1.qubits == c2.qubits && getCode(c1) == getCode(c2)
Base.hash(c::QCircuit, h::UInt) = hash((c.qubits, getCode(code)), h)

"Add gate to circuit"
function add!(qc::QCircuit, gate::QuantumGate)
    add!(qc.dcg, gate)
    qc.has_code = false
    nothing
end
function add!(qc::QCircuit, gate::Type{T}, qubit::Qubit) where T <: QuantumGate
    @assert qubit in qc.vqubits "Qubits is out of circuit qubits range."

    add!(qc.dcg, gate(qubit))
    qc.has_code = false
    nothing
end
function add!(qc::QCircuit, gate::Type{T}, args::Tuple) where T <: QuantumGate
    add!(qc.dcg, gate(args...))
    qc.has_code = false
    nothing
end
function add!(qc::QCircuit, gate::Type{T}, qubit::Integer) where T <: QuantumGate
    @assert qubit < qc.qubits "Qubits $qubit is out of circuit qubits range $(qc.qubits)."

    add!(qc, gate, qc.vqubits[qubit + 1])

    nothing
end
function add!(qc::QCircuit, gate::Type{T}, qubits::Vector) where T <: QuantumGate
    for qubit in qubits
        add!(qc, gate, qubit)
    end
end
function add!(qc::QCircuit, gate::Type{T}, qubits::QuantumRegister) where T <: QuantumGate
    for qubit in qubits
        add!(qc, gate, qubit)
    end
end
function add!(qc::QCircuit, gate::Type{T}, qubit1::Qubit, qubit2::Qubit) where T <: QuantumGate
    @assert qubit1 in qc.vqubits "Qubits is out of circuit qubits range."
    @assert qubit2 in qc.vqubits "Qubits is out of circuit qubits range."

    add!(qc.dcg, gate(qubit1, qubit2))
    qc.has_code = false
    nothing
end
function add!(qc::QCircuit, gate::Type{T}, qubit1::Integer, qubit2::Integer) where T <: QuantumGate
    @assert qubit1 < qc.qubits "Qubits $qubit1 is out of circuit qubits range $(qc.qubits)."
    @assert qubit2 < qc.qubits "Qubits $qubit2 is out of circuit qubits range $(qc.qubits)."

    add!(qc, gate, qc.vqubits[qubit1 + 1], qc.vqubits[qubit2 + 1])

    nothing
end
function add!(qc::QCircuit, gate::Type{T}, qubit1::Qubit, qubit2::Qubit, params) where T <: QuantumGate
    @assert qubit1 in qc.vqubits "Qubits is out of circuit qubits range."
    @assert qubit2 in qc.vqubits "Qubits is out of circuit qubits range."

    if params == nothing
        add!(qc.dcg, gate(qubit1, qubit2))
    else
        add!(qc.dcg, gate(qubit1, qubit2, params))
    end
    qc.has_code = false
    nothing
end
function add!(qc::QCircuit, gate::Type{T}, qubit1::Integer, qubit2::Integer, params) where T <: QuantumGate
    @assert qubit1 < qc.qubits "Qubits $qubit1 is out of circuit qubits range $(qc.qubits)."
    @assert qubit2 < qc.qubits "Qubits $qubit2 is out of circuit qubits range $(qc.qubits)."

    add!(qc, gate, qc.vqubits[qubit1 + 1], qc.vqubits[qubit2 + 1], params)

    nothing
end


function add!(qc::QCircuit, gate::Type{T}, qubit::Qubit, θ) where T <: QuantumGate
    @assert qubit in qc.vqubits "Qubits is out of circuit qubits range."

    add!(qc.dcg, gate(qubit, θ))
    qc.has_code = false
    nothing
end
function add!(qc::QCircuit, gate::Type{T}, qubit::Integer, θ) where T <: QuantumGate
    @assert qubit < qc.qubits "Qubits $qubit is out of circuit qubits range $(qc.qubits)."

    add!(qc, gate, qc.vqubits[qubit + 1], θ)

    nothing
end
function add!(qc::QCircuit, gate::Type{T}, qubits::Vector, θ) where T <: QuantumGate
    for qubit in qubits
        add!(qc, gate, qubit, θ)
    end
end

function add!(qc::QCircuit, gate::Type{T}, qubit::Qubit, θ, ϕ, λ) where T <: QuantumGate
    @assert qubit in qc.vqubits "Qubits is out of circuit qubits range."

    add!(qc.dcg, gate(qubit, θ, ϕ, λ))
    qc.has_code = false
    nothing
end
function add!(qc::QCircuit, gate::Type{T}, qubit::Integer, θ, ϕ, λ) where T <: QuantumGate
    @assert qubit < qc.qubits "Qubits $qubit is out of circuit qubits range $(qc.qubits)."

    add!(qc, gate, qc.vqubits[qubit + 1], θ, ϕ, λ)

    nothing
end
function add!(qc::QCircuit, gate::Type{T}, qubits::QuantumRegister, θ, ϕ, λ) where T <: QuantumGate
    for qubit in qubits
        add!(qc, gate, qubit, ParameterT(θ), ParameterT(ϕ), ParameterT(λ))
    end
end
function add!(qc::QCircuit, gate::Type{T}, qubits::Vector, θ, ϕ, λ) where T <: QuantumGate
    for qubit in qubits
        add!(qc, gate, qubit, ParameterT(θ), ParameterT(ϕ), ParameterT(λ))
    end
end

function add!(qc::QCircuit, gates::Vector{T}) where T <: QuantumGate
    for g in gates
        add!(qc, g)
    end
    nothing
end

"Add measuers to circuit"
function measure!(qc::QCircuit, qubit::Qubit, cbit::Cbit, setMatrix::Bool=true)
    push!(qc.measures, qubit => cbit)

    if setMatrix
        setMeasureMatrix!(qc)
    end
end
"Add measuers to circuit"
function measure!(qc::QCircuit, qubit::Integer, cbit::Integer, setMatrix::Bool=true)
    @assert qubit < qc.qubits "Qubits $qubit is out of circuit qubits range $(qc.qubits)."
    @assert cbit < length(qc.vcbits) "Qubits $qubit is out of circuit qubits range $(length(qc.vcbits))."

    measure!(qc, qc.vqubits[qubit + 1], qc.vcbits[cbit + 1], setMatrix)
end
function measure!(qc::QCircuit, qubits::Vector{T}, cbits::Vector{T}) where T <: Integer
    n = length(qubits)
    @assert n == length(cbits) "The length of qubits and classical bits should be equal."

    for i in 1:n
        measure!(qc, qubits[i], cbits[i], false)
    end

    setMeasureMatrix!(qc)
end

function setMeasureMatrix!(qc::QCircuit)
    measuered = [(i-1, v) for (i, v) in enumerate(sort([getid(q)  for (q, c) in qc.measures]))]
    measuered_qubits = length(measuered)
    mes_matrix = zeros(2^measuered_qubits, 2^qc.qubits)
    for i in 0:(2^qc.qubits-1)
        bin = digits(i, base=2, pad=qc.qubits)

        # calculate new index
        ni = 0
        for (idx, m) in measuered
            if bin[m + 1] == 1
                ni += 2^idx
            end
        end
        mes_matrix[ni+1, i+1] = 1.0
    end

    qc.measures_matrix = mes_matrix
end

"Add the classical register to circuit."
function setClassicalRegister!(qc::QCircuit, cr::ClassicalRegister)
    @assert qc.cRegisters == ClassicalRegister[] "Unable to set register."
    @assert qc.vcbits == Cbit[] "Internal error :P"

    push!(qc.cRegisters, cr)
    append!(qc.vcbits, cr.bits)
    for (i, q) in enumerate(qc.vcbits)
        setid(q, i - 1)
    end
end

"Function to convert QCircuit to Qiskit circuit"
function toQiskit(circuit::QCircuit)
    if isempty(circuit.cRegisters)
        qc = QiskitCircuit(circuit.qRegisters, nothing)
    else
        qc = QiskitCircuit(circuit.qRegisters, circuit.cRegisters)
    end

    for gate in getCode(circuit)
        add!(qc, gate)
    end

    for (q, c) in circuit.measures
        measure!(qc, q, c)
    end

    return qc
end

"Show method"
function Base.show(io::IO, qc::QCircuit)
    print(io, toQiskit(qc))
end

"Convert circuit to string"
function toString(qc::QCircuit)
    io = IOBuffer()
    show(io, qc)
    return String(take!(io))
end

"Function create new parameters vector"
function getparameters(qc::QCircuit)
    params = ParamT[]
    for gate in getCode(qc)
        appendparams!(params, gate)
    end

    return params
end

"Function create new random parameters vector"
function getRandParameters(qc::QCircuit)
    params = ParamT[]
    for gate in getCode(qc)
        appendRandParams!(params, gate)
    end

    return params
end

"Save the parameters to the quantum circuits."
function setparameters!(qc::QCircuit, params)
    i = 1
    for gate in getCode(qc)
        l = length(gate)
        if l > 0
            @views setparameters!(gate, params[i:(i+l-1)])
            i = i + l
        end
    end

    @assert i - 1 == length(params) "There is different number of parameters in circuit $(i-1) and in vector $(length(params))."
end

"Function bind the parameters, so now there only values."
function bindparameters!(qc::QCircuit)
    for gate in getCode(qc)
        bindparameters!(gate)
    end
end

function tomatrix(qc::QCircuit, params=nothing)
    i = 1
    unitary = I * (1.0 + 0im)
    for gate in getCode(qc)
        l = length(gate)
        if params != nothing && l > 0
            @views gate_unitary = tomatrix(qc.qubits, gate, params[i:(i+l-1)])
            # TODO gate_unitary = tomatrix(qc.qubits, gate, params[i:(i+l-1)])
            i = i + l
        else
            gate_unitary = tomatrix(qc.qubits, gate)
        end

        #println("GU: $gate_unitary")
        unitary = gate_unitary * unitary
        #println("U: $unitary")
    end
    @assert i - 1 == length(params) "There is different number of parameters in circuit $(i-1) and in vector $(length(params))."

    return unitary
end

function standardGateError(qc::QCircuit, params=nothing)
    err = 0.0
    i = 1
    for gate in getCode(qc)
        l = length(gate)
        if params != nothing && l > 0
            @views gerr = standardGateError(gate, params[i:(i+l-1)])
            i = i + l
        else
            gerr = standardGateError(gate)
        end

        err +=  gerr
    end

    #return err / (100 * length(qc.code))
    return err
end

function Base.inv(qc::QCircuit)
    newqc = QCircuit(qc)

    for gate in reverse(getCode(qc))
        newgate = inv(gate)
        add!(newqc, typeof(newgate), getArgs(newgate)...)
    end

    return newqc
end

################################################################################

function decompose(qc::QCircuit)
    newqc = QCircuit(qc)

    for gate in getCode(qc)
        newgates = decompose(gate)
        for newgate in newgates
            add!(newqc, typeof(newgate), getArgs(newgate)...)
        end
    end

    return newqc
end

function simplify(qc::QCircuit)
    newqc = QCircuit(qc)

    for gate in getCode(qc)
        newgates = simplify(gate)
        for newgate in newgates
            add!(newqc, typeof(newgate), getArgs(newgate)...)
        end
    end

    return newqc
end

################################################################################

function Base.append!(qc::QCircuit, apqc::QCircuit)
    @assert length(qc.measures) == 0 "Unable append to measured circuit."
    @assert length(apqc.measures) == 0 "Unable append measured circuit."
    @assert length(qc.vqubits) == length(apqc.vqubits) "Unable append, the qubits are different"

    for (q, aq) in zip(qc.vqubits, apqc.vqubits)
        @assert q == aq "There are difference in qubits."
    end

    for c in getCode(apqc)
        add!(qc, typeof(c), getArgs(c)...)
    end
end

end  # module Circuit
