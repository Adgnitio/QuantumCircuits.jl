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
import QuantumCircuits.QCircuits.Registers as Reg
using QuantumCircuits.QCircuits.Qiskit
using QuantumCircuits.QCircuits.Graph
using QuantumCircuits.QCircuits.Math
using QuantumCircuits.QCircuits.Qiskit: IS_QISKIT, NoQiskitError

using MacroTools

using LinearAlgebra

using CBOOCall: @cbooify

import QuantumCircuits.QCircuits.QBase: add!, tomatrix, setparameters!, simplify,
       standardGateError, decompose, measure!, bindparameters!, needbedecompsed
import Base.show

export QCircuit, getparameters, getRandParameters, toString, setClassicalRegister!, toPythonQiskitCode, ismeasured

"Nothing function"
const nop = () -> nothing

@doc raw"""
    QCircuit - Quantum Circuit structure
The fundamental element of quantum computing is the quantum circuit. A quantum circuit is a computational routine consisting of coherent quantum operations on quantum data, such as qubits. It is an ordered sequence of quantum gates, measurements and resets, which may be conditioned on real-time classical computation. A set of quantum gates is said to be universal if any unitary transformation of the quantum data can be efficiently approximated arbitrarily well as as sequence of gates in the set. Any quantum program can be represented by a sequence of quantum circuits and classical near-time computation.

In QuantumCircuits, this core element is represented by the QCircuit struct.
Below is an example of a quantum circuit that makes a three-qubit GHZ state defined as:
```math
|\psi\rangle = \frac{|000\rangle + |111\rangle}{\sqrt{2}}
```

```julia
using QuantumCircuits
qc = QCircuit(3)
qc.h(0)
qc.cx(0, 1)
qc.cx(1, 2)
```

```
      ┌───┐          
q0_0: ┤ H ├──■───────
      └───┘┌─┴─┐     
q0_1: ─────┤ X ├──■──
           └───┘┌─┴─┐
q0_2: ──────────┤ X ├
                └───┘
c0: 3/═══════════════
                     
```

"""
mutable struct QCircuit <: QuantumCircuit
    qubits::Int
    qRegisters::Vector{QuantumAbstractRegister}
    cRegisters::Vector{ClassicalRegister}
    dcg::DirectedGraph{Int}
    has_code::Bool
    code::Vector{QuantumGate}
    vqubits::Vector{Qubit}
    vcbits::Vector{Cbit}
    measures::Vector{Pair{Qubit, Cbit}}
    measures_matrix::Matrix{Float64}

    function QCircuit(qRegs::Vector{<:QuantumAbstractRegister}, cRegs::Vector{ClassicalRegister}; inline_optimization=true)
        n = sum([length(i) for i in qRegs])

        qc = new(n, qRegs, cRegs, DirectedGraph{Int}(n, inline_optimization), false, QuantumGate[], Qubit[], ClassicalRegister[], Pair{Qubit, Cbit}[], eye(2^n))

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

        return qc
    end
end
QCircuit(n::Integer; inline_optimization=true) = QCircuit([QuantumRegister(n)], [ClassicalRegister(n)], inline_optimization=inline_optimization)
QCircuit(qReg::QuantumAbstractRegister, cReg::ClassicalRegister; inline_optimization=true) = QCircuit([qReg], [cReg], inline_optimization=inline_optimization)
QCircuit(reg::QuantumAbstractRegister; inline_optimization=true) = QCircuit([reg], ClassicalRegister[ClassicalRegister(reg.tomeasure ? length(reg) : 0)], inline_optimization=inline_optimization)
QCircuit(regs::Vector{<:QuantumAbstractRegister}; inline_optimization=true) = QCircuit(regs, [ClassicalRegister(sum([(i.tomeasure ? length(i) : 0) for i in regs]))], inline_optimization=inline_optimization)
QCircuit(regs::Vector{<:QuantumAbstractRegister}, cReg::ClassicalRegister; inline_optimization=true) = QCircuit(regs, [cReg], inline_optimization=inline_optimization)
function QCircuit(qc::QCircuit)
    qregs = [QuantumRegister(length(r), r.name) for r in qc.qRegisters]
    cregs = [ClassicalRegister(length(r), r.name) for r in qc.cRegisters]

    QCircuit(qregs, cregs)
end

###################################################################################
function Base.getindex(qc::QCircuit, idx::Int)
    @assert idx >=0 && idx < qc.qubits "The $idx is out of bound for $(qc.qubits) circuits."

    return qc.vqubits[idx + 1]
end

"The size of register."
Base.length(qc::QCircuit) = length(qc.vqubits)


###################################################################################
@cbooify QCircuit (x, sx, y, z, h, cx, s, sdg, t, tdg, u, u3, rx, ry, rz, p, cp, swap, rzx,
                   u4, barrier, measure, add!, set!)

"Add function macro"
macro addfunction(name, gate)
    eval(quote
        $name(qc::QCircuit, args...) = add!(qc, $gate, args...)
    end)
end
macro addfunction1param(name, gate)
    eval(quote
        $name(qc::QCircuit, q, θ=ParameterT(rand()*2π)) = add!(qc, $gate, q, θ)
    end)
end
macro addfunction3param(name, gate)
    eval(quote
        $name(qc::QCircuit, q, θ=ParameterT(rand()*2π), ϕ=ParameterT(rand()*2π), λ=ParameterT(rand()*2π)) = add!(qc, $gate, q, θ, ϕ, λ)
    end)
end

@addfunction(x, X)
@addfunction(sx, Sx)
@addfunction(y, Y)
@addfunction(z, Z)
@addfunction(h, H)
@addfunction(cx, CX)
@addfunction(swap, Swap)
@addfunction(s, S)
@addfunction(sdg, Sd)
@addfunction(t, T)
@addfunction(tdg, Td)
@addfunction(u, U)

@addfunction1param(rx, Rx)
@addfunction1param(ry, Ry)
@addfunction1param(rz, Rz)
@addfunction1param(p, P)

@addfunction3param(u3, U3)

rzx(qc::QCircuit, q1, q2, θ=ParameterT(rand()*2π)) = add!(qc, Rzx, q1, q2, θ)
u4(qc::QCircuit, q1, q2, params=[ParameterT(rand() * 2π) for i in 1:U4_params]) = add!(qc, U4, q1, q2, params)
cp(qc::QCircuit, q1, q2, λ=ParameterT(rand()*2π)) = add!(qc, CP, q1, q2, λ)

barrier(qc::QCircuit) = add!(qc, Barrier(qc.vqubits))
measure(qc::QCircuit, q, c) = measure!(qc, q, c)
measure(qc::QCircuit) = measure!(qc)


function set!(qc::QCircuit, reg::QuantumInteger, num::Integer)
    @assert reg.state == Reg.Empty "Unable to set used register."
    @assert num < 2^reg.integer "The number $num is out of the register with $(reg.integer) qubits."

    # Set the value
    for (i, v) in enumerate(reverse(bitstring(num)))
        if v == '1'
            qc.x(reg[i-1])
        end
    end

    # The number was set
    reg.state = Reg.SettedNumber
    nothing
end

function set!(qc::QCircuit, reg::QuantumFloat, num::AbstractFloat)
    @assert reg.state == Reg.Empty "Unable to set used register."

    integer_part = Int(floor(num))
    fractional_part = (num - integer_part) * 2^reg.fractional
    fractional_part_int = Int(fractional_part)

    @assert integer_part < 2^reg.integer "The number $num is out of the register with $(reg.integer) qubits."
    @assert abs(fractional_part_int - fractional_part) < 1e-10 "The number $num is out of the register fractional part with $(reg.fractional) qubits."


    # Set the value
    for (i, v) in enumerate(reverse(bitstring(integer_part * 2^reg.fractional + fractional_part_int)))
        if v == '1'
            qc.x(reg[i-1])
        end
    end

    # The number was set
    reg.state = Reg.SettedNumber
    nothing
end

###################################################################################

function getCode(c::QCircuit)
    if c.has_code
        return c.code
    else
        c.code = to_vector(c.dcg)
        c.has_code = true

        return c.code
    end
end

Base.:(==)(c1::QCircuit, c2::QCircuit) = c1.qubits == c2.qubits && getCode(c1) == getCode(c2) && c1.measures_matrix == c2.measures_matrix
Base.hash(c::QCircuit, h::UInt) = hash((c.qubits, getCode(code), c.measures_matrix), h)

"Add gate to circuit"
function add!(qc::QCircuit, gate::QuantumGate)
    add!(qc.dcg, gate)
    qc.has_code = false
    nothing
end
function add!(qc::QCircuit, gate::Type{<:QuantumGate}, qubit::Qubit)
    @assert qubit in qc.vqubits "Qubits is out of circuit qubits range."

    add!(qc.dcg, gate(qubit))
    qc.has_code = false
    nothing
end
function add!(qc::QCircuit, gate::Type{<:QuantumGate}, args::Tuple)
    add!(qc.dcg, gate(args...))
    qc.has_code = false
    nothing
end
function add!(qc::QCircuit, gate::Type{<:QuantumGate}, qubit::Integer)
    @assert qubit < qc.qubits "Qubits $qubit is out of circuit qubits range $(qc.qubits)."

    add!(qc, gate, qc.vqubits[qubit + 1])

    nothing
end
function add!(qc::QCircuit, gate::Type{<:QuantumGate}, qubits::AbstractVector)
    for qubit in qubits
        add!(qc, gate, qubit)
    end
end
function add!(qc::QCircuit, gate::Type{<:QuantumGate}, qubit1::Qubit, qubit2::Qubit)
    @assert qubit1 in qc.vqubits "Qubits is out of circuit qubits range."
    @assert qubit2 in qc.vqubits "Qubits is out of circuit qubits range."
    @assert qubit1 != qubit2 "Qubits have to be different: $qubit1, $qubit2."

    add!(qc.dcg, gate(qubit1, qubit2))
    qc.has_code = false
    nothing
end
function add!(qc::QCircuit, gate::Type{<:QuantumGate}, qubit1::Integer, qubit2::Integer)
    @assert qubit1 < qc.qubits "Qubits $qubit1 is out of circuit qubits range $(qc.qubits)."
    @assert qubit2 < qc.qubits "Qubits $qubit2 is out of circuit qubits range $(qc.qubits)."
    @assert qubit1 != qubit2 "Qubits have to be different: $qubit1, $qubit2."

    add!(qc, gate, qc.vqubits[qubit1 + 1], qc.vqubits[qubit2 + 1])

    nothing
end
function add!(qc::QCircuit, gate::Type{<:QuantumGate}, qubit1::Qubit, qubit2::Qubit, params)
    @assert qubit1 in qc.vqubits "Qubits is out of circuit qubits range."
    @assert qubit2 in qc.vqubits "Qubits is out of circuit qubits range."
    @assert qubit1 != qubit2 "Qubits have to be different: $qubit1, $qubit2."

    if isnothing(params)
        add!(qc.dcg, gate(qubit1, qubit2))
    else
        add!(qc.dcg, gate(qubit1, qubit2, params))
    end
    qc.has_code = false
    nothing
end
function add!(qc::QCircuit, gate::Type{<:QuantumGate}, qubit1::Integer, qubit2::Integer, params)
    @assert qubit1 < qc.qubits "Qubits $qubit1 is out of circuit qubits range $(qc.qubits)."
    @assert qubit2 < qc.qubits "Qubits $qubit2 is out of circuit qubits range $(qc.qubits)."
    @assert qubit1 != qubit2 "Qubits have to be different: $qubit1, $qubit2."

    add!(qc, gate, qc.vqubits[qubit1 + 1], qc.vqubits[qubit2 + 1], params)

    nothing
end


function add!(qc::QCircuit, gate::Type{<:QuantumGate}, qubit::Qubit, θ)
    @assert qubit in qc.vqubits "Qubits is out of circuit qubits range."

    add!(qc.dcg, gate(qubit, θ))
    qc.has_code = false
    nothing
end
function add!(qc::QCircuit, gate::Type{<:QuantumGate}, qubit::Integer, θ)
    @assert qubit < qc.qubits "Qubits $qubit is out of circuit qubits range $(qc.qubits)."

    add!(qc, gate, qc.vqubits[qubit + 1], θ)

    nothing
end

function add!(qc::QCircuit, gate::Type{<:QuantumGate}, qubits::Vector, θ)
    for qubit in qubits
        add!(qc, gate, qubit, θ)
    end
end

function add!(qc::QCircuit, gate::Type{<:QuantumGate}, qubit::Qubit, θ, ϕ, λ)
    @assert qubit in qc.vqubits "Qubits is out of circuit qubits range."

    add!(qc.dcg, gate(qubit, θ, ϕ, λ))
    qc.has_code = false
    nothing
end
function add!(qc::QCircuit, gate::Type{<:QuantumGate}, qubit::Integer, θ, ϕ, λ)
    @assert qubit < qc.qubits "Qubits $qubit is out of circuit qubits range $(qc.qubits)."

    add!(qc, gate, qc.vqubits[qubit + 1], θ, ϕ, λ)

    nothing
end
function add!(qc::QCircuit, gate::Type{<:QuantumGate}, qubits::AbstractVector, θ, ϕ, λ)
    for qubit in qubits
        add!(qc, gate, qubit, ParameterT(θ), ParameterT(ϕ), ParameterT(λ))
    end
end
function add!(qc::QCircuit, gates::Vector{<:QuantumGate})
    for g in gates
        add!(qc, g)
    end
    nothing
end

"Add measures to circuit"
function measure!(qc::QCircuit, qubit::Qubit, cbit::Cbit, setMatrix::Bool=true)
    push!(qc.measures, qubit => cbit)

    if setMatrix
        setMeasureMatrix!(qc)
    end
    nothing
end

"Add measures to circuit"
function measure!(qc::QCircuit, qubit::Integer, cbit::Integer, setMatrix::Bool=true)
    @assert qubit < qc.qubits "Qubits $qubit is out of circuit qubits range $(qc.qubits)."
    @assert cbit < length(qc.vcbits) "Qubits $qubit is out of circuit qubits range $(length(qc.vcbits))."

    measure!(qc, qc.vqubits[qubit + 1], qc.vcbits[cbit + 1], setMatrix)
end
function measure!(qc::QCircuit, qubits::AbstractVector{T}, cbits::AbstractVector{V}) where {T, V}
    @assert length(qubits) == length(cbits) "The length of qubits and classical bits should be equal."

    for (i, j) in zip(eachindex(qubits), eachindex(cbits))
        measure!(qc, qubits[i], cbits[j], false)
    end

    setMeasureMatrix!(qc)
    nothing
end

"Add measures to circuit"
function measure!(qc::QCircuit)
    @assert length(qc.cRegisters) == 1 "Something go wrong."
    cr = qc.cRegisters[1]

    i = 0
    for r in qc.qRegisters
        if r.tomeasure
            l = length(r)
            qc.measure(r, cr[i:i+l-1])
            i = i + l
        end
    end
end


function setMeasureMatrix!(qc::QCircuit)
    measured = [(i-1, v) for (i, v) in enumerate(sort!([getid(q)  for (q, c) in qc.measures]))]
    measured_qubits = length(measured)
    mes_matrix = zeros(2^measured_qubits, 2^qc.qubits)
    for i in 0:(2^qc.qubits-1)
        bin = digits(i, base=2, pad=qc.qubits)

        # calculate new index
        ni = 0
        for (idx, m) in measured
            if bin[m + 1] == 1
                ni += 2^idx
            end
        end
        mes_matrix[ni+1, i+1] = 1.0
    end

    qc.measures_matrix = mes_matrix
end

function ismeasured(qc, reg)
    @assert reg in qc.qRegisters "Regiester has to be in circuits."

    # measured qubits
    mes_qubits = Set([getid(k[1]) for k in qc.measures])

    # check if all qubits were measured
    return all(getid(b) in mes_qubits for b in reg.bits)
end

###################################################################################

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
    if !IS_QISKIT
        throw(NoQiskitError("Unable to convert to Qiskit, qiskit is not installed."))
    end

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

"Function to convert QCircuit to Python code using Qiskit library."
function toPythonQiskitCode(circuit::QCircuit)
    code = ""

    args = ""
    for (i, r) in enumerate(circuit.qRegisters)
        if isnothing(r.name)
            code *= "qr$i = QuantumRegister($(length(r.bits)))\n"
        else
            code *= "qr$i = QuantumRegister($(length(r.bits)), \"$(r.name)\")\n"
        end
        if i > 1
            args *= ", "
        end
        args *= "qr$i"
    end

    if !isempty(circuit.cRegisters)
        args *= ", "
        for (i, r) in enumerate(circuit.cRegisters)
            if isnothing(r.name)
                code *= "cr$i = ClassicalRegister($(length(r.bits)))\n"
            else
                code *= "cr$i = ClassicalRegister($(length(r.bits)), \"$(r.name)\")\n"
            end
            if i > 1
                args *= ", "
            end
            args *= "cr$i"
        end
    end

    code *= "qc = QuantumCircuit($args)\n"

    for gate in getCode(circuit)
        code *= getPythonCode("qc", gate)
    end

    for (q, c) in circuit.measures
        code *= "qc.measure($q, $(getid(c)))\n"
    end

    return code
end

"Show method"
function Base.show(io::IO, qc::QCircuit)
    if IS_QISKIT
        print(io, toQiskit(qc))
    else
        print(io, "QCircuit")
    end
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
    unitary = LinearAlgebra.I * (1.0 + 0im)
    for gate in getCode(qc)
        l = length(gate)
        if !isnothing(params) && l > 0
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
        if !isnothing(params) && l > 0
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

    # add measrments
    for (q, c) in qc.measures        
        newqc.measure(getid(q), getid(c))
   end

    return newqc
end

function decompose(gate::CP)
#      ┌────────┐
# q_0: ┤ P(λ/2) ├──■───────────────■────────────
#      └────────┘┌─┴─┐┌─────────┐┌─┴─┐┌────────┐
# q_1: ──────────┤ X ├┤ P(-λ/2) ├┤ X ├┤ P(λ/2) ├
#                └───┘└─────────┘└───┘└────────┘   
    return [
        P(gate.control, gate.λ/2),
        CX(gate.control, gate.target),
        P(gate.target, -gate.λ/2),
        CX(gate.control, gate.target),
        P(gate.target, gate.λ/2)
    ]
end

function decompose(gate::Swap)
#           ┌───┐     
# q_0: ──■──┤ X ├──■──
#      ┌─┴─┐└─┬─┘┌─┴─┐
# q_1: ┤ X ├──■──┤ X ├
#      └───┘     └───┘ 
    return [
        CX(gate.control, gate.target),
        CX(gate.target, gate.control),
        CX(gate.control, gate.target)        
    ]
end


needbedecompsed(qc::QCircuit) = any(needbedecompsed(gate) for gate in getCode(qc))
needbedecompsed(::CP) = true
needbedecompsed(::Swap) = true


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
