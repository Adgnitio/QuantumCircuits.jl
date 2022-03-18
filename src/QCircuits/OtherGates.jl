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

module OtherGates

using QuantumCircuits.QCircuits.QBase
using QuantumCircuits.QCircuits.Math
using QuantumCircuits.QCircuits.Gates
using QuantumCircuits.QCircuits.Gates: Parameter, Xmatrix, Zmatrix
using QuantumCircuits.QCircuits.Registers

import Base: show, length, inv
import QuantumCircuits.QCircuits.QBase: tomatrix, setparameters!, simplify,
       standardGateError, decompose

import QuantumCircuits.QCircuits.Gates: getqubits, getqubitsids, appendparams!,
       appendRandParams!, setparameters!, getArgs, decompose

export Rzx



"Rotation gates"
abstract type TwoQubitsRotationGate <: QuantumGate end

"Rotation macro"
macro twoQubitsRotationGate(name)
    eval(quote
        mutable struct $name <: TwoQubitsRotationGate
            qubit1::Qubit
            qubit2::Qubit
            θ::Parameter
        end
        $name(qubit1::Qubit, qubit2::Qubit) = $name(qubit1, qubit2, nothing)

        Base.:(==)(g1::$name, g2::$name) = g1.qubit1 == g2.qubit1 && g1.qubit2 == g2.qubit2 && g1.θ == g2.θ
        Base.hash(g::$name, h::UInt) = hash((g.qubit1, g.qubit2, g.θ), h)

        #$name(qubis::Vector{T}, θ::ParamT) where T <: Integer = [$name(q, θ) for q in qubis]
        getArgs(g::$name) = (getid(g.qubit1), getid(g.qubit2), g.θ)

        inv(g::$name) = $name(g.qubit1, g.qubit2, -g.θ)
    end)
end


@twoQubitsRotationGate(Rzx)


"Return the qubits on which operate the gate"
getqubits(gate::TwoQubitsRotationGate) = (gate.qubit1, gate.qubit2)

"Return the qubits ids on which operate the gate"
getqubitsids(gate::TwoQubitsRotationGate) = (getid(gate.qubit1), getid(gate.qubit2))

Base.show(io::IO, gate::TwoQubitsRotationGate) = print(io, "$(typeof(gate))($(gate.qubit1), $(gate.qubit2), $(gate.θ))")

function appendparams!(param, gate::TwoQubitsRotationGate)
    @assert gate.θ != nothing "Unable to append param without set it first."
    push!(param, gate.θ)
end

appendRandParams!(param, gate::TwoQubitsRotationGate) = append!(param, rand(1) * 2π)

Base.length(gate::TwoQubitsRotationGate) = 1

function setparameters!(gate::TwoQubitsRotationGate, params)
    @assert length(params) == 1 "Rzx has 1 parameters but $(length(param)) was provided."
    @assert params[1] != nothing "Unable to set nothing as parameters"
    gate.θ = params[1]
end

function tomatrix(g::Rzx, param=nothing)
    if param == nothing
        θ = g.θ
    else
        @assert length(param) == 1 "Rzx has 1 parameters but $(length(param)) was provided."
        θ = param[1]
    end

    q1 = getid(g.qubit1)
    q2 = getid(g.qubit2)
    Δ = abs(q1 - q2)

    if q1 < q2
        mat = kron(Xmatrix, eye(2^(Δ-1)), Zmatrix)
    else
        mat = kron(Zmatrix, eye(2^(Δ-1)), Xmatrix)
    end

    return exp(-im * (θ/2) * mat)
end


# function decompose(gate::U4)
#     return [
#         U3(gate.qubit1, gate.params[1], gate.params[2], gate.params[3]),
#         U3(gate.qubit2, gate.params[4], gate.params[5], gate.params[6]),
#         CX(gate.qubit2, gate.qubit1),
#         Rz(gate.qubit1, gate.params[7]),
#         Ry(gate.qubit2, gate.params[8]),
#         CX(gate.qubit1, gate.qubit2),
#         Ry(gate.qubit2, gate.params[9]),
#         CX(gate.qubit2, gate.qubit1),
#         U3(gate.qubit1, gate.params[10], gate.params[11], gate.params[12]),
#         U3(gate.qubit2, gate.params[13], gate.params[14], gate.params[15])
#     ]
# end


end  # module OtherGates
