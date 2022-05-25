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

module ComplexGates

using QuantumCircuits.QCircuits.QBase
using QuantumCircuits.QCircuits.Gates
using QuantumCircuits.QCircuits.Gates: Parameter
using QuantumCircuits.QCircuits.Registers

import QuantumCircuits.QCircuits.Gates: getqubits, getqubitsids, appendparams!,
       appendRandParams!, setparameters!, getArgs, decompose, bindparameters!

export U4

const U4_params = 4*3 + 3

"Universary two qubits gate."
struct U4 <: QuantumGate
    qubit1::Qubit
    qubit2::Qubit
    params::Vector{Parameter}
end
U4(qubit1::Qubit, qubit2::Qubit) = U4(qubit1, qubit2, [ParameterT(rand() * 2π) for i in 1:U4_params])

Base.:(==)(g1::U4, g2::U4) = g1.qubit1 == g2.qubit1 && g1.qubit2 == g2.qubit2 && g1.params == g2.params
Base.hash(g::U4, h::UInt) = hash((g.qubit1, g.qubit2, g.params), h)

getArgs(g::U4) = (getid(g.qubit1), getid(g.qubit2), params)

"Return the qubits on which operate the gate"
getqubits(gate::U4) = (gate.qubit1, gate.qubit2)

"Return the qubits ids on which operate the gate"
getqubitsids(gate::U4) = (getid(gate.qubit1), getid(gate.qubit2))

Base.show(io::IO, gate::U4) = print(io, "U4($(gate.qubit1), $(gate.qubit2))")

function appendparams!(param, gate::U4)
    for p in gate.params
        @assert !isnothing(p) "Unable to append param without set it first."
        appendparams!(param, p)
    end
end

# function bindparameters!(gate::U4)
#     for p in gate.params
#         @assert p != nothing "Unable to append param without set it first."
#         appendparams!(param, p)
#     end
#     gate.θ = getvalue(gate.θ)
#     gate.ϕ = getvalue(gate.ϕ)
#     gate.λ = getvalue(gate.λ)
# end
#

appendRandParams!(param, gate::U4) = append!(param, rand(U4_params) * 2π)

Base.length(::U4) = U4_params

function setparameters!(gate::U4, params)
    @assert length(params) == U4_params "U4 has $(U4_params) parameters but $(length(param)) was provided."
    for i in 1:U4_params
        gate.params[i] = ParameterT(params[i])
    end
end


function decompose(gate::U4)
    return [
        U3(gate.qubit1, gate.params[1], gate.params[2], gate.params[3]),
        U3(gate.qubit2, gate.params[4], gate.params[5], gate.params[6]),
        CX(gate.qubit2, gate.qubit1),
        Rz(gate.qubit1, gate.params[7]),
        Ry(gate.qubit2, gate.params[8]),
        CX(gate.qubit1, gate.qubit2),
        Ry(gate.qubit2, gate.params[9]),
        CX(gate.qubit2, gate.qubit1),
        U3(gate.qubit1, gate.params[10], gate.params[11], gate.params[12]),
        U3(gate.qubit2, gate.params[13], gate.params[14], gate.params[15])
    ]
end


end  # module ComplexGates
