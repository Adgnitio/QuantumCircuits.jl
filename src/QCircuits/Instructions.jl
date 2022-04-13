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

module Instructions

using QuantumCircuits.QCircuits.QBase
using QuantumCircuits.QCircuits.Gates
using QuantumCircuits.QCircuits.Gates: ParamT
using QuantumCircuits.QCircuits.Registers
using QuantumCircuits.QCircuits.Math

import Base: show, length, inv

import QuantumCircuits.QCircuits.Gates: getqubits, getqubitsids, appendparams!,
       appendRandParams!, setparameters!, getArgs

import QuantumCircuits.QCircuits.QBase: tomatrix, simplify, standardGateError,
       decompose

export Barrier

"Universary two qubits gate."
struct Barrier <: QuantumGate
    qubits::Vector{Qubit}
end

getArgs(g::Barrier) = ((g.qubits, ), )

"Return the qubits on which operate the gate"
getqubits(::Barrier) = ([getid(qubit) for qubit in g.qubits]..., )

"Return the qubits ids on which operate the gate"
getqubitsids(gate::Barrier) = [getid(qubit) for qubit in gate.qubits]

Base.show(io::IO, ::Barrier) = print(io, "Barrier")

decompose(gate::Barrier) = [gate]

tomatrix(gate::Barrier) = eye(2^length(gate.qubits))

inv(gate::Barrier) = gate

"The simplify the gates."
simplify(gate::Barrier) = [gate]

end  # module Instructions
