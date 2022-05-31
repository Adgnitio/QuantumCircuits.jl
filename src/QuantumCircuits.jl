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

module QuantumCircuits

include("QCircuits/QCircuits.jl")
include("QML/QML.jl")
include("Execute/Execute.jl")
include("Simulation/Simulation.jl")

using QuantumCircuits.QCircuits.QBase
using QuantumCircuits.QCircuits.Registers
using QuantumCircuits.QCircuits.Circuit

export QuantumCircuit, QuantumDevice, QuantumGate, add!, QCircuit,
       QuantumRegister, ClassicalRegister, tomatrix, decompose, simplify, @gate, @circ

end # module
