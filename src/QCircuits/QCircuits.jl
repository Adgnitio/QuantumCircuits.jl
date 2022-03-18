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

module QCircuits

include("Math.jl")
include("QBase.jl")
include("Registers.jl")
include("Gates.jl")
include("Instructions.jl")
include("OtherGates.jl")
include("ComplexGates.jl")
include("Graph.jl")
include("Qiskit.jl")
include("Circuit.jl")

# using QGen.Quantum.QBase
#
# export QuantumCircuit, QuantumDevice, QuantumGate, add!

end # module
