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

println("Execute tests from: MathTest")
include("MathTest.jl")
println("Execute tests from: RegistersTest")
include("RegistersTest.jl")
println("Execute tests from: CircuitTest")
include("CircuitTest.jl")
println("Execute tests from: CircuitMacroTest")
include("CircuitMacroTest.jl")
println("Execute tests from: InstructionsTest")
include("InstructionsTest.jl")
println("Execute tests from: CircuitSimplifyTest")
include("CircuitSimplifyTest.jl")
println("Execute tests from: CircuitUnitaryTest")
include("CircuitUnitaryTest.jl")
println("Execute tests from: GatesTest")
include("GatesTest.jl")
println("Execute tests from: OtherGatesTest")
include("OtherGatesTest.jl")
println("Execute tests from: CircuitMeasureTest")
include("CircuitMeasureTest.jl")
println("Execute tests from: GraphTest")
include("GraphTest.jl")
println("Execute tests from: CircuitPythonTest")
include("CircuitPythonTest.jl")
println("Execute tests from: CircuitProbTest")
include("CircuitProbTest.jl")

