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


using Test

using QuantumCircuits
using QuantumCircuits.QML
using QuantumCircuits.QML.Optimization
using QuantumCircuits.Execute



# using QGen.Quantum
# using QGen.Quantum.Gates
# using QGen.Quantum.QBase
# using QGen.Quantum.Circuit
# using QGen.Quantum.Math
# using QGen.Genetic.CircOpt
# using QGen.Math.Stat
# using QGen.Quantum.Execute
# using QGen.Quantum.Execute: str2state, state2probability

const SingOptBackend = QuantumSimulator()

# Start in minimum and singular point
n = 2
qc = QCircuit(n)
qc.ry(0, ParameterT(π/2))
qc.ry(1, ParameterT(π/4))
qc.cx(0, 1)
qc.measure([0, 1], [0, 1])

exp = execute(SingOptBackend, qc)
loss(params) = sum(abs.(exp - execute(SingOptBackend, qc, params)))
dloss(params) = real(loss'(params))

#####################

pa = getparameters(qc)
pa[1] += 0.1
val, _, iter = gradientDescent(loss, dloss, pa, α=1.0, maxItr=1000, isExpectedZero=true, argsArePeriodic=true, ϵ=1e-12)
#val, _, iter = gradientDescent(loss, dloss, pa, α=1.0, maxItr=10, isExpectedZero=true, argsArePeriodic=true, ϵ=1e-12, debug=true)
@test abs(val) < 1e-10


#####################
pa = getRandParameters(qc)
val, _, iter = gradientDescent(loss, dloss, pa, α=1.0, maxItr=1000, isExpectedZero=true, argsArePeriodic=true, ϵ=1e-12)
@test abs(val) < 1e-10
