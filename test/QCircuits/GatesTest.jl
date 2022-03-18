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
using QuantumCircuits.QCircuits.Gates: appendparams!, ParamT, appendRandParams!,
                                       U3

# using QGen.Quantum.QBase
# using QGen.Quantum.Gates
# using QGen.Quantum.Registers
# using QGen.Quantum.Circuit
# using QGen.Quantum.Gates: appendparams!, ParamT, ParameterT, appendRandParams!,
#                           setparameters!

function getparam(g)
    param = ParamT[]
    appendparams!(param, g)

    return param
end
function getRandParam(g)
    param = ParamT[]
    appendRandParams!(param, g)

    return param
end


qr = QuantumRegister(1)
g = U3(qr[0], 1.0, 2.0, 3.0)
@test getparam(g) == ParamT[]
@test getRandParam(g) == ParamT[]
@test length(g) == 0
@test_throws AssertionError setparameters!(g, [2.0, 1.0, 4.0])
@test_throws AssertionError setparameters!(g, [2.0])
setparameters!(g, [])


g = U3(qr[0], ParameterT(1.0), ParameterT(2.0), ParameterT(3.0))
@test getparam(g) == ParamT[1.0, 2.0, 3.0]
@test length(getRandParam(g)) == 3
@test getparam(g) != getRandParam(g)
@test length(g) == 3
setparameters!(g, [2.0, 1.0, 4.0])
@test getparam(g) == ParamT[2.0, 1.0, 4.0]
@test_throws AssertionError setparameters!(g, [2.0, 1.0, 4.0, 4.0])
@test_throws AssertionError setparameters!(g, [2.0, 1.0])
setparameters!(g, [-2.0, -1.0, -4.0])
@test getparam(g) == ParamT[-2.0, -1.0, -4.0]


g = U3(qr[0])
@test length(getparam(g)) == 3
@test length(getRandParam(g)) == 3
@test getparam(g) != getRandParam(g)
@test length(g) == 3
setparameters!(g, [2.0, 1.0, 4.0])
@test getparam(g) == ParamT[2.0, 1.0, 4.0]


qr = QuantumRegister(1)
g = U3(qr[0], ParameterT(1.0), ParameterT(2.0), 3.0)
@test getparam(g) == ParamT[1.0, 2.0]
@test length(getRandParam(g)) == 2
@test getparam(g) != getRandParam(g)
@test length(g) == 2
setparameters!(g, [2.0, 1.0])
@test getparam(g) == ParamT[2.0, 1.0]
@test_throws AssertionError setparameters!(g, [2.0, 1.0, 4.0])


################################################################################
#                Complex Params                                                #
################################################################################
using QuantumCircuits.QCircuits.ComplexGates
using QuantumCircuits.QCircuits.ComplexGates: U4_params

qr = QuantumRegister(2)
g = U4(qr[0], qr[1])
p = getparam(g)
@test length(p) == U4_params
rp = getRandParam(g)
@test length(rp) == U4_params
@test p != rp
@test length(g) == U4_params
np = [i for i in 1:U4_params]
@test p != np
setparameters!(g, np)
@test getparam(g) == np


qc = QCircuit(2)
qc.u4(0, 1)
p = getparameters(qc)
@test length(p) == U4_params
rp = getRandParameters(qc)
@test length(rp) == U4_params
@test p != rp
np = [i for i in 1:U4_params]
@test p != np
setparameters!(qc, np)
@test getparameters(qc) == np


dqc = decompose(qc)
@test length(getparameters(dqc)) == U4_params
@test getparameters(dqc) == np
np = [-i for i in 1:U4_params]
setparameters!(dqc, np)
@test getparameters(dqc) == np
bindparameters!(dqc)
@test getparameters(dqc) == Number[]

# @test length(getparameters(qc)) == U4_params
# bindparameters!(qc)
# @test getparameters(qc) == Number[]


################################################################################

qr = QuantumRegister(3, "q")
qc = QCircuit(qr)
qc.u3(qr)
n = 3*3
p = getparameters(qc)
@test length(p) == n
rp = getRandParameters(qc)
@test length(rp) == n
@test p != rp
setparameters!(qc, rp)
@test getparameters(qc) == rp

np = [i for i in 1:n]
@test rp != np
setparameters!(qc, np)
@test getparameters(qc) == np

bindparameters!(qc)
@test getparameters(qc) == Number[]
