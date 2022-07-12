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
using QuantumCircuits.Circuit: getCode
using QuantumCircuits.QML

################################################################################
#  U3                                                                          #
################################################################################
# Test the graph u3 optimization
expqc = QCircuit(1)
expqc.u3(0, ParameterT(1.0), ParameterT(2.0), ParameterT(3.0))


qc = QCircuit(1)
qc.u3(0, ParameterT(1.0), ParameterT(2.0), ParameterT(3.0))
qc.u3(0, ParameterT(4.0), ParameterT(5.0), ParameterT(6.0))
@test qc == expqc
@test length(getCode(qc)) == 1

# with random parameters
qc = QCircuit(1)
qc.u3(0)
qc.u3(0)
@test length(getparameters(qc)) == 3
@test length(getCode(qc)) == 1

# No optimization in case wthout parameters
qc = QCircuit(1)
qc.u3(0, 1.0, 2.0, 3.0)
qc.u3(0, 4.0, 5.0, 6.0)
@test length(getparameters(qc)) == 0
@test length(getCode(qc)) == 2



################################################################################
#  Remove gate                                                                 #
################################################################################
using QuantumCircuits.QCircuits.Graph: remove!

function dcg_rm!(qc, idx)
    remove!(qc.dcg, getCode(qc)[idx])
    qc.has_code = false
end

qc = QCircuit(1)
qc.x(0)
qc.h(0)

dcg_rm!(qc, 2)
expqc = QCircuit(1, inline_optimization=false)
expqc.x(0)
@test expqc == qc
@test length(getCode(qc)) == 1

dcg_rm!(qc, 1)
expqc = QCircuit(1, inline_optimization=false)
@test expqc == qc
@test length(getCode(qc)) == 0

qc.y(0)
qc.z(0)
expqc = QCircuit(1, inline_optimization=false)
expqc.y(0)
expqc.z(0)
@test expqc == qc
@test length(getCode(qc)) == 2

dcg_rm!(qc, 2)
expqc = QCircuit(1, inline_optimization=false)
expqc.y(0)
@test expqc == qc
@test length(getCode(qc)) == 1

qc.s(0)
expqc = QCircuit(1, inline_optimization=false)
expqc.y(0)
expqc.s(0)
@test expqc == qc
@test length(getCode(qc)) == 2


###################################################

qc = QCircuit(4)
qc.cx(0, 1)
qc.cx(2, 3)
qc.cx(1, 2)
qc.cx(0, 1)
qc.cx(2, 3)
@test length(getCode(qc)) == 5

dcg_rm!(qc, 5)
expqc = QCircuit(4, inline_optimization=false)
expqc.cx(0, 1)
expqc.cx(2, 3)
expqc.cx(1, 2)
expqc.cx(0, 1)
@test expqc == qc
@test length(getCode(qc)) == 4

dcg_rm!(qc, 4)
expqc = QCircuit(4, inline_optimization=false)
expqc.cx(0, 1)
expqc.cx(2, 3)
expqc.cx(1, 2)
@test expqc == qc
@test length(getCode(qc)) == 3

dcg_rm!(qc, 3)
expqc = QCircuit(4, inline_optimization=false)
expqc.cx(0, 1)
expqc.cx(2, 3)
@test expqc == qc
@test length(getCode(qc)) == 2

dcg_rm!(qc, 2)
dcg_rm!(qc, 1)
expqc = QCircuit(4, inline_optimization=false)
@test expqc == qc
@test length(getCode(qc)) == 0

qc.cx(1, 0)
qc.cx(3, 2)
qc.cx(2, 1)
qc.cx(1, 0)
qc.cx(3, 2)

expqc = QCircuit(4, inline_optimization=false)
expqc.cx(1, 0)
expqc.cx(3, 2)
expqc.cx(2, 1)
expqc.cx(1, 0)
expqc.cx(3, 2)
@test expqc == qc
@test length(getCode(qc)) == 5

################################################################################
#  CX                                                                          #
################################################################################
qc = QCircuit(2)
qc.cx(0, 1)

expqc = QCircuit(2, inline_optimization=false)
expqc.cx(0, 1)
@test expqc == qc

# Add next CX
qc.cx(0, 1)
expqc = QCircuit(2, inline_optimization=false)
@test expqc == qc

qc.cx(1, 0)
expqc = QCircuit(2, inline_optimization=false)
expqc.cx(1, 0)
@test expqc == qc

qc.cx(1, 0)
expqc = QCircuit(2, inline_optimization=false)
@test expqc == qc


qc.cx(0, 1)
qc.x(0)
qc.cx(0, 1)

expqc = QCircuit(2, inline_optimization=false)
expqc.cx(0, 1)
expqc.x(0)
expqc.cx(0, 1)
@test expqc == qc

################################################################################
qc = QCircuit(2)
qc.cx(0, 1)
qc.cx(1, 0)

expqc = QCircuit(2, inline_optimization=false)
expqc.cx(0, 1)
expqc.cx(1, 0)
@test expqc == qc
@test length(getCode(qc)) == 2

