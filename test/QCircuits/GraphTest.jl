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
