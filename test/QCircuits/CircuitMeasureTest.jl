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
using QuantumCircuits.Execute

qc1 = QCircuit(4)
qc1.x(2)
qc1.h(1)
qc1.cx(3, 2)
qc1.measure([1,2,3], [1,2,3])

qc2 = QCircuit(4)
qc2.x(2)
qc2.h(1)
qc2.cx(3, 2)
qc2.measure(1:3, 1:3)

@test qc1 == qc2

# Create circuit with barrier
qc1 = QCircuit(4)
qc1.x(2)
qc1.h(1)
qc1.barrier()
qc1.cx(3, 2)
