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
using QuantumCircuits.QCircuits.Circuit: getCode, toQiskit
using QuantumCircuits.QCircuits.Math


n = 7
qr = QuantumRegister(n)
qc = QCircuit(qr)
qc.h(qr)
qc.barrier()
qc.x(qr)

@test length(getCode(qc)) == 2*n+1

eqr = QuantumRegister(n)
eqc = QCircuit(eqr)
eqc.h(eqr)
eqc.x(eqr)


toQiskit(qc)
@test unitary_error(tomatrix(qc), tomatrix(eqc)) < 1e-16
@test unitary_error(tomatrix(inv(qc)), tomatrix(inv(eqc))) < 1e-16
@test unitary_error(tomatrix(decompose(qc)), tomatrix(decompose(eqc))) < 1e-16
@test unitary_error(tomatrix(simplify(qc)), tomatrix(simplify(eqc))) < 1e-16
