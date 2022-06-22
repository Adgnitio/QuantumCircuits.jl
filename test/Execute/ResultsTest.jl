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
using QuantumCircuits.Execute.Results: Result, ResultsSet, getresults


# X
qc = QCircuit(2)
qc.x(0)
rs = getresults(qc)
@test abs(rs["01"].p - 1.00) < 1e-8


# Hadamards
qc = QCircuit(2)
qc.h([0, 1])
rs = getresults(qc)
@test abs(rs["00"].p - 0.25) < 1e-8
@test abs(rs["01"].p - 0.25) < 1e-8
@test abs(rs["10"].p - 0.25) < 1e-8
@test abs(rs["11"].p - 0.25) < 1e-8


# Registers
qr1 = QuantumRegister(2)
qr2 = QuantumRegister(2)
cr = ClassicalRegister(2+2)
qc = QCircuit([qr1, qr2], cr)
qc.x(qr1[1])
qc.h(qr2)
qc.measure(qr1, cr[2:3])
qc.measure(qr2, cr[0:1])
rs = getresults(qc)
@test abs(rs["1000"].p - 0.25) < 1e-8
@test abs(rs["1001"].p - 0.25) < 1e-8
@test abs(rs["1010"].p - 0.25) < 1e-8
@test abs(rs["1011"].p - 0.25) < 1e-8
rs1 = rs[qr1]
@test abs(rs1["00"].p - 0.25) < 1e-8
@test abs(rs1["01"].p - 0.25) < 1e-8
@test abs(rs1["10"].p - 0.25) < 1e-8
@test abs(rs1["11"].p - 0.25) < 1e-8
rs2 = rs[qr2]
@test abs(rs2["10"].p - 1.0) < 1e-8

# Registers 2
qr1 = QuantumRegister(2)
qr2 = QuantumRegister(2)
qr3 = QuantumRegister(2)
cr = ClassicalRegister(2+2)
qc = QCircuit([qr1, qr2, qr3], cr)
qc.x(qr1[1])
qc.h(qr2)
qc.measure(qr1, cr[2:3])
qc.measure(qr2, cr[0:1])
rs = getresults(qc)
@test abs(rs["1000"].p - 0.25) < 1e-8
@test abs(rs["1001"].p - 0.25) < 1e-8
@test abs(rs["1010"].p - 0.25) < 1e-8
@test abs(rs["1011"].p - 0.25) < 1e-8
rs1 = rs[qr1]
@test abs(rs1["00"].p - 0.25) < 1e-8
@test abs(rs1["01"].p - 0.25) < 1e-8
@test abs(rs1["10"].p - 0.25) < 1e-8
@test abs(rs1["11"].p - 0.25) < 1e-8
rs2 = rs[qr2]
@test abs(rs2["10"].p - 1.0) < 1e-8

# Registers 3 
qr1 = QuantumRegister(2)
qr2 = QuantumRegister(2)
qr3 = QuantumRegister(2)
cr = ClassicalRegister(2+2)
qc = QCircuit([qr1, qr2, qr3], cr)
qc.x(qr1[1])
qc.h(qr2)
qc.h(qr3)
qc.measure(qr1, cr[2:3])
qc.measure(qr2, cr[0:1])
rs = getresults(qc)
@test abs(rs["1000"].p - 0.25) < 1e-8
@test abs(rs["1001"].p - 0.25) < 1e-8
@test abs(rs["1010"].p - 0.25) < 1e-8
@test abs(rs["1011"].p - 0.25) < 1e-8
rs1 = rs[qr1]
@test abs(rs1["00"].p - 0.25) < 1e-8
@test abs(rs1["01"].p - 0.25) < 1e-8
@test abs(rs1["10"].p - 0.25) < 1e-8
@test abs(rs1["11"].p - 0.25) < 1e-8
rs2 = rs[qr2]
@test abs(rs2["10"].p - 1.0) < 1e-8



