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
using QuantumCircuits.Execute.Results: ResultsSet, getresults


# X
qc = QCircuit(2)
qc.x(0)
rs = getresults(qc)
@test abs(rs["01"] - 1.00) < 1e-8
## X
qr = QuantumRegister(2)
cr = ClassicalRegister(2)
qc = QCircuit(qr, cr)
qc.x(0)
qc.measure(qr, cr)
rs = getresults(qc)
@test abs(rs["01"] - 1.00) < 1e-8




# Hadamards
qc = QCircuit(2)
qc.h([0, 1])
rs = getresults(qc)
@test abs(rs["00"] - 0.25) < 1e-8
@test abs(rs["01"] - 0.25) < 1e-8
@test abs(rs["10"] - 0.25) < 1e-8
@test abs(rs["11"] - 0.25) < 1e-8


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
@test abs(rs["1000"] - 0.25) < 1e-8
@test abs(rs["1001"] - 0.25) < 1e-8
@test abs(rs["1010"] - 0.25) < 1e-8
@test abs(rs["1011"] - 0.25) < 1e-8
rs1 = rs[qr2]
@test abs(rs1["00"] - 0.25) < 1e-8
@test abs(rs1["01"] - 0.25) < 1e-8
@test abs(rs1["10"] - 0.25) < 1e-8
@test abs(rs1["11"] - 0.25) < 1e-8
rs2 = rs[qr1]
@test abs(rs2["10"] - 1.0) < 1e-8

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
@test abs(rs["1000"] - 0.25) < 1e-8
@test abs(rs["1001"] - 0.25) < 1e-8
@test abs(rs["1010"] - 0.25) < 1e-8
@test abs(rs["1011"] - 0.25) < 1e-8
rs1 = rs[qr2]
@test abs(rs1["00"] - 0.25) < 1e-8
@test abs(rs1["01"] - 0.25) < 1e-8
@test abs(rs1["10"] - 0.25) < 1e-8
@test abs(rs1["11"] - 0.25) < 1e-8
rs2 = rs[qr1]
@test abs(rs2["10"] - 1.0) < 1e-8


# Registers 2a
qr1 = QuantumRegister(2)
qr2 = QuantumRegister(2)
qr3 = QuantumRegister(2)
cr = ClassicalRegister(2+2)
qc = QCircuit([qr3, qr1, qr2], cr)
qc.x(qr1[1])
qc.h(qr2)
qc.measure(qr1, cr[2:3])
qc.measure(qr2, cr[0:1])
rs = getresults(qc)
@test abs(rs["1000"] - 0.25) < 1e-8
@test abs(rs["1001"] - 0.25) < 1e-8
@test abs(rs["1010"] - 0.25) < 1e-8
@test abs(rs["1011"] - 0.25) < 1e-8
rs1 = rs[qr2]
@test abs(rs1["00"] - 0.25) < 1e-8
@test abs(rs1["01"] - 0.25) < 1e-8
@test abs(rs1["10"] - 0.25) < 1e-8
@test abs(rs1["11"] - 0.25) < 1e-8
rs2 = rs[qr1]
@test abs(rs2["10"] - 1.0) < 1e-8


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
@test abs(rs["1000"] - 0.25) < 1e-8
@test abs(rs["1001"] - 0.25) < 1e-8
@test abs(rs["1010"] - 0.25) < 1e-8
@test abs(rs["1011"] - 0.25) < 1e-8
rs1 = rs[qr2]
@test abs(rs1["00"] - 0.25) < 1e-8
@test abs(rs1["01"] - 0.25) < 1e-8
@test abs(rs1["10"] - 0.25) < 1e-8
@test abs(rs1["11"] - 0.25) < 1e-8
rs2 = rs[qr1]
@test abs(rs2["10"] - 1.0) < 1e-8


# Registers =======================
qr1 = QuantumRegister(3)
qr2 = QuantumRegister(3)
qr3 = QuantumRegister(3)
cr = ClassicalRegister(3+3+3)
qc = QCircuit([qr1, qr2, qr3], cr)
qc.x(qr1[1])
qc.h(qr2)
qc.x(qr3[0])
qc.measure(qr1, cr[0:2])
qc.measure(qr2, cr[3:5])
qc.measure(qr3, cr[6:8])
rs = getresults(qc)
@test abs(rs["001000010"] - 0.125) < 1e-8
@test abs(rs["001001010"] - 0.125) < 1e-8
@test abs(rs["001010010"] - 0.125) < 1e-8
@test abs(rs["001011010"] - 0.125) < 1e-8
@test abs(rs["001100010"] - 0.125) < 1e-8
@test abs(rs["001101010"] - 0.125) < 1e-8
@test abs(rs["001110010"] - 0.125) < 1e-8
@test abs(rs["001111010"] - 0.125) < 1e-8
rs1 = rs[qr1]
@test abs(rs1["010"] - 1.0) < 1e-8
rs2 = rs[qr2]
@test abs(rs2["000"] - 0.125) < 1e-8
@test abs(rs2["001"] - 0.125) < 1e-8
@test abs(rs2["010"] - 0.125) < 1e-8
@test abs(rs2["011"] - 0.125) < 1e-8
@test abs(rs2["100"] - 0.125) < 1e-8
@test abs(rs2["101"] - 0.125) < 1e-8
@test abs(rs2["110"] - 0.125) < 1e-8
@test abs(rs2["111"] - 0.125) < 1e-8
rs3 = rs[qr3]
@test abs(rs3["001"] - 1.0) < 1e-8

# Registers =======================
qr1 = QuantumRegister(3)
qr2 = QuantumRegister(3)
qr3 = QuantumRegister(3)
cr = ClassicalRegister(3+3+3)
qc = QCircuit([qr1, qr2, qr3], cr)
qc.x(qr1[1])
qc.h(qr2)
qc.x(qr3[0])
qc.measure(qr1, cr[3:5])
qc.measure(qr2, cr[0:2])
qc.measure(qr3, cr[6:8])
rs = getresults(qc)
@test abs(rs["001010000"] - 0.125) < 1e-8
@test abs(rs["001010001"] - 0.125) < 1e-8
@test abs(rs["001010010"] - 0.125) < 1e-8
@test abs(rs["001010011"] - 0.125) < 1e-8
@test abs(rs["001010100"] - 0.125) < 1e-8
@test abs(rs["001010101"] - 0.125) < 1e-8
@test abs(rs["001010110"] - 0.125) < 1e-8
@test abs(rs["001010111"] - 0.125) < 1e-8
rs1 = rs[qr1]
@test abs(rs1["010"] - 1.0) < 1e-8
rs2 = rs[qr2]
@test abs(rs2["000"] - 0.125) < 1e-8
@test abs(rs2["001"] - 0.125) < 1e-8
@test abs(rs2["010"] - 0.125) < 1e-8
@test abs(rs2["011"] - 0.125) < 1e-8
@test abs(rs2["100"] - 0.125) < 1e-8
@test abs(rs2["101"] - 0.125) < 1e-8
@test abs(rs2["110"] - 0.125) < 1e-8
@test abs(rs2["111"] - 0.125) < 1e-8
rs3 = rs[qr3]
@test abs(rs3["001"] - 1.0) < 1e-8

# Registers =======================
qr1 = QuantumRegister(2)
qr2 = QuantumRegister(3)
qr3 = QuantumRegister(4)
cr = ClassicalRegister(2+3+4)
qc = QCircuit([qr1, qr2, qr3], cr)
qc.h(qr1)
qc.x(qr2[1])
qc.x(qr3[0])
qc.x(qr3[2])
qc.measure(qr1, cr[4:5])
qc.measure(qr2, cr[6:8])
qc.measure(qr3, cr[0:3])
rs = getresults(qc)
@test abs(rs["010000101"] - 0.25) < 1e-8
@test abs(rs["010010101"] - 0.25) < 1e-8
@test abs(rs["010100101"] - 0.25) < 1e-8
@test abs(rs["010110101"] - 0.25) < 1e-8

