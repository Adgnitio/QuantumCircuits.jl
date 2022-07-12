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
using QuantumCircuits.QCircuits.QLib: qft_rotations!, swap_registers!, qft!


################################################################################
#  qft_rotations                                                               #
################################################################################
expqc = QCircuit(3)
expqc.h(2)
expqc.cp(1, 2, π/2)
expqc.cp(0, 2, π/4)
expqc.h(1)
expqc.cp(0, 1, π/2)
expqc.h(0)

qc = QCircuit(3)
code = QuantumGate[]
qft_rotations!(qc, code, qc.qubits)
qc.add!(code)
@test qc == expqc



################################################################################
#  qft                                                                         #
################################################################################
qc = QCircuit(3)
qft!(qc, doswap=false)
@test qc == expqc

## with swap
expqc.swap(0, 2)

qc = QCircuit(3)
qft!(qc)
@test qc == expqc

# inverse
expqc = QCircuit(3)
expqc.swap(0, 2)
expqc.h(0)
expqc.cp(0, 1, -π/2)
expqc.h(1)
expqc.cp(0, 2, -π/4)
expqc.cp(1, 2, -π/2)
expqc.h(2)

qc = QCircuit(3)
qft!(qc, inverse=true)
@test qc == expqc

################################################################################
#  qft - registers                                                             #
################################################################################
expqr1 = QuantumRegister(3)
expqr2 = QuantumRegister(3)
expqc = QCircuit([expqr1, expqr2])
expqc.h(expqr1[2])
expqc.cp(expqr1[1], expqr1[2], π/2)
expqc.cp(expqr1[0], expqr1[2], π/4)
expqc.h(expqr1[1])
expqc.cp(expqr1[0], expqr1[1], π/2)
expqc.h(expqr1[0])


qr1 = QuantumRegister(3)
qr2 = QuantumRegister(3)
qc = QCircuit([qr1, qr2])
code = QuantumGate[]
qft_rotations!(qr1, code, length(qr1))
qc.add!(code)
@test qc == expqc


#####
qr1 = QuantumRegister(3)
qr2 = QuantumRegister(3)
qc = QCircuit([qr1, qr2])
qft!(qc, qr1, doswap=false)
@test qc == expqc

## with swap
expqc.swap(qr1[0], qr1[2])

qr1 = QuantumRegister(3)
qr2 = QuantumRegister(3)
qc = QCircuit([qr1, qr2])
qft!(qc, qr1)
@test qc == expqc



# inverse
expqr1 = QuantumRegister(3)
expqr2 = QuantumRegister(3)
expqc = QCircuit([expqr1, expqr2])
expqc.swap(qr1[0], qr1[2])
expqc.h(expqr1[0])
expqc.cp(expqr1[0], expqr1[1], -π/2)
expqc.h(expqr1[1])
expqc.cp(expqr1[0], expqr1[2], -π/4)
expqc.cp(expqr1[1], expqr1[2], -π/2)
expqc.h(expqr1[2])


qr1 = QuantumRegister(3)
qr2 = QuantumRegister(3)
qc = QCircuit([qr1, qr2])
qft!(qc, qr1, inverse=true)
@test qc == expqc


################################################################################
#  qft - identity                                                              #
################################################################################
using QuantumCircuits.Execute

qr = QuantumRegister(6)
qc = QCircuit(qr)
qft!(qc, qr)
qft!(qc, qr, inverse=true)
qc.measure()

rs = getresults(qc)
@test abs(rs[qr]["000000"] - 1.0) < 1e-8
