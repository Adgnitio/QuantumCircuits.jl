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
expqc.cp(0, 2, π/4)
expqc.cp(1, 2, π/2)
expqc.h(1)
expqc.cp(0, 1, π/2)
expqc.h(0)

qc = QCircuit(3)
qft_rotations!(qc, 3)
@test qc == expqc


################################################################################
#  qft                                                                         #
################################################################################
expqc = QCircuit(3)
expqc.h(2)
expqc.cp(0, 2, π/4)
expqc.cp(1, 2, π/2)
expqc.h(1)
expqc.cp(0, 1, π/2)
expqc.h(0)

qc = QCircuit(3)
qft!(qc, 3, doswap=false)
@test qc == expqc

###
expqc = QCircuit(3)
expqc.h(2)
expqc.cp(0, 2, π/4)
expqc.cp(1, 2, π/2)
expqc.h(1)
expqc.cp(0, 1, π/2)
expqc.h(0)
expqc.swap(0, 2)

qc = QCircuit(3)
qft!(qc, 3)
@test qc == expqc





################################################################################


# n = 10
# qc = QCircuit(n)
# qft!(qc, 3)


# @test qc == expqc


# decompose(qc)