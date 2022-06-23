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
using QuantumCircuits.Circuit

###############################################################
qra = QuantumRegister(2, "a")
qrb = QuantumRegister(2, "b")
cr = ClassicalRegister(4, "c")
qc = QCircuit([qra, qrb], cr)
qc.h(qra[0])
qc.z(qrb[0])
qc.cx(qra[0], qra[1])
qc.cx(qrb[0], qrb[1])
qc.cx(qra[0], qrb[1])

ex = "     ┌───┐          \na_0: ┤ H ├──■────■──
     └───┘┌─┴─┐  │  \na_1: ─────┤ X ├──┼──
     ┌───┐└───┘  │  \nb_0: ┤ Z ├──■────┼──
     └───┘┌─┴─┐┌─┴─┐
b_1: ─────┤ X ├┤ X ├
          └───┘└───┘
c: 4/═══════════════
                    "

@test ex == string(qc)

@test !ismeasured(qc, qra)
@test !ismeasured(qc, qrb)
qc.measure(qra, cr[0:1])
@test ismeasured(qc, qra)
@test !ismeasured(qc, qrb)
qc.measure(qrb, cr[2:3])
@test ismeasured(qc, qra)
@test ismeasured(qc, qrb)

###############################################################

using QuantumCircuits.QCircuits.Gates

qr = QuantumRegister(2, "q")
cr = ClassicalRegister(2, "c")
qc = QCircuit(qr, cr, inline_optimization=false)
add!(qc, X, qr[0])
add!(qc, H, qr[1])
add!(qc, Y, 0)
add!(qc, Z, 1)
add!(qc, CX, qr[0], qr[1])
add!(qc, CX, 1, 0)
qc.x(0)
qc.x(qr[1])

ex = "     ┌───┐┌───┐     ┌───┐┌───┐
q_0: ┤ X ├┤ Y ├──■──┤ X ├┤ X ├
     ├───┤├───┤┌─┴─┐└─┬─┘├───┤
q_1: ┤ H ├┤ Z ├┤ X ├──■──┤ X ├
     └───┘└───┘└───┘     └───┘
c: 2/═════════════════════════
                              "

@test ex == string(qc)

###############################################################
using QuantumCircuits.QCircuits.Registers

a = Qubit(1)
b = Qubit(2)
setid(a, 0)
setid(b, 0)

@test a == b

a = Qubit(1)
b = Qubit(1)
setid(a, 0)
setid(b, 1)

@test a != b

###############################################################
using QuantumCircuits.QCircuits.Circuit: getCode

qr = QuantumRegister(3, "a")
qc = QCircuit(qr)
qc.h([0, 1])
qc.x([qr[0], qr[2]])
qc.y(qr)

@test length(getCode(qc)) == 7


###############################################################

using QuantumCircuits
using QuantumCircuits.Execute

for i in 0:2^3-1   
     for j in 0:2^3-1
          nqr = QuantumRegister(3, "a", tomeasure=false)
          num_a = QuantumInteger(3)
          num_b = QuantumInteger(3)
          nqc = QCircuit([nqr, num_a, num_b])
          nqc.set!(num_a, i)
          nqc.set!(num_b, j)
          nqc.h(qr)
          nqc.measure()

          rs = getresults(nqc)
          @test abs(rs[num_a][i] - 1.0) < 1e-8
          @test abs(rs[num_b][j] - 1.0) < 1e-8
     end
end



qc.measures


# # Add 2 to register num_a
# qc.add!(num_a, 2)

###

#num_a.state

#qc.cp(pi/2, 1, 2)