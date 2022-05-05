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

using LinearAlgebra: I
using QuantumCircuits.QCircuits.Gates: Rx, Ry, Rz
using QuantumCircuits.QCircuits.Math


################################################################################
#  Get save, parameters                                                        #
################################################################################
qc = QCircuit(1)
@gate qc.u3(0, ParameterT(1.0), ParameterT(2.0), ParameterT(3.0))
@gate qc.u3(0, ParameterT(4.0), ParameterT(5.0), ParameterT(6.0))
@test getparameters(qc) == [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]


setparameters!(qc, [1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
@test getparameters(qc) == [1.5, 2.5, 3.5, 4.5, 5.5, 6.5]

################################################################################
# Rotations
gates = [Rx, Ry, Rz]
θs = [0, π/4, π/2, π]

for g in gates
    for θ in θs
        circ = QCircuit(1)
        add!(circ, g, 0, ParameterT(θ))

        @test getparameters(circ) == [θ]
        setparameters!(circ, [1.0])
        @test getparameters(circ) == [1.0]
    end
end

################################################################################
#  Compare                                                                     #
################################################################################
qc1 = QCircuit(2)
@circ qc1 begin
    x(0)
    h(1)
    cx(0, 1)
end

qc2 = QCircuit(2)
@circ qc2 begin
    h(1)
    x(0)
    cx(0, 1)
end

@test qc1 == qc2


qr11 = QuantumRegister(1)
qr12 = QuantumRegister(1)
qc1 = QCircuit([qr11, qr12])
@gate qc1.x(qr11[0])
@gate qc1.h(qr12[0])

qc2 = QCircuit(2)
qc2.h(1)
qc2.x(0)

@test qc1 == qc2

qr11 = QuantumRegister(1)
qr12 = QuantumRegister(1)
qc1 = QCircuit([qr11, qr12])
@circ qc1 begin
    x(qr11[0])
    h(qr12[0])
end


################################################################################
function check_inv(qc)
    err = matrix_norm(tomatrix(inv(qc)) * tomatrix(qc) - I)
    @test abs(err) < 1e-16
end

# Inverse
qc = QCircuit(1)
qc.x(0)
qc.h(0)
check_inv(qc)

qc = QCircuit(1)
qc.y(0)
qc.z(0)
qc.x(0)
qc.h(0)
check_inv(qc)

qc = QCircuit(1)
qc.rx(0, π/2)
check_inv(qc)

qc = QCircuit(1)
qc.ry(0, π/2)
check_inv(qc)

qc = QCircuit(1)
qc.rz(0, π/2)
check_inv(qc)

qc = QCircuit(2)
qc.ry(0, π)
qc.cx(0, 1)
qc.rz(1, π)
qc.rx(0, π/2)
check_inv(qc)

qc = QCircuit(2)
qc.s(0)
qc.cx(1, 0)
qc.t(1)
check_inv(qc)

qc = QCircuit(2)
qc.sdg(0)
qc.cx(1, 0)
qc.tdg(1)
check_inv(qc)

qc = QCircuit(1)
qc.u3(0, π/4, π/2, π)
check_inv(qc)

qc = QCircuit(2)
qc.u3(0, π/4, π/2, π)
qc.cx(1, 0)
qc.u3(1, π/2, π/4, π/8)
check_inv(qc)


################################################################################
#  Append                                                                     #
################################################################################
qc1 = QCircuit(2)
qc1.x(0)
qc1.h(1)
qc1.cx(0, 1)
qc1.u3(1, 1, 2, 3)

qc2 = QCircuit(2)
qc2.cx(0, 1)
qc2.rx(0, 4.0)
qc2.ry(1, 5.0)
qc2.rz(1, 6.0)

testqc = QCircuit(2)
testqc.x(0)
testqc.h(1)
testqc.cx(0, 1)
testqc.u3(1, 1, 2, 3)
testqc.cx(0, 1)
testqc.rx(0, 4.0)
testqc.ry(1, 5.0)
testqc.rz(1, 6.0)

append!(qc1, qc2)

@test qc1 == testqc
