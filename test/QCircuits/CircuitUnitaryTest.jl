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

using QuantumCircuits.QCircuits.Math
using QuantumCircuits.QCircuits.Circuit: toQiskit
using QuantumCircuits.QCircuits.Gates: Rx, Ry, Rz, X, Y, Z, S, Sd, T, Td, H, CX, U3, P, CP, Swap

import Zygote

function test_derivate(qc, explen)
    pa = getRandParameters(qc)
    @test length(pa) == explen
    loss(params) = real(matrix_norm(tomatrix(qc, params)))
    loss(pa)
    Zygote.gradient(loss, pa)[1]
end

################################################################################
#  1 qubit, rotation gate                                                      #
################################################################################
gates = [Rx, Ry, Rz, P]
θs = [0, π/4, π/2, 3π/4, π, 5π/4, 3π/2, 7π/4, 2π]

for g in gates
    for θ in θs
        circ = QCircuit(1)
        add!(circ, g, 0, ParameterT(θ))
        err = unitary_error(tomatrix(toQiskit(circ)), tomatrix(circ))
        if err >= 1e-8
            println(circ)
        end
        @test err < 1e-8
        test_derivate(circ, 1)
    end
end


################################################################################
#  1 qubit, 1 gate unitary                                                     #
################################################################################
gates = [X, Y, Z, S, Sd, T, Td, H]

for g in gates
    circ = QCircuit(1)
    add!(circ, g, 0)
    @test unitary_error(tomatrix(toQiskit(circ)), tomatrix(circ)) < 1e-8
    test_derivate(circ, 0)
end

################################################################################
#  1 qubit, 2/3 gate unitary                                                   #
################################################################################
for g in gates
    for g2 in gates
        circ = QCircuit(1)
        add!(circ, g, 0)
        add!(circ, g2, 0)
        @test unitary_error(tomatrix(toQiskit(circ)), tomatrix(circ)) < 1e-8
        test_derivate(circ, 0)
    end
end
for g in gates
    for g2 in gates
        for g3 in gates
            circ = QCircuit(1)
            add!(circ, g, 0)
            add!(circ, g2, 0)
            add!(circ, g3, 0)
            @test unitary_error(tomatrix(toQiskit(circ)), tomatrix(circ)) < 1e-8
            test_derivate(circ, 0)
        end
    end
end

qc = QCircuit(1)
qc.u3(0, ParameterT(0), ParameterT(1.0), ParameterT(2.0))
qc.x(0)
qc.u3(0, ParameterT(2.0), ParameterT(3.0), ParameterT(4.0))
qiskitmatrix = tomatrix(toQiskit(qc))
@test unitary_error(qiskitmatrix, tomatrix(qc)) < 1e-8
params = [4.0, 3.0, 2.0, 1.0, 0.0, -1.0]
gcparammatrix = tomatrix(qc, params)
@test unitary_error(qiskitmatrix, gcparammatrix) > 1.0
setparameters!(qc, params)
qiskitmatrix = tomatrix(toQiskit(qc))
@test unitary_error(qiskitmatrix, gcparammatrix) < 1e-8

################################################################################
#  2 qubit, 2 gate unitary                                                     #
################################################################################
for g in gates
    for g2 in gates
        circ = QCircuit(2)
        add!(circ, g, 0)
        add!(circ, g2, 1)
        @test unitary_error(tomatrix(toQiskit(circ)), tomatrix(circ)) < 1e-8
        test_derivate(circ, 0)
    end
end

for g in gates
    for g2 in gates
        for g3 in gates
            circ = QCircuit(2)
            add!(circ, g, 0)
            add!(circ, g2, 1)
            add!(circ, g3, 1)
            @test unitary_error(tomatrix(toQiskit(circ)), tomatrix(circ)) < 1e-8
            test_derivate(circ, 0)
        end
    end
end

for g in gates
    for g2 in gates
        for g3 in gates
            circ = QCircuit(3)
            add!(circ, g, 0)
            add!(circ, g2, 1)
            add!(circ, g3, 2)
            @test unitary_error(tomatrix(toQiskit(circ)), tomatrix(circ)) < 1e-8
            test_derivate(circ, 0)
        end
    end
end

################################################################################
#  CX                                                                          #
################################################################################
function test_unitary(qubits, gates, expected_param=0; inline_optimization=true, usedecompose=false)
    circ = QCircuit(qubits, inline_optimization=inline_optimization)
    for (g, args) in gates
        add!(circ, g, args...)
    end
    if usedecompose
        dcirc = decompose(circ)
    else
        dcirc = circ
    end


    @test unitary_error(tomatrix(toQiskit(circ)), tomatrix(dcirc)) < 1e-8
    test_derivate(dcirc, expected_param)
end

test_unitary(2, [(CX, [0, 1])])
test_unitary(2, [(CX, [1, 0])])
test_unitary(2, [(H, [0]), (H, [1]), (CX, [0, 1]), (X, [0]), (Z, [1])])
test_unitary(2, [(H, [0]), (H, [1]), (CX, [1, 0]), (X, [0]), (Z, [1])])
test_unitary(6, [(H, [0]), (H, [1]), (Z, [2]), (CX, [2, 3]), (S, [0]), (Y, [1]),
                 (Td, [5])])
test_unitary(6, [(H, [0]), (H, [1]), (CX, [3, 2]), (X, [0]), (Sd, [1]), (H, [4])])

test_unitary(3, [(CX, [1, 2])])

test_unitary(4, [(CX, [0, 3])])
test_unitary(4, [(CX, [0, 3]), (H, [0]), (Y, [1])])
test_unitary(4, [(CX, [3, 0])])
test_unitary(5, [(CX, [1, 3])])
test_unitary(5, [(CX, [3, 1])])

test_unitary(6, [(H, [0]), (H, [1]), (CX, [2, 3]), (S, [0]), (Y, [1]),
                 (Td, [5]), (Z, [2]), (CX, [0, 3]), (H, [1]),
                 (Z, [2]), (S, [5]), (CX, [0, 3]), (Td, [3]), (H, [4]),
                 (Y, [1]), (Y, [2]), (CX, [1, 3]), (H, [5]), (H, [3]),
                 (CX, [3, 1]), (Z, [2])])


################################################################################
#  U3                                                                          #
################################################################################
test_unitary(1, [(U3, (0, ParameterT(π/4), ParameterT(π/2), ParameterT(π)))], 3)
test_unitary(1, [(U3, (0, ParameterT(-π/4), ParameterT(-π/2), ParameterT(-π)))], 3)
test_unitary(1, [(U3, (0, ParameterT(π/4), ParameterT(π/2), ParameterT(π))), (U3, (0, ParameterT(-π/4), ParameterT(-π/2), ParameterT(-π)))], 6, inline_optimization=false)
test_unitary(2, [(U3, (0, ParameterT(π/4), ParameterT(π/2), ParameterT(π))), (CX, [1, 0]), (U3, (1, ParameterT(-π/4), ParameterT(-π/2), ParameterT(-π)))], 6)

################################################################################
#  CP                                                                          #
################################################################################
λs = [0, π/4, π/2, 3π/4, π, 5π/4, 3π/2, 7π/4, 2π]

for λ in λs
    test_unitary(2, [(X, [0]), (X, [1]), (CP, (0, 1, λ))], usedecompose=true)
    test_unitary(2, [(CP, (1, 0, λ))], usedecompose=true)
    test_unitary(2, [(X, [0]), (X, [1]), (CP, (1, 0, λ))], usedecompose=true)
    test_unitary(2, [(H, [1]), (CP, (1, 0, λ), (H, [0]))], usedecompose=true)

    # minus
    test_unitary(2, [(X, [0]), (X, [1]), (CP, (0, 1, -λ))], usedecompose=true)
    test_unitary(2, [(CP, (1, 0, -λ))], usedecompose=true)
    test_unitary(2, [(X, [0]), (X, [1]), (CP, (1, 0, -λ))], usedecompose=true)
    test_unitary(2, [(H, [1]), (CP, (1, 0, -λ), (H, [0]))], usedecompose=true)
end


################################################################################
#  Swap                                                                        #
################################################################################
test_unitary(2, [(X, [0]), (X, [1]), (Swap, [0, 1])], usedecompose=true)
test_unitary(2, [(Swap, [0, 1])], usedecompose=true)
test_unitary(2, [(X, [0]), (Swap, [0, 1])], usedecompose=true)
