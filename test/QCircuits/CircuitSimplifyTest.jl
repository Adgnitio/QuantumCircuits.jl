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
using QuantumCircuits.QCircuits.Gates
using QuantumCircuits.QCircuits.Math
using QuantumCircuits.QCircuits.Circuit: toQiskit, getCode


################################################################################
#  Simplify - one gate                                                         #
################################################################################
const u2(qubit, ϕ, λ) = (U3, (qubit, π/2, ϕ, λ))
const u1(qubit, λ) = (U3, (qubit, 0, 0, λ))

function test_simplify(u, args, gate)
    qc = QCircuit(1)
    add!(qc, u, args...)
    sqc = simplify(qc)
    code = getCode(sqc)
    @test length(code) == 1
    @test typeof(code[1]) == gate
end
test_simplify(U3, (0, π, 0, π), X)
test_simplify(U3, (0, π, π/2, π/2), Y)
test_simplify(u1(0, π)..., Z)
test_simplify(u2(0, 0, π)..., H)
test_simplify(u1(0, π/2)..., S)
test_simplify(u1(0, -π/2)..., Sd)
test_simplify(u1(0, π/4)..., T)
test_simplify(u1(0, -π/4)..., Td)

test_simplify(U3, (0, -π, 3.1415, 0.004643), X)
# qc = QCircuit(1)
# qc.u3(0, -π, 3.1415, 0.004643)
# println(simplify(qc))
# using QGen.Quantum.Gates: Xmatrix
# unitary_error(tomatrix(qc), Xmatrix)

################################################################################
#  Simplify - identify                                                         #
################################################################################
qc = QCircuit(1)
qc.u3(0, 0, -3.678862306482877, -2.6043230006967093)
sqc = simplify(qc)
@test length(getCode(sqc)) == 0

################################################################################
#  Simplify - 2 gates                                                          #
################################################################################
function test_simplify_gates(u3, args, gates, orGate=nothing)
    qc = QCircuit(1)
    add!(qc, u3, args...)
    sqc = simplify(qc)
    code = getCode(sqc)
    @test length(code) == length(gates)
    for (i, g) in enumerate(gates)
        if isnothing(orGate)
            @test typeof(code[i]) == g
        else
            @test typeof(code[i]) == g || typeof(code[i]) == orGate[i]
        end
    end

    @test unitary_error(tomatrix(qc), tomatrix(sqc)) < 1e-8
end

test_simplify_gates(U3, (0, -π/2, π, π), [H, X])
test_simplify_gates(U3, (0, π, 0, π/2), [Y, Sd], [Sd, X])
test_simplify_gates(U3, (0, π, 0, 5π/4), [T, X])
test_simplify_gates(U3, (0, π, 2π, 0), [Z, X])
test_simplify_gates(U3, (0, 0, 2.3114387473071725, 2.400950233077517), [Sd])


# qc = QCircuit(1)
# add!(qc, U3(0, 0, 2.3114387473071725, 2.400950233077517))
# tomatrix(qc)
#
#
# sqc = simplify(qc)
#
# qc1 = QCircuit(1)
# add!(qc1, Z(0))
# add!(qc1, S(0))
# tomatrix(qc1)
# unitary_error(tomatrix(qc), tomatrix(qc1))
