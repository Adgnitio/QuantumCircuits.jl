
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

using QuantumCircuits.QCircuits.Math


function test_unitary_error(exp_unit, mat_unit)
    dif = exp_unit - mat_unit
    return matrix_norm(dif)
end

A = [1+im 2+im; 3+im 4+im]
B = [1-im 2-im; 3-im 4-im]

@test abs(test_unitary_error(A, B) - unitary_error(A, B)) < 1e-6

A = rand(64, 64)
B = rand(64, 64)
@test abs(test_unitary_error(A, B) - unitary_error(A, B)) < 1e-6
