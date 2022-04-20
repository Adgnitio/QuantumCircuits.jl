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

module Math

using LinearAlgebra

export unitary_error, matrix_norm, eye, observe_unitary_error,
       min_observe_unitary_error

"Calculate the unitary matrix error"
unitary_error(exp_unit, mat_unit) = _unitary_error(exp_unit, mat_unit)
function unitary_error(exp_unit, mat_unit, exp_scale, mat_scale)
    return _unitary_error(exp_unit, mat_unit, x -> x * exp_scale, x -> x * mat_scale)
end
function _unitary_error(exp_unit, mat_unit, fa=identity, fb=identity)
    return sum(((a, b),) -> abs2(fa(a) - fb(b)), zip(exp_unit, mat_unit))
end

"Calculate the observe unitary matrix error, the error which we can obserwe
from quantum state"
function observe_unitary_error(exp_unit, mat_unit, index=1)
    # Remove phase
    exp_scale = cis(-angle(exp_unit[index]))
    mat_scale = cis(-angle(mat_unit[index]))

    return unitary_error(exp_unit, mat_unit, exp_scale, mat_scale)
end

"Calculate the minimum observe unitary matrix error."
function min_observe_unitary_error(exp_unit, mat_unit)
    l = length(exp_unit)
    return minimum(i -> observe_unitary_error(exp_unit, mat_unit, i), 1:l)
end

"Calculate the matrix norm"
function matrix_norm(mat)
    return abs(tr(adjoint(mat) * mat))
end

"Return the identity matrix."
eye(n) = Matrix(I, n, n)

end  # module Math
