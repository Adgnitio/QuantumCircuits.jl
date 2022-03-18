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
using QuantumCircuits.QCircuits.OtherGates: Rzx
using QuantumCircuits.QCircuits.Math
using QuantumCircuits.QCircuits.Circuit: toQiskit


################################################################################
#  2 qubit, rotation gate                                                      #
################################################################################
gates = [Rzx]
θs = [0, π/4, π/2, 3π/4, π, 5π/4, 3π/2, 7π/4, 2π]

for g in gates
    for θ in θs
        circ = QCircuit(5)
        add!(circ, g, 1, 3, θ)
        err = unitary_error(tomatrix(toQiskit(circ)), tomatrix(circ))
        if err >= 1e-8
            println(circ)
        end
        @test err < 1e-8
        test_derivate(circ, 1)
    end
end

for g in gates
    for θ in θs
        circ = QCircuit(5)
        add!(circ, g, 3, 1, θ)
        err = unitary_error(tomatrix(toQiskit(circ)), tomatrix(circ))
        if err >= 1e-8
            println(circ)
        end
        @test err < 1e-8
        test_derivate(circ, 1)
    end
end
