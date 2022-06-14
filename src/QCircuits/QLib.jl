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

module QLib

using QuantumCircuits.QCircuits.QBase
using QuantumCircuits.QCircuits.Circuit

"Performs qft on the first n qubits in circuit (without swaps)."
function qft_rotations!(circuit::QCircuit, n::Integer)
    # Break the recurence
    n == 0 && return nothing

    n -= 1
    circuit.h(n)
    for qubit in 0:(n-1)
        circuit.cp(qubit, n, π/2^(n-qubit))
    end

    # At the end of our function, we call the same function again on
    # the next qubits (we reduced n by one earlier in the function)
    qft_rotations!(circuit, n)
end

# "Performs inverse qft on the first n qubits in circuit (without swaps)."
# function inv_qft_rotations!(circuit::QCircuit, n::Integer)
#     # Break the recurence
#     n == 0 && return nothing

#     n -= 1
#     circuit.h(n)
#     for qubit in 0:(n-1)
#         circuit.cp(qubit, n, π/2^(n-qubit))
#     end

#     # At the end of our function, we call the same function again on
#     # the next qubits (we reduced n by one earlier in the function)
#     qft_rotations!(circuit, n)
# end

function swap_registers!(circuit::QCircuit, n::Integer)
    for qubit in 0:Int(floor(n / 2))-1
        @show qubit, n-qubit-1
        circuit.swap(qubit, n-qubit-1)
    end

    return nothing
end

"QFT on the first n qubits in circuit"
function qft!(circuit::QCircuit, n::Integer; doswap=true)
    qft_rotations!(circuit, n)
    doswap && swap_registers!(circuit, n)
    
    return nothing
end 

# "Does the inverse QFT on the first n qubits in circuit"
# function inverse_qft!(circuit::QCircuit, n::Integer; doswap=true)
#     # First we create a QFT circuit of the correct size:
#     qft_circ = qft(QCircuit(n), n)

#     # Then we take the inverse of this circuit
#     invqft_circ = qft_circ.inverse()
#     # And add it to the first n qubits in our existing circuit
#     circuit.append(invqft_circ, circuit.qubits[:n])
#     #return circuit.decompose() # .decompose() allows us to see the individual gates

#     return nothing
# end 

end  # module QLib