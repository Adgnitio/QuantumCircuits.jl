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
using QuantumCircuits.QCircuits.Registers
using QuantumCircuits.QCircuits.Circuit
using QuantumCircuits.QCircuits.Gates: H, CP


"Performs qft on the first n qubits in circuit (without swaps)."
function qft_rotations!(qc::Union{QCircuit, QuantumRegister}, code::Vector{T}, n::Integer; inverse=false) where T <: QuantumGate
    # Break the recurence
    n == 0 && return nothing

    n -= 1
    #circuit.h(n)
    push!(code, H(qc[n]))
    for qubit in 0:(n-1)
        #circuit.cp(qubit, n, π/2^(n-qubit))
        push!(code, CP(qc[qubit], qc[n], π/2^(n-qubit)))
    end

    # At the end of our function, we call the same function again on
    # the next qubits (we reduced n by one earlier in the function)
    qft_rotations!(qc, code, n)

    if inverse
        return reverse!(code)
    else
        return code
    end
end

function swap_registers!(circuit::QCircuit, n::Integer)
    for qubit in 0:Int(floor(n / 2))-1
        circuit.swap(qubit, n-qubit-1)
    end

    return nothing
end

"QFT on the first n qubits in circuit"
function qft!(circuit::QCircuit; doswap=true, inverse=false)
    doswap && inverse && swap_registers!(circuit, circuit.qubits)

    code = QuantumGate[]
    qft_rotations!(circuit, code, circuit.qubits, inverse=inverse)
    circuit.add!(code)

    doswap && !inverse && swap_registers!(circuit, circuit.qubits)
    
    return nothing
end 

end  # module QLib