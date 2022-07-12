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

export qft!

"Performs qft on the first n qubits in circuit (without swaps)."
function qft_rotations!(qc::Union{QCircuit, QuantumAbstractRegister}, code::Vector{T}, n::Integer) where T <: QuantumGate
    n -= 1
    push!(code, H(qc[n]))
    for i in (n-1):-1:0
        for j in n:-1:(i+1)
            push!(code, CP(qc[i], qc[j], π/2^(j-i)))
        end
        push!(code, H(qc[i]))
    end

    return code
end

function swap_registers!(circuit::QCircuit, reg::Union{QCircuit, QuantumAbstractRegister})
    n = length(reg)

    for qubit in 0:Int(floor(n / 2))-1
        circuit.swap(qubit, n-qubit-1)
    end

    return nothing
end

"QFT on the first n qubits in circuit"
function qft!(circuit::QCircuit, reg::Union{QCircuit, QuantumAbstractRegister} = circuit; doswap=true, inverse=false)
    doswap && inverse && swap_registers!(circuit, reg)

    code = QuantumGate[]
    qft_rotations!(reg, code, length(reg))
    
    # Inverse
    if inverse
        code = [inv(gate) for gate in reverse!(code)]
    end
    circuit.add!(code)

    doswap && !inverse && swap_registers!(circuit, reg)
    
    return nothing
end 

end  # module QLib