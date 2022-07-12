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

module Arithmetic

using QuantumCircuits.QCircuits
using QuantumCircuits.QCircuits.Circuit
using QuantumCircuits.QCircuits.Registers
using QuantumCircuits.QCircuits.QLib

import QuantumCircuits.QCircuits.Registers as Reg
import QuantumCircuits.QCircuits.QBase: qft!, add!, changebase!

###############################################################
#  Register - Change base                                     #
###############################################################
function changebase!(qc::QCircuit, reg::T, state::Reg.QuantumNumberState) where {T <: QuantumNumber}
    if state == Reg.QFTbase && (reg.state == Reg.NormalBase || reg.state == Reg.Empty)
        qft!(qc, reg, doswap=false)

        reg.state = Reg.QFTbase
        return nothing
    elseif state == Reg.NormalBase && reg.state == Reg.QFTbase 
        qft!(qc, reg, doswap=false, inverse=true)

        reg.state = Reg.NormalBase
        return nothing
    end
end

changebase!(::QCircuit, ::T, ::Reg.QuantumNumberState) where {T <: QuantumAbstractRegister} = nothing

function innerAdd!(qc::QCircuit, reg::QuantumNumber, num::Number)
    # QFT
    changebase!(qc, reg, Reg.QFTbase)

    # Addition
    for i in 0:length(reg)-1
        qc.p(reg[i], num * 2π/2^(i+1))
    end

    # Inverse QFT
    # we, don't do the inverse. We try to do this in leasy maner.
    # qft!(qc, reg, doswap=false, inverse=true)
end

function add!(qc::QCircuit, reg::QuantumInteger, num::Number)
    @assert num < 2^reg.integer "The number $num is out of the register with $(reg.integer) qubits."

    innerAdd!(qc, reg, num)
end

function add!(qc::QCircuit, reg::QuantumFloat, num::Number)
    val = getIntValue(reg, num)

    innerAdd!(qc, reg, val)
end

end  # module Arithmetic
