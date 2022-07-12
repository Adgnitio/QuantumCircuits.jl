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

changebase!(qc::QCircuit, reg::T, state::Reg.QuantumNumberState) where {T <: QuantumAbstractRegister} = nothing

function add!(qc::QCircuit, reg::QuantumInteger, num::Number)
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


end  # module Arithmetic
