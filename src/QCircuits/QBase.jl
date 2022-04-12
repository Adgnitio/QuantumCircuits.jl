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

module QBase

export QuantumCircuit, QuantumDevice, QuantumGate, QuantumObject,
       add!, tomatrix, setparameters!, simplify, measure!,
       standardGateError, bindparameters!, decompose

"Quantum object :)"
abstract type QuantumObject end

"Quantum Gate"
abstract type QuantumGate <: QuantumObject end

"Quantum Circuit object"
abstract type QuantumCircuit <: QuantumObject end

"Quantum davice, it may be simulator or real hardware."
abstract type QuantumDevice end


"The quantum gates"
add!(qc::QuantumCircuit, ::QuantumGate) = error("Please implement add! function used to adding gate to circuit for ", typeof(qc))

"Get the matrix version of the quantum object."
tomatrix(qc::QuantumObject) = error("Please implement tomatrix function for ", typeof(qc))
tomatrix(qc::QuantumObject, ::Any) = error("Please implement tomatrix with params function for ", typeof(qc))

"Method return the error of standard gate"
standardGateError(qc::QuantumObject) = error("Please implement standardGateError function for ", typeof(qc))
standardGateError(qc::QuantumObject, ::Any) = error("Please implement standardGateError with params function for ", typeof(qc))

"Set the parameters to the quantum object"
setparameters!(qc::QuantumObject, ::Any) = error("Please implement setparameters! function for ", typeof(qc))

"Bind the parameters to the quantum object"
bindparameters!(qc::QuantumObject) = error("Please implement bindparameters! function for ", typeof(qc))

"Simplify the quantum object"
simplify(qc::QuantumObject) = error("Please implement simplify function for ", typeof(qc))

decompose(qc::QuantumObject) = error("Please implement decompose function for ", typeof(qc))

Base.inv(qc::QuantumObject) = error("Please implement inv function for ", typeof(qc))

measure!(qc::QuantumObject) = error("Please implement measure! function for ", typeof(qc))

end  # module QBase
