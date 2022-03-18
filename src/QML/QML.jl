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

"Quantum machine learning"
module QML


using QuantumCircuits.QCircuits.Gates: ParameterT
using QuantumCircuits.QCircuits.Circuit: getparameters, getRandParameters
using QuantumCircuits.QCircuits.QBase: setparameters!, bindparameters!

export ParameterT, getparameters, getRandParameters, setparameters!,
       bindparameters!

end # module
