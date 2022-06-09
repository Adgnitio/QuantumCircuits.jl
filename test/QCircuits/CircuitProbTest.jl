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
using QuantumCircuits.QML

using QuantumCircuits.QCircuits.Math
using QuantumCircuits.QCircuits.Circuit
using QuantumCircuits.QCircuits.Circuit: toQiskit
using QuantumCircuits.QCircuits.Qiskit
using QuantumCircuits.QCircuits.Qiskit: qiskit
using QuantumCircuits.Execute

statevectorsim = qiskit.Aer.get_backend("statevector_simulator")

function calcQiskitStateVector(qc)
    qcirc = toQiskit(qc)

    statevector = qiskit.execute(qcirc.qc, statevectorsim).result().get_statevector()
    return statevector.probabilities()    
end

################################################################################
#  U3                                                                          #
################################################################################
circ_N = 1
circ_qr = QuantumRegister(circ_N)
circ_cr = ClassicalRegister(circ_N)
circ = QCircuit(circ_qr, circ_cr)
circ.u3(0, 1.3537, 2.5652, 5.5142)

@test sum(abs.(execute(circ) - calcQiskitStateVector(circ))) < 1e-6

