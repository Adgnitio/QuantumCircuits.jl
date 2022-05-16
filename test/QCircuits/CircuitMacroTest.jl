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
using QuantumCircuits.Circuit: getCode 

using QuantumCircuits.Circuit # TODO

# Expected circuits
expc = QCircuit(2)
add!(expc, QuantumCircuits.Circuit.X, 0)
add!(expc, QuantumCircuits.Circuit.H, 1)
add!(expc, QuantumCircuits.Circuit.CX, 1, 0)

# Get the circuit from the only expression
acirc = QCircuit(2)
acirc.x(0)
acirc.h(1)
acirc.cx(1, 0)
acirc == expc


# capture the local variable
i=0
j=1
acirc = QCircuit(2)
acirc.x(i)
acirc.h(j)
acirc.cx(j, i)
@test acirc == expc
@test acirc.qubits == 2



######################################################################################
expc = QCircuit(5)
n = expc.qubits-1
add!(expc, QuantumCircuits.Circuit.H, 0:n)
foreach(i -> add!(expc, QuantumCircuits.Circuit.CX, i, i+1), 0:n-1)


# With loops
acirc = QCircuit(5)
n = acirc.qubits - 1
foreach(i -> acirc.h(i), 0:n)
foreach(i -> acirc.cx(i, i+1), 0:n-1)
@test acirc == expc


######################################################################################

# Expected circuits
expc = QCircuit(2)
add!(expc, QuantumCircuits.Circuit.Sx, 0)
add!(expc, QuantumCircuits.Circuit.Y, 1)
add!(expc, QuantumCircuits.Circuit.Z, 0)
add!(expc, QuantumCircuits.Circuit.S, 1)
add!(expc, QuantumCircuits.Circuit.Sd, 0)
add!(expc, QuantumCircuits.Circuit.T, 1)
add!(expc, QuantumCircuits.Circuit.Td, 0)
add!(expc, QuantumCircuits.Circuit.U, 1, π/5, π/6, π/7)
add!(expc, QuantumCircuits.Circuit.U3, 0, π/5, π/6, π/7)
add!(expc, QuantumCircuits.Circuit.Rx, 1, π/2)
add!(expc, QuantumCircuits.Circuit.Ry, 0, π/3)
add!(expc, QuantumCircuits.Circuit.Rz, 1, π/4)
add!(expc, QuantumCircuits.Circuit.Rzx, 0, 1, π/2)
add!(expc, QuantumCircuits.Circuit.U4, 0, 1, [i for i in 1:15])

# @gate
acirc = QCircuit(2)
acirc.sx(0)
acirc.y(1)
acirc.z(0)
acirc.s(1)
acirc.sdg(0)
acirc.t(1)
acirc.tdg(0)
acirc.u(1, π/5, π/6, π/7)
acirc.u3(0, π/5, π/6, π/7)
acirc.rx(1, π/2)
acirc.ry(0, π/3)
acirc.rz(1, π/4)
acirc.rzx(0, 1, π/2)
acirc.u4(0, 1, [i for i in 1:15])
@test acirc == expc


# Test if we can property createt with random parameters
acirc = QCircuit(2)
acirc.u3(0)
acirc.u4(0, 1)
@test length(getCode(acirc)) == 2
@test typeof(getCode(acirc)[1]) == QuantumCircuits.Circuit.U3
@test typeof(getCode(acirc)[2]) == QuantumCircuits.Circuit.U4

