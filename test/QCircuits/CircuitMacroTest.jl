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
@gate acirc.x(0)
@gate acirc.h(1)
@gate acirc.cx(1, 0)
@test acirc == expc


# capture the local variable
i=0
j=1
acirc = QCircuit(2)
@gate acirc.x(i)
@gate acirc.h(j)
@gate acirc.cx(j, i)
@test acirc == expc

# This is not transformed because it doesn't match a function call
# So it works as intended.
@gate acirc.qubits
@test acirc.qubits == 2


######################################################################################
# qasm-like macro
acirc = QCircuit(2)
@circ acirc begin
    x(0)
    h(1)
    cx(1, 0)
end
@test acirc == expc

# qasm-like macro with variables
i = -1
j = -1
acirc = QCircuit(2)
@circ acirc begin
    i = 0
    j = 1

    x(i)
    h(j)
    cx(j, i)
end
@test acirc == expc

# qasm-like with function
acirc = QCircuit(2)
@circ acirc begin
    i = 0
    j = i + 1

    x(i)
    h(j)
    cx(j, i)
end
@test acirc == expc



######################################################################################
expc = QCircuit(5)
n = expc.qubits-1
add!(expc, QuantumCircuits.Circuit.H, 0:n)
foreach(i -> add!(expc, QuantumCircuits.Circuit.CX, i, i+1), 0:n-1)


# With loops
acirc = QCircuit(5)
@circ acirc begin
    #TODO this doesn't works
    n = 5-1 #acirc.qubits 
    foreach(i -> h(i), 0:n)
    foreach(i -> cx(i, i+1), 0:n-1)
end
@test acirc == expc


# The same with gates
acirc = QCircuit(5)
n = acirc.qubits-1
@gate acirc.h(0:n)
for i in 0:n-1
    @gate acirc.cx(i, i+1)
end
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
@gate acirc.sx(0)
@gate acirc.y(1)
@gate acirc.z(0)
@gate acirc.s(1)
@gate acirc.sdg(0)
@gate acirc.t(1)
@gate acirc.tdg(0)
@gate acirc.u(1, π/5, π/6, π/7)
@gate acirc.u3(0, π/5, π/6, π/7)
@gate acirc.rx(1, π/2)
@gate acirc.ry(0, π/3)
@gate acirc.rz(1, π/4)
@gate acirc.rzx(0, 1, π/2)
@gate acirc.u4(0, 1, [i for i in 1:15])
@test acirc == expc

# @circ
acirc = QCircuit(2)
@circ acirc begin
    sx(0)
    y(1)
    z(0)
    s(1)
    sdg(0)
    t(1)
    tdg(0)
    u(1, π/5, π/6, π/7)
    u3(0, π/5, π/6, π/7)
    rx(1, π/2)
    ry(0, π/3)
    rz(1, π/4)
    rzx(0, 1, π/2)
    u4(0, 1, [i for i in 1:15])
end
@test acirc == expc



# Test if we can property createt with random parameters
acirc = QCircuit(2)
@circ acirc begin
    u3(0)
    u4(0, 1)
end
@test length(getCode(acirc)) == 2
@test typeof(getCode(acirc)[1]) == QuantumCircuits.Circuit.U3
@test typeof(getCode(acirc)[2]) == QuantumCircuits.Circuit.U4

