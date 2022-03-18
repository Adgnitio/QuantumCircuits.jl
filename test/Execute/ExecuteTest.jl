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
using QuantumCircuits.Execute
using QuantumCircuits.Execute: _0,  _1, qjacobian
using QuantumCircuits.QCircuits.Math
using QuantumCircuits.QCircuits.Circuit: toQiskit
using QuantumCircuits.QCircuits.Circuit
#using QuantumCircuits.QCircuits.Registers

using Zygote
using Zygote: jacobian

# using QGen.Quantum.Execute
# using QGen.Quantum.Qiskit
# using QGen.Quantum.Gates
# using QGen.Quantum.Circuit
# using QGen.Quantum.Circuit: toQiskit
# using QGen.Quantum.Math
# using QGen.Quantum.Qiskit
# using QGen.Quantum.Registers



@test ket"0" == _0
@test ket"1" == _1

@test ket"00" == kron(_0, _0)
@test ket"01" == kron(_0, _1)
@test ket"10" == kron(_1, _0)
@test ket"11" == kron(_1, _1)

@test ket"0101" == kron(_0, kron(_1, kron(_0, _1)))

const backend = QuantumSimulator()
const qiskitBackend = QiskitQuantum()

#############################################################################

err(x, y) = sum(abs.(x - y))

function doTests(qc, expected, measuere_mat, checkQiskit=true)
    @test err(execute(backend, qc), expected) < 1e-5
    @test qc.measures_matrix == measuere_mat
    if checkQiskit
        @test err(execute(qiskitBackend, toQiskit(qc)), expected) < 0.05
    end
end

circuits = QCircuit[]

qc = QCircuit(2)
qc.h(0)
qc.x(1)
qc.measure([0, 1], [0, 1])
doTests(qc, [0.0, 0.0, 0.5, 0.5], eye(4))
push!(circuits, qc)

###

qc = QCircuit(2)
qc.h(0)
qc.x(1)
qc.measure(0, 0)
doTests(qc, [0.5, 0.5], [1 0 1 0; 0 1 0 1], false)

###

qr = QuantumRegister(2)
cr = ClassicalRegister(1)
qc = QCircuit(qr, cr)
qc.h(0)
qc.x(1)
qc.measure(0, 0)
doTests(qc, [0.5, 0.5], [1 0 1 0; 0 1 0 1])
push!(circuits, qc)

###

qc = QCircuit(2)
qc.h(0)
qc.x(1)
qc.measure(1, 0)
doTests(qc, [0.0, 1.0], [1 1 0 0; 0 0 1 1], false)

###

qc = QCircuit(2)
qc.h(0)
qc.x(1)
qc.measure(1, 1)
doTests(qc, [0.0, 1.0], [1 1 0 0; 0 0 1 1], false)

###

qr = QuantumRegister(2)
cr = ClassicalRegister(1)
qc = QCircuit(qr, cr)
qc.h(0)
qc.x(1)
qc.measure(1, 0)
doTests(qc, [0.0, 1.0], [1 1 0 0; 0 0 1 1])
push!(circuits, qc)

###

qr = QuantumRegister(3)
cr = ClassicalRegister(2)
qc = QCircuit(qr, cr)
qc.h(0)
qc.x(1)
qc.measure([0, 2], [0, 1])
expected = [0.5, 0.5, 0.0, 0.0]
@test err(execute(backend, qc), expected) < 1e-5
@test err(execute(qiskitBackend, toQiskit(qc)), expected) < 0.05
push!(circuits, qc)

###

qr = QuantumRegister(3)
cr = ClassicalRegister(2)
qc = QCircuit(qr, cr)
qc.h(0)
qc.x(1)
qc.measure([0, 1], [0, 1])
expected = [0.0, 0.0, 0.5, 0.5]
@test err(execute(backend, qc), expected) < 1e-5
@test err(execute(qiskitBackend, toQiskit(qc)), expected) < 0.05
push!(circuits, qc)

#############################################################################
res_sim = execute(backend, circuits)
res_dev = execute(qiskitBackend, [toQiskit(qc) for qc in circuits], debug=false)

@test sum([err(s, d) for (s, d) in zip(res_sim, res_dev)]) < 0.5

#############################################################################

qr = QuantumRegister(3)
cr = ClassicalRegister(2)
qc = QCircuit(qr, cr)
qc.rx(0)
qc.rx(1)
qc.rx(2)
qc.measure([0, 1], [0, 1])


pa = [1.0, 1.0, 1.0]
execute(backend, qc, pa)

loss(params::Vector{Float64}) = loss_expected_zero_state(execute(backend, qc, params))
loss(pa)
loss'(pa)


#############################################################################


@test loss_expected_zero_state(ket"0") < 1e-8
@test loss_expected_zero_state(ket"00") < 1e-8
@test loss_expected_zero_state(ket"000") < 1e-8
@test loss_expected_zero_state(ket"0000") < 1e-8

val1 = loss_expected_zero_state([1/2 1/6 1/6 1/6])
val2 = loss_expected_zero_state([1/2 1/2 0.0 0.0])
@test val1 + 0.1 < val2

val1 = loss_expected_zero_state([0.0 1/3 1/3 1/3])
val2 = loss_expected_zero_state([0.0 1.0 0.0 0.0])
@test val1 + 0.1 < val2

#############################################################################

qr = QuantumRegister(3)
cr = ClassicalRegister(3)
qc = QCircuit(qr, cr)
qc.rx(0)
qc.rx(1)
qc.rx(2)
qc.measure([0, 1, 2], [0, 1, 2])

pa = [1.0, 0.2, 0.6]

# Simulator
executeSim(params::Vector{Float64}) = execute(backend, qc, params)
lossAD(params::Vector{Float64}) = loss_expected_zero_state(executeSim(params))

res_sim = executeSim(pa)
res_sim_los = lossAD(pa)
res_sim_der = real.(lossAD'(pa))
res_sim_jak = jacobian(executeSim, pa)[1]

res_dev, res_dev_jak = qjacobian(qiskitBackend, qc, pa)
@test unitary_error(res_dev, res_sim) < 0.001
@test unitary_error(res_dev_jak, res_sim_jak) < 0.001


fn(params::Vector) = qderivative(qiskitBackend, qc, loss_expected_zero_state, params)
res_dev_los, res_dev_der = fn(pa)

@test unitary_error(res_dev_los, res_sim_los) < 0.001
@test unitary_error(res_dev_der, res_sim_der) < 0.001

#############################################################################

qr = QuantumRegister(2)
qc = QCircuit(qr)
qc.h(0)
qc.x(1)
@test err(execute(backend, qc), [0.0, 0.0, 0.5, 0.5]) < 1e-5
@test qc.measures_matrix == eye(4)

# Add measures
cr = ClassicalRegister(1)
setClassicalRegister!(qc, cr)
qc.measure([0], [0])
@test err(execute(backend, qc), [0.5, 0.5]) < 1e-5
@test qc.measures_matrix == [1 0 1 0; 0 1 0 1]

#############################################################################
n = 3
qr = QuantumRegister(7, "q")
cr = ClassicalRegister(3)
qc = QCircuit(qr, cr)
qr = [qr[1], qr[3], qr[5]]

qc.u3(qr)
for i in (n-2):-1:0
    i = i+1
    qc.cx(qr[i], qr[i+1])
    #qc.rzx(qr[i], qr[i+1])
end
qc.u3(qr)
for i in 0:(n-2)
    i = i+1
    qc.cx(qr[i], qr[i+1])
    #qc.rzx(qr[i], qr[i+1])
end
qc.u3(qr)
qc.measure([1, 3, 5], [0, 1, 2])


pa = [3.762539403251488, 0.001127890377849541, 2.1658488014045303, 6.174052504545824, 4.9629009414952305, 0.38255621799191375, 1.7274354859392205, 1.0857210633046896, 3.103638820041333, 3.758391851901853, 0.5765051731246599, 0.5798377641373272, 3.8759799027809825, 0.6188604600285638, 4.2218535852028385, 2.2390389868469462, 2.6674508252859623, 5.878739321230936, 1.340082725285039, 5.916066560675415, 3.625237411363865, 5.4397023221566725, 6.212355261857649, 2.5222234263059673, 5.8952667846227405, 6.113119466717698, 0.9810581564033299]

loss(params) = loss_expected_zero_state(execute(backend, qc, params))
dloss(params) = real(loss'(params))

ml = loss(pa)
mdl = dloss(pa)

sl, sdl = qderivative(qiskitBackend, qc, loss_expected_zero_state, pa, 4)
#println(join(pa, ", "))
@test abs(ml - sl) < 0.2
@test maximum(abs.(mdl - sdl)) < 2.0

# sl2, sdl2 = qderivative(qiskitBackend, ansact, loss_expected_zero_state, pa, 1)
# @test abs(ml - sl) - abs(ml - sl2) < 0.05
# @test maximum(abs.(mdl - sdl)) - maximum(abs.(mdl - sdl2)) < 0.1
