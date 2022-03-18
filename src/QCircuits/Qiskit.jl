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

module Qiskit

using PyCall

using QuantumCircuits.QCircuits.QBase
using QuantumCircuits.QCircuits.Gates
using QuantumCircuits.QCircuits.OtherGates
using QuantumCircuits.QCircuits.Instructions
using QuantumCircuits.QCircuits.Registers

import QuantumCircuits.QCircuits.QBase: add!, tomatrix, measure!
import Base.show

export QiskitSimulator, QiskitCircuit, getQRegister

const qiskit = PyNULL() # pyimport("qiskit")
const plt = PyNULL()

function __init__()
    copy!(qiskit, pyimport("qiskit"))
    copy!(plt, pyimport("matplotlib.pyplot"))

    # imports
    pyimport("qiskit.ignis")
    pyimport("qiskit.ignis.verification")
end

"Quantum device provided by Qiskit."
abstract type QiskitDevice <: QuantumDevice end

"Nothing function"
const nop = () -> nothing

"Quantum circuit provided by Qiskit."
mutable struct QiskitCircuit <: QuantumCircuit
    qc::PyObject
    qubits::Dict{Qubit, PyObject}
    cbits::Dict{Cbit, PyObject}
    qRegs::Dict{String, PyObject}
    cRegs::Dict{String, PyObject}
    draw::Function

    function QiskitCircuit(qRegs::Vector{QuantumRegister},
                           cRegs::Union{Vector{ClassicalRegister}, Nothing} = nothing)
        qregs = [qiskit.QuantumRegister(length(r.bits), r.name) for r in qRegs]
        qubits = Dict{Qubit, PyObject}()
        for (qr, qqr) in zip(qRegs, qregs)
            for (k, v) in zip(qr, qqr)
                push!(qubits, k => v)
            end
        end
        qRegs = Dict{String, PyObject}()
        for r in qregs
            push!(qRegs, r.name => r)
        end

        if cRegs != nothing
            cregs = [qiskit.ClassicalRegister(length(r.bits), r.name) for r in cRegs]

            cbits = Dict{Cbit, PyObject}()
            for (cr, ccr) in zip(cRegs, cregs)
                for (k, v) in zip(cr, ccr)
                    push!(cbits, k => v)
                end
            end

            cRegs = Dict{String, PyObject}()
            for r in cregs
                push!(cRegs, r.name => r)
            end

            # , name="A"
            qc = new(qiskit.QuantumCircuit(qregs..., cregs...), qubits, cbits, qRegs, cRegs, nop)
        else
            cRegs = Dict{String, PyObject}()
            cbits = Dict{Cbit, PyObject}()
            qc = new(qiskit.QuantumCircuit(qregs...), qubits, cbits, qRegs, cRegs, nop)
        end

        qc.draw = () -> draw(qc)

        return qc
    end
end

"Draw the circuit."
function draw(qc::QiskitCircuit)
    qc.qc.draw()
    plt.show()
end

"Function return the register"
getQRegister(qc::QiskitCircuit, name::String) = qc.qRegs[name]


"Add quantum gate to Qiskit circuit"
addQiskitCode(qc::QiskitCircuit, gate::X) = qc.qc.x(qc.qubits[gate.qubit])
addQiskitCode(qc::QiskitCircuit, gate::Sx) = qc.qc.sx(qc.qubits[gate.qubit])
addQiskitCode(qc::QiskitCircuit, gate::Y) = qc.qc.y(qc.qubits[gate.qubit])
addQiskitCode(qc::QiskitCircuit, gate::Z) = qc.qc.z(qc.qubits[gate.qubit])
addQiskitCode(qc::QiskitCircuit, gate::H) = qc.qc.h(qc.qubits[gate.qubit])
addQiskitCode(qc::QiskitCircuit, gate::S) = qc.qc.s(qc.qubits[gate.qubit])
addQiskitCode(qc::QiskitCircuit, gate::Sd) = qc.qc.sdg(qc.qubits[gate.qubit])
addQiskitCode(qc::QiskitCircuit, gate::T) = qc.qc.t(qc.qubits[gate.qubit])
addQiskitCode(qc::QiskitCircuit, gate::Td) = qc.qc.tdg(qc.qubits[gate.qubit])
addQiskitCode(qc::QiskitCircuit, cx::CX) = qc.qc.cx(qc.qubits[cx.control], qc.qubits[cx.target])
addQiskitCode(qc::QiskitCircuit, gate::U3) = qc.qc.u3(getvalue(gate.θ), getvalue(gate.ϕ), getvalue(gate.λ), qc.qubits[gate.qubit])
addQiskitCode(qc::QiskitCircuit, gate::Rx) = qc.qc.rx(getvalue(gate.θ), qc.qubits[gate.qubit])
addQiskitCode(qc::QiskitCircuit, gate::Ry) = qc.qc.ry(getvalue(gate.θ), qc.qubits[gate.qubit])
addQiskitCode(qc::QiskitCircuit, gate::Rz) = qc.qc.rz(getvalue(gate.θ), qc.qubits[gate.qubit])
addQiskitCode(qc::QiskitCircuit, gate::Rzx) = qc.qc.rzx(getvalue(gate.θ), qc.qubits[gate.qubit1], qc.qubits[gate.qubit2])
addQiskitCode(qc::QiskitCircuit, gate::Barrier) = qc.qc.barrier()
function addQiskitCode(qc::QiskitCircuit, gate::U)
    qc.qc.rz(gate.β, qc.qubits[gate.qubit])
    qc.qc.ry(gate.γ, qc.qubits[gate.qubit])
    qc.qc.rz(gate.δ, qc.qubits[gate.qubit])
end

"Add gate to circuit"
function add!(qc::QiskitCircuit, gate::QuantumGate)
    addQiskitCode(qc, gate)
    nothing
end

function measure!(qc::QiskitCircuit, qubit::Qubit, cbit::Cbit)
    qc.qc.measure(qc.qubits[qubit], qc.cbits[cbit])
    nothing
end

"Qiskit simulator"
struct QiskitSimulator <: QiskitDevice
end

"Show method"
function Base.show(io::IO, qc::QiskitCircuit)
    s = string(qc.qc.draw(output="text"))
    print(io, s[10:end])
end

"Get the matrix version of the cirquit."
function tomatrix(qc::QiskitCircuit)
    backend = qiskit.BasicAer.get_backend("unitary_simulator")

    return qiskit.execute(qc.qc, backend).result().get_unitary()
end

end  # module Qiskit
