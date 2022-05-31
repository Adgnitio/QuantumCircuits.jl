# 1 Introduction
We want to implement the trotter step in an effective way, but a naive implementation of ZZ, XX, and YY gates requires 6 CX gates. We will show how to do the same using only the 3CX gate by Cartan's KAK decomposition.

# 2 ZZ, XX, and YY gates
Now we have to create ZZ, XX, and YY gates for use in simulation.

```julia
using QuantumCircuits
using QuantumCircuits.QML
using QuantumCircuits.QCircuits.Qiskit
using QuantumCircuits.QCircuits.Qiskit: qiskit
using QuantumCircuits.QCircuits.Circuit
using QuantumCircuits.QCircuits.Circuit: toQiskit, getCode
using QuantumCircuits.QCircuits.Gates: CX
using QuantumCircuits.Execute: generate_mesuere_circuits, extractProbability, correctMeasures
```

```julia
# This code com from using QuantumCircuits.Simulation.Gates module.
function ZZ(qc, q0, q1, t, usePulse=false)
    if usePulse
        qc.h(q1)
        qc.rzx(q0, q1, t)
        qc.x(q0)
        qc.rzx(q0, q1, -t)
        qc.x(q0)
        qc.h(q1)
    else
        qc.cx(q0, q1)
        qc.rz(q1, 2*t)
        qc.cx(q0, q1)
    end
end

function YY(qc, q0, q1, t, usePulse=false)
    if usePulse
        qc.sdg([q0, q1])
        qc.h(q0)
        qc.rzx(q0, q1, t)
        qc.x(q0)
        qc.rzx(q0, q1, -t)
        qc.x(q0)
        qc.h(q0)
        qc.s([q0, q1])
    else
        qc.rx([q0, q1], π/2)
        qc.cx(q0, q1)
        qc.rz(q1, 2*t)
        qc.cx(q0, q1)
        qc.rx([q0, q1], -π/2)
    end
end

function XX(qc, q0, q1, t, usePulse=false)
    if usePulse
        qc.h(q0)
        qc.rzx(q0, q1, t)
        qc.x(q0)
        qc.rzx(q0, q1, -t)
        qc.x(q0)
        qc.h(q0)
    else
        qc.ry([q0, q1], π/2)
        qc.cx(q0, q1)
        qc.rz(q1, 2*t)
        qc.cx(q0, q1)
        qc.ry([q0, q1], -π/2)
    end
end
```

# 3 U4 - Cartan's KAK decomposition
Any 2 qubits unitary can be decomposed used 3 cx gate (see the paper _"Minimal Universal Two-qubit Quantum Circuits"_ https://arxiv.org/abs/quant-ph/0308033). In QuantumCircuits library we have defined the U4 gate.

```julia
qc = QCircuit(2)
# Create u4 gate with random parameter
qc.u4(0, 1) 

# Decompose circuit to use base gates
qc2 = decompose(qc)
qc2
```

```
      ┌──────────────────────────┐┌───┐┌────────────┐                   ┌───┐»
q0_0: ┤ U3(2.2832,4.3098,1.6032) ├┤ X ├┤ Rz(3.6536) ├──■────────────────┤ X ├»
      ├─────────────────────────┬┘└─┬─┘├────────────┤┌─┴─┐┌────────────┐└─┬─┘»
q0_1: ┤ U3(4.263,1.3415,2.3539) ├───■──┤ Ry(5.2857) ├┤ X ├┤ Ry(4.1467) ├──■──»
      └─────────────────────────┘      └────────────┘└───┘└────────────┘     »
c0: 2/═══════════════════════════════════════════════════════════════════════»
                                                                             »
«      ┌────────────────────────────┐
«q0_0: ┤ U3(3.4384,0.79627,0.14553) ├
«      └┬──────────────────────────┬┘
«q0_1: ─┤ U3(4.9937,3.1561,5.6151) ├─
«       └──────────────────────────┘ 
«c0: 2/══════════════════════════════
«                                    
```

Let us assume that we would like to find parameters of the U4 gate implementing exactly ZZ, YY, and XX gate combinations for a given time $t$.

```julia
t = π/2

qc = QCircuit(2)
ZZ(qc, 0, 1, t)
YY(qc, 0, 1, t)
XX(qc, 0, 1, t)
expmat = tomatrix(qc)
qc
```

```
                         ┌─────────┐                   ┌──────────┐┌─────────┐»
q1_0: ──■─────────────■──┤ Rx(π/2) ├──■─────────────■──┤ Rx(-π/2) ├┤ Ry(π/2) ├»
      ┌─┴─┐┌───────┐┌─┴─┐├─────────┤┌─┴─┐┌───────┐┌─┴─┐├──────────┤├─────────┤»
q1_1: ┤ X ├┤ Rz(π) ├┤ X ├┤ Rx(π/2) ├┤ X ├┤ Rz(π) ├┤ X ├┤ Rx(-π/2) ├┤ Ry(π/2) ├»
      └───┘└───────┘└───┘└─────────┘└───┘└───────┘└───┘└──────────┘└─────────┘»
c1: 2/════════════════════════════════════════════════════════════════════════»
                                                                              »
«                         ┌──────────┐
«q1_0: ──■─────────────■──┤ Ry(-π/2) ├
«      ┌─┴─┐┌───────┐┌─┴─┐├──────────┤
«q1_1: ┤ X ├┤ Rz(π) ├┤ X ├┤ Ry(-π/2) ├
«      └───┘└───────┘└───┘└──────────┘
«c1: 2/═══════════════════════════════
«                                     
```

We have to choose the ansact, in our case this will be a U4 gate and in that case, we are sure that we always can find the correct parameters.

```julia
qr = QuantumRegister(2)
qc = QCircuit(qr)
qc.u4(qr[0], qr[1])

params = getRandParameters(qc)
setparameters!(qc, params)
qc = decompose(qc)
```

```
      ┌───────────────────────────┐┌───┐┌────────────┐                   ┌───┐»
q2_0: ┤ U3(1.7624,0.77116,3.9005) ├┤ X ├┤ Rz(5.2775) ├──■────────────────┤ X ├»
      └┬──────────────────────────┤└─┬─┘├────────────┤┌─┴─┐┌────────────┐└─┬─┘»
q2_1: ─┤ U3(5.0677,5.7902,6.2463) ├──■──┤ Ry(4.4635) ├┤ X ├┤ Ry(2.5877) ├──■──»
       └──────────────────────────┘     └────────────┘└───┘└────────────┘     »
«      ┌──────────────────────────┐
«q2_0: ┤ U3(0.89282,1.8969,3.728) ├
«      ├──────────────────────────┤
«q2_1: ┤ U3(4.2926,2.4936,2.6366) ├
«      └──────────────────────────┘
```

Now we can find the parameter of our ansact that perfectly fit our expected unitary matrix.


```julia
params, _, err, _  = findparam(expmat, qc, debug=false, trystandard=false)
err
```

```
8.988143676440324e-8
```

```julia
qc
```

```
      ┌──────────────────────────┐┌───┐┌──────────┐                ┌───┐»
q3_0: ┤ U3(1.4791,0.48683,4.394) ├┤ X ├┤ Rz(3π/2) ├──■─────────────┤ X ├»
      └─┬─────────────────────┬──┘└─┬─┘├──────────┤┌─┴─┐┌─────────┐└─┬─┘»
q3_1: ──┤ U3(π,5.7747,5.1031) ├─────■──┤ Ry(3π/2) ├┤ X ├┤ Ry(π/2) ├──■──»
        └─────────────────────┘        └──────────┘└───┘└─────────┘     »
«      ┌───────────────────────────┐
«q3_0: ┤ U3(-1.4791,1.8894,4.2257) ├
«      └──┬─────────────────────┬──┘
«q3_1: ───┤ U3(π,1.9652,2.8646) ├───
«         └─────────────────────┘   
```

A description of the optimization method used is available in the notebook [Eva](@ref).


```julia
```

```
```

```julia
```

```
```

