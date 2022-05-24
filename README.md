Status: 
[![CI](https://github.com/Adgnitio/QuantumCircuits.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Adgnitio/QuantumCircuits.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/Adgnitio/QuantumCircuits.jl/branch/master/graph/badge.svg?token=KGJWIV6QF4)](https://codecov.io/gh/Adgnitio/QuantumCircuits.jl)

<!-- Stable version: [![Documentation](--)](--/stable/) -->
<!-- Dev version: [![Documentation](https://github.com/Adgnitio/QuantumCircuits.jl/actions/workflows/documentation.yml/badge.svg)](https://Adgnitio.github.io/QuantumCircuits.jl/dev/)-->



# QuantumCircuits.jl
QuantumCircuits is an open-source library for working with quantum computers at the application level.


## Installation
`QuantumCircuits` is in the general registry so you can install it by:
```julia
julia> import Pkg
julia> Pkg.add("QuantumCircuits")
```
Note: The library require  Qiskit installed. You can do this running the 'pip install qiskit qiskit.ignis matplotlib' command.

## Usage

```julia
using QuantumCircuits
using QuantumCircuits.Execute
```

## Introduction

We have various types of backends at our disposal, the simulator written in Julia, the Qiskit simulator, or real device available by Qiskit.

```julia
# We use the simulator written in Julia
const backend = QuantumSimulator()

# Let's create an example circuit.
qc1 = QCircuit(2)
qc1.x(0)
qc1.h(1)
qc1.cx(0, 1)
qc1
```

```
      ┌───┐     
q0_0: ┤ X ├──■──
      ├───┤┌─┴─┐
q0_1: ┤ H ├┤ X ├
      └───┘└───┘
c0: 2/══════════
```

Now, we can execute it. Because there is no measurement, we measure all qubits.

```julia
execute(backend, qc1)
```
```
4-element Vector{Float64}:
 0.0
 0.4999999999999999
 0.0
 0.5000000000000001
```

We can also add measurement explicitly.
```julia
qc1.measure(1, 1)
qc1
```
```
      ┌───┐        
q1_0: ┤ X ├──■─────
      ├───┤┌─┴─┐┌─┐
q1_1: ┤ H ├┤ X ├┤M├
      └───┘└───┘└╥┘
c1: 2/═══════════╩═
                 1 
```

```julia
execute(backend, qc1)
```
```
2-element Vector{Float64}:
 0.4999999999999999
 0.5000000000000001
```

## Registers
We can also create circuit with registers directly.

```julia
qr = QuantumRegister(3)
cr = ClassicalRegister(2)
qc = QCircuit(qr, cr)
qc.h(0)
qc.x(1)
qc.x(2)
qc.measure([0, 1], [0, 1])
qc
```
```
      ┌───┐┌─┐   
q2_0: ┤ H ├┤M├───
      ├───┤└╥┘┌─┐
q2_1: ┤ X ├─╫─┤M├
      ├───┤ ║ └╥┘
q2_2: ┤ X ├─╫──╫─
      └───┘ ║  ║ 
c2: 2/══════╩══╩═
            0  1 
```

```julia
execute(backend, qc)
```
```
4-element Vector{Float64}:
 0.0
 0.0
 0.5000000000000001
 0.4999999999999999
```
## Bug reports and Contributing
Please report any issues via the Github **[issue tracker](https://github.com/Adgnitio/QuantumCircuits.jl/issues)**. All types of issues are welcome and encouraged; this includes bug reports, documentation typos, feature requests, etc. 

QuantumCircuits is being actively developed and suggestions or other forms of contributions are encouraged. 
