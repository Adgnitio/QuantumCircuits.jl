using Test


using QuantumCircuits
using QuantumCircuits.QCircuits.Circuit


qc1 = QCircuit(2)
qc1.x(0)
qc1.h(1)
qc1.cx(0, 1)
@test toPythonQiskitCode(qc1) == "qr1 = QuantumRegister(2)\ncr1 = ClassicalRegister(2)\nqc = QuantumCircuit(qr1, cr1)\nqc.x(0)\nqc.h(1)\nqc.cx(0, 1)\n"


qr1 = QuantumRegister(3, "q")
cr1 = ClassicalRegister(2, "c")
qc1 = QCircuit(qr1, cr1)
qc1.x(0)
qc1.h(1)
qc1.cx(0, 1)
qc1.measure([0, 1], [0, 1])
@test toPythonQiskitCode(qc1) == "qr1 = QuantumRegister(3, \"q\")\ncr1 = ClassicalRegister(2, \"c\")\nqc = QuantumCircuit(qr1, cr1)\nqc.x(0)\nqc.h(1)\nqc.cx(0, 1)\nqc.measure(0, 0)\nqc.measure(1, 1)\n"
