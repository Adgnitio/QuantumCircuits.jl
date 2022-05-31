using Test

using QuantumCircuits
using QuantumCircuits.QCircuits.Math

using QuantumCircuits.Simulation.Gates


gates = [ZZ, XX, YY]
θs = [0, π/4, π/2, 3π/4, π, 5π/4, 3π/2, 7π/4, 2π]

for g in gates
    for θ in θs
        simqc = QCircuit(2)
        g(simqc, 0, 1, θ, true)
        m1 = tomatrix(simqc)

        simqc = QCircuit(2)
        g(simqc, 0, 1, θ, false)
        m2 = tomatrix(simqc)

        @test unitary_error(m1, m2) < 1e-8
    end
end
