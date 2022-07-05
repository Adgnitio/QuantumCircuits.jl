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

using QuantumCircuits
using QuantumCircuits.QML
using QuantumCircuits.QML.Optimization
using QuantumCircuits.Execute
using QuantumCircuits.Execute: state2probability
using QuantumCircuits.QCircuits.Circuit

using Distributions: LogNormal, quantile, cdf, pdf

μ = 0
σ = 0.1

"Cumulative normal probability distribution"
cldf(val) = cdf(LogNormal(μ, σ), val)
mydiff(v::AbstractVector) = @views v[(begin+1):end] .- v[begin:(end-1)]
pldf(val) = pdf(LogNormal(μ, σ), val)


# Qubits
N = 9    

####
npoints = 2^N
tail_value = 1e-5
start_val = quantile(LogNormal(μ, σ), tail_value)
end_val = quantile(LogNormal(μ, σ), 1 - tail_value)


# Generate expected distribution manualy
diss_x = LinRange(start_val, end_val, npoints+1)
diss_p = mydiff(cldf.(diss_x))
diss_p = diss_p ./ sum(diss_p)


# Generate ansact
qr = QuantumRegister(N)
cr = ClassicalRegister(N)
qc = QCircuit(qr, cr)
#add_ent_gate(i, j) = qc.cx(i, j)
add_ent_gate(i, j) = qc.u4(i, j)

#qc.x(4)
for i in 0:3
    add_ent_gate(4-i, 4-i-1)
    add_ent_gate(4+i, 4+i+1)
end

# decompose
qc = decompose(qc)
qc.measure(0:N-1, 0:N-1)

# Random parameter
params = getRandParameters(qc)
setparameters!(qc, params)


mse_error_loss(state_end) = sum((diss_p - state_end).^2)^0.5

function loss(params)
    start = ket"000000000"
    mat = tomatrix(qc, params)
    state_end = state2probability(mat * start)
    
    # mse error
    return mse_error_loss(state_end)
end


# check if loss function works
loss(params)
# check if derivative of loss function works
loss'(params)

# copy start params
start_params = params[:]
# run the optimization to find the best parameters
@time val, x, itr = gradientDescent(loss, loss', params, α=0.1, maxItr=30,
                              useBigValInc=true, argsArePeriodic=true)

loss(params)   
# loss start params -> 0.16722889593223778
# loss at the end, after 30 itr -> 0.025258502809051248
# 30 itrerations take 51.608174 seconds



#############################################################################
#                  Compare with Qiskit                                      #
#############################################################################
using QuantumCircuits.QCircuits.Qiskit
using QuantumCircuits.QCircuits.Qiskit: qiskit

const qiskitBackendSim = QiskitQuantum()

# number of parameters
length(params) # 99
# time of julia derivatives
@time loss'(params)
# 1.446867 seconds (61.62 k allocations: 1.857 GiB, 13.19% gc time)

@time qderivative(qiskitBackendSim, qc, mse_error_loss, params, corrMes=false)
# Job ID: 2374456f-ddda-4d4a-a8de-5f29a55d998c run 199 circuits.
# Job Status: job has successfully run
# Job 2374456f-ddda-4d4a-a8de-5f29a55d998c took 7 562 milliseconds.

