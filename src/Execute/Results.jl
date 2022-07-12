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

module Results

using QuantumCircuits.QCircuits
using QuantumCircuits.QCircuits.Circuit
using QuantumCircuits.QCircuits.Registers
using QuantumCircuits.Execute.Devices

export ResultsSet

"The results set."
struct ResultsSet <: AbstractDict{String, Float64}
    results::Dict{String, Float64}
    mapping::Vector{Int64}
    circuit::QCircuit
end

"Iterable method"
Base.iterate(rs::ResultsSet) = Base.iterate(rs.results)
Base.iterate(rs::ResultsSet, state) = Base.iterate(rs.results, state)

"Show methods"
Base.show(io::IO, rs::ResultsSet) = show(io::IO, rs.results)
Base.show(io::IO, ::MIME{Symbol("text/plain")}, rs::ResultsSet) = show(io::IO, rs.results)

"Function execute the quantum circuit and extract the results"
function getresults(qc::QCircuit, params=nothing)
    ret = execute(qc::QCircuit, params)

    N = qc.qubits
    M = length(qc.measures)
    if M > 0
        mapping = zeros(Int64, N)
        q2c = zeros(Int64, N)        
        for (q, c) in qc.measures        
            q2c[getid(c) + 1] = getid(q) + 1
            mapping[getid(q) + 1] = getid(c) + 1
        end
    else
        q2c = [i for i in 1:N]
        mapping = [i for i in 1:M]
    end
    # remove all zeros
    deleteat!(q2c, findall(x->x==0, q2c))

    measured_qubits = Set([getid(q) + 1 for (q, c) in qc.measures])
    results = Dict{String, Float64}()
    # if there is no measurmend, we measure all
    if M == 0
        nM = N
    else
        nM = M
    end
    nor_sum = 0.0
    for (i, v) in enumerate(ret)
        if v > 1e-10                        
            bs = reverse(bitstring(i-1)[end-nM+1:end])
            if N > M && M != 0
                next = 1
                nbs = ""
                for j in 1:N                    
                    if j in measured_qubits
                        nbs = nbs * bs[next]
                        next += 1
                    else
                        nbs = nbs * "X"
                    end
                end
                bs = nbs
            end

            bs = reverse(String([bs[i] for i in q2c]))
            push!(results, bs => v)

            nor_sum += v
        end
    end

    # Normalize
    for (v, p) in results
        results[v] = p/nor_sum
    end

    return ResultsSet(results, mapping, qc)  
end

function Base.get(rs::ResultsSet, key::String, _)
    return rs.results[key]
end

function innerget(rs::ResultsSet, qr::QuantumAbstractRegister)
    @assert ismeasured(rs.circuit, qr) "The register has to be measured."

    results = Dict{String, Float64}()
    for (k, v) in rs.results
        bs = reverse(String([reverse(k)[rs.mapping[getid(i)+1]] for i in qr.bits]))
        p = get(results, bs, 0)
        results[bs] = p + v
    end

    return results
end

Base.get(rs::ResultsSet, qr::QuantumRegister, _) = innerget(rs, qr)

function Base.get(rs::ResultsSet, qr::QuantumInteger, _)    
    results = innerget(rs, qr)

    res = Dict{Int64, Float64}()
    for (k, p) in results
        push!(res, parse(Int64, k; base=2) => p)
    end

    return res
end

function Base.get(rs::ResultsSet, qr::QuantumFloat, _)    
    results = innerget(rs, qr)

    res = Dict{Float64, Float64}()
    for (k, p) in results
        push!(res, parse(Int64, k; base=2) / 2^qr.fractional => p)
    end

    return res
end


end  # module Results
