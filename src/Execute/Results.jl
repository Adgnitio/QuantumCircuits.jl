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

export Result, ResultsSet

"The result row."
struct Result
    bits::String
    p::Float64
end


"The results set."
struct ResultsSet <: AbstractDict{String, Result}
    results::Dict{String, Result}
    mapping::Vector{Int64}
    circuit::QCircuit
end

"Iterable method"
Base.iterate(rs::ResultsSet) = Base.iterate(rs.results)
Base.iterate(rs::ResultsSet, state) = Base.iterate(rs.results, state)

"Show methods"
Base.show(io::IO, r::Result) = show(io::IO, r.bits)
Base.show(io::IO, rs::ResultsSet) = show(io::IO, rs.results)
Base.show(io::IO, ::MIME{Symbol("text/plain")}, rs::ResultsSet) = show(io::IO, rs.results)

"Function execute the quantum circuit and extract the results"
function getresults(qc::QCircuit, params=nothing)
    ret = execute(qc::QCircuit, params)

    N = qc.qubits
    M = length(qc.measures)
    if M > 0
        mapping = zeros(Int64, M)
        q2c = zeros(Int64, M)
        for (q, c) in qc.measures
            mapping[getid(c)+1] = getid(q) + 1
            q2c[M - getid(q)] = N - getid(c)
        end
    else
        q2c = [i for i in 1:N]
        mapping = [i for i in 1:N]
    end

    results = Dict{String, Result}()
    for (i, v) in enumerate(ret)
        if v > 0            
            bs = bitstring(i-1)[end-qc.qubits+1:end]
            bs = String([bs[i] for i in q2c])
            push!(results, bs => Result(bs, v))
        end
    end

    return ResultsSet(results, mapping, qc)    
end

function Base.get(rs::ResultsSet, key::String, _)
    return rs.results[key]
end

function Base.get(rs::ResultsSet, qr::QuantumRegister, _)
    @assert ismeasured(rs.circuit, qr) "The register has to be measured."

    results = Dict{String, Float64}()
    for (k, v) in rs.results
        bs = String([k[rs.mapping[getid(i)+1]] for i in qr.bits])
        p = get(results, bs, 0)
        results[bs] = p + v.p
    end

    res= Dict{String, Result}()
    for (bs, p) in results
        push!(res, bs => Result(bs, p))
    end

    return res
end


end  # module Results
