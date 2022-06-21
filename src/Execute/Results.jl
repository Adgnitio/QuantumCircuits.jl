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

    results = Dict{String, Result}()
    for (i, v) in enumerate(ret)
        if v > 0
            bs = bitstring(i-1)[end-qc.qubits+1:end]
            push!(results, bs => Result(bs, v))
        end
    end

    return ResultsSet(results, qc)
end

function Base.get(rs::ResultsSet, key::String, _)
    return rs.results[key]
end

# function Base.get(rs::ResultsSet, qr::QuantumRegister, _)

#     return rs.results[key]
# end


end  # module Results
