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

module Graph

using QuantumCircuits.QCircuits.QBase
using QuantumCircuits.QCircuits.Gates
import QuantumCircuits.QCircuits.QBase: add!

export DirectedGraph, to_vector

"Circuit direct graph structure"
struct DirectedGraph{IndexT<:Integer}
    qubits::Int
    start_nodes::Vector{IndexT}
    end_nodes::Vector{IndexT}
    vertices::Vector{QuantumGate}
    vertices_map::Dict{QuantumGate, IndexT}
    edges::Vector{Pair{IndexT, IndexT}}
    #sources::Vector{QuantumGate}
    #targets::Vector{QuantumGate}
    in_edges::Vector{Set{IndexT}}
    out_edges::Vector{Set{Pair{IndexT, IndexT}}}
end

"Create new empyt graph"
function DirectedGraph{IndexT}(qubits::Int) where IndexT<:Integer
    g = DirectedGraph{IndexT}(
                      qubits,
                      zeros(IndexT, qubits),
                      zeros(IndexT, qubits),
                      QuantumGate[],
                      Dict{QuantumGate, IndexT}(),
                      Vector{Pair{IndexT, IndexT}}(),
                      Set{IndexT}[],
                      Vector{Set{Pair{Integer, IndexT}}}[])

    return g
end

"The index of start and end node."
const StartEndNode = 0
"Mapping from qubit number to line number."
const qubitToLine(q) = q + 1

"Add gate to circuit"
function add!(g::DirectedGraph{IndexT}, gate::QuantumGate) where IndexT<:Integer
    qubits = getqubitsids(gate)
    @assert all(q < g.qubits for q in qubits) "Qubits in gate $gate is out of DAG qubits range $(g.qubits)."

    @assert length(g.vertices) == length(g.in_edges) "Inexpected inconsistency between length of vertices and in edges"
    @assert length(g.vertices) == length(g.out_edges) "Inexpected inconsistency between length of vertices and out edges"

    # calculate gate index
    push!(g.vertices, gate)
    idx = length(g.vertices)
    g.vertices_map[gate] = idx

    # add in/out edges sets
    push!(g.in_edges, Set{IndexT}())
    push!(g.out_edges, Set{Pair{Integer, IndexT}}())

    # actualize edge
    for q in qubits
        q = qubitToLine(q) # qiskit qubit start indexes from 0

        last_end = g.end_nodes[q]

        # check if this is first gate on given qubits
        if last_end == StartEndNode
            g.start_nodes[q] = idx
        else
            # Add the output edge
            push!(g.out_edges[last_end], q => idx)
            # TODO perfomance
            if (last_end => idx) ∉ g.edges
                push!(g.edges, last_end => idx)
            end
        end

        g.end_nodes[q] = idx
        push!(g.in_edges[idx], last_end)
    end

    nothing
end

function to_vector(g::DirectedGraph{IndexT})::Vector{QuantumGate} where IndexT<:Integer
    gates = Set{IndexT}([x for x in 1:length(g.vertices)])

    lines = g.start_nodes[1:end]
    code = QuantumGate[]

    while any(lines .!= 0)
        for l in lines
            # if there is somethink to do :)
            if l != StartEndNode && isempty(intersect(gates, g.in_edges[l]))
                # we process all source edges, so we put the gate
                push!(code, g.vertices[l])
                # remove vertice from gates
                pop!(gates, l)
                # Update lines
                expected_qubits = Set(qubitToLine(n) for n in getqubitsids(g.vertices[l]))
                for (q, idx) in g.out_edges[l]
                    lines[q] = idx
                    pop!(expected_qubits, q)
                end
                # cleare the qubits without out edge
                for q in expected_qubits
                    lines[q] = StartEndNode
                end

                break
            end
        end
    end

    return code
end


# "Get index of the vertex in graph"
# vertex_index(vertex::QuantumGate, g::DirectedGraph) = g.verticesmap[vertex]

# AbstractGraph interface
# vertices(g::DirectedGraph) = g.vertices
# edges(g::DirectedGraph) = g.edges
# in_edges(vertex::QuantumGate, g::DirectedGraph{V, E}) where {V, E} = g.inedges[vertex_index(vertex)]
# out_edges(vertex::QuantumGate, g::DirectedGraph{V, E}) where {V, E} = g.outedges[vertex_index(vertex)]


end  # module Graph
