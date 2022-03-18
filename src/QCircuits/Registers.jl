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

module Registers

using QuantumCircuits.QCircuits.QBase

export Bit, Qubit, Cbit, Register, ClassicalRegister, QuantumRegister,
       getid, setid

"The abstract bit"
abstract type Bit end

"The abstract register"
abstract type Register end

"Classical bit"
mutable struct Cbit <: Bit
    index::Union{Integer, Nothing}
    regIndex::Integer

    Cbit(idx::Integer) = new(nothing, idx)
end

"Quantum bit"
mutable struct Qubit <: Bit
    index::Union{Integer, Nothing}
    regIndex::Integer

    Qubit(idx::Integer) = new(nothing, idx)
end

function Base.show(io::IO, q::Qubit)
    if q.index == nothing
        print(io, "$(q.regIndex)")
    else
        print(io, "$(q.index)")
    end
end

function Base.:(==)(q1::Qubit, q2::Qubit)
    @assert q1.index != nothing "Unable compare un assigned register."
    @assert q2.index != nothing "Unable compare un assigned register."

    return q1.index == q2.index
end

function Base.hash(q::Qubit, h::UInt)
    @assert q.index != nothing "Unable compare un assigned register."

    return hash(q.index, h)
end


"Get bit index"
function getid(b::Bit)
    @assert b.index != nothing "The bit is not assigned to circuit."

    return b.index
end
"Get bit index"
function setid(b::Bit, id::Integer)
    @assert b.index == nothing "Unbale to assig bit to circuit twice."

    b.index = id
end

"Register macro"
macro register(name, bit)
    eval(quote
        struct $name <: Register
            name::Union{String, Nothing}
            bits::Vector{$bit}

            $name(n::Integer, name::Union{String, Nothing}) = new(name, [$bit(i) for i in 0:(n-1)])
        end
        $name(n::Integer) = $name(n, nothing)
    end)
end

# Classical Register
@register(ClassicalRegister, Cbit)
# Quantum Register
@register(QuantumRegister, Qubit)


"The size of register."
Base.length(reg::Register) = length(reg.bits)

"Get bit from register on given index."
function Base.getindex(reg::Register, idx::Integer)
    @assert idx < length(reg) "The index $idx is out of bound for register $(reg.name)."

    return reg.bits[idx + 1]
end


"Iterable method"
Base.iterate(reg::Register) = Base.iterate(reg.bits)
Base.iterate(reg::Register, state) = Base.iterate(reg.bits, state)

end  # module Registers
