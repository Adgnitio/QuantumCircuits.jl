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
       QuantumNumber, QuantumInteger, QuantumFloat, QuantumAbstractRegister, getid, setid,
       getIntValue

"The abstract bit"
abstract type Bit end

"The abstract register"
abstract type Register{T} <: AbstractVector{T} end

"Classical bit"
mutable struct Cbit <: Bit
    index::Union{Int, Nothing}
    regIndex::Int

    Cbit(idx::Integer) = new(nothing, idx)
end

"Quantum bit"
mutable struct Qubit <: Bit
    index::Union{Int, Nothing}
    regIndex::Int

    Qubit(idx::Integer) = new(nothing, idx)
end

"The quantum abstract register"
abstract type QuantumAbstractRegister <: Register{Qubit} end


"The show method"
function Base.show(io::IO, q::Qubit)
    if isnothing(q.index)
        print(io, "$(q.regIndex)")
    else
        print(io, "$(q.index)")
    end
end

function Base.:(==)(q1::Qubit, q2::Qubit)
    @assert !isnothing(q1.index) "Unable compare un assigned register."
    @assert !isnothing(q2.index) "Unable compare un assigned register."

    return q1.index == q2.index
end

function Base.hash(q::Qubit, h::UInt)
    @assert !isnothing(q.index) "Unable compare un assigned register."

    return hash(q.index, h)
end


"Get bit index"
function getid(b::Bit)
    @assert !isnothing(b.index) "The bit is not assigned to circuit."

    return b.index
end
"Get bit index"
function setid(b::Bit, id::Integer)
    @assert isnothing(b.index) "Unbale to assig bit to circuit twice."

    b.index = id
end

"Register macro"
macro register(name, bit, basetype)
    eval(quote
        struct $name <: $basetype
            name::Union{String, Nothing}
            bits::Vector{$bit}
            tomeasure::Bool

            $name(n::Integer, name::Union{String, Nothing}; tomeasure::Bool=true) = new(name, [$bit(i) for i in 0:(n-1)], tomeasure)
        end
        $name(n::Integer; tomeasure::Bool=true) = $name(n, nothing, tomeasure=tomeasure)
        function Base.show(io::IO, reg::$name)
            if isnothing(reg.name)
                Base.show(io, string($name) * "($(length(reg)))")
            else
                Base.show(io, string($name) * "($(reg.name), $(length(reg)))")
            end
        end
    end)
end
Base.show(io::IO, ::MIME{Symbol("text/plain")}, reg::Register) = show(io::IO, reg)

# Classical Register
@register(ClassicalRegister, Cbit, Register{Cbit})
# Quantum Register
@register(QuantumRegister, Qubit, QuantumAbstractRegister)

# Quantum numbers state of register
@enum QuantumNumberState Empty QFTbase NormalBase

"The quantum number register"
abstract type QuantumNumber <: QuantumAbstractRegister end


# Quantum Register for store Numbers
mutable struct QuantumInteger <: QuantumNumber
    name::Union{String, Nothing}
    integer::Integer # integer part precision
    bits::Vector{Qubit}
    state::QuantumNumberState
    tomeasure::Bool

    QuantumInteger(integer::Integer, name::Union{String, Nothing}; tomeasure::Bool=true) = 
        new(name,
            integer,
            [Qubit(i) for i in 0:(integer  - 1)],
            Empty,
            tomeasure)
end
QuantumInteger(integer::Integer; tomeasure::Bool=true) = QuantumInteger(integer, nothing, tomeasure=tomeasure)
function Base.show(io::IO, reg::QuantumInteger)
    if isnothing(reg.name)
        Base.show(io, "QuantumInteger($(reg.integer))")
    else
        Base.show(io, "QuantumInteger($(reg.name), $(reg.integer))")
    end
end

ClassicalRegister(reg::QuantumInteger) = ClassicalRegister(length(reg))



# Quantum Register for store Numbers
mutable struct QuantumFloat <: QuantumNumber
    name::Union{String, Nothing}
    integer::Integer # integer part precision
    fractional::Integer # fractional part precision
    bits::Vector{Qubit}
    state::QuantumNumberState
    tomeasure::Bool

    QuantumFloat(integer::Integer, fractional::Integer, name::Union{String, Nothing}; tomeasure::Bool=true) = 
        new(name,
            integer,
            fractional,
            [Qubit(i) for i in 0:(integer + fractional - 1)], Empty, tomeasure)
end
QuantumFloat(integer::Integer, fractional::Integer; tomeasure::Bool=true) = QuantumFloat(integer, fractional, nothing, tomeasure=tomeasure)
QuantumFloat(integer::Integer; tomeasure::Bool=true) = QuantumFloat(integer, 0, nothing, tomeasure=tomeasure)
function Base.show(io::IO, reg::QuantumFloat)
    if isnothing(reg.name)
        Base.show(io, "QuantumNumber($(reg.integer), $(reg.fractional))")
    else
        Base.show(io, "QuantumNumber($(reg.name), $(reg.integer), $(reg.fractional))")
    end
end

function getIntValue(reg::QuantumFloat, num::Number)
    integer_part = Int(floor(num))
    fractional_part = (num - integer_part) * 2^reg.fractional
    fractional_part_int = Int(fractional_part)

    @assert integer_part < 2^reg.integer "The number $num is out of the register with $(reg.integer) qubits."
    @assert abs(fractional_part_int - fractional_part) < 1e-10 "The number $num is out of the register fractional part with $(reg.fractional) qubits."

    return integer_part * 2^reg.fractional + fractional_part_int
end

"The size of register."
Base.length(reg::Register) = length(reg.bits)

"Size"
Base.size(reg::Register) = (length(reg),)

"Each index method"
Base.eachindex(reg::Register) = 0:(length(reg)-1)

Base.checkbounds(reg::Register, I...) = length(I) == 1 && all(i >= 0 && i < length(reg) for i in I[1])

"Get bit from register on given index."
function Base.getindex(reg::Register, idx::Integer)
    @assert idx < length(reg) "The index $idx is out of bound for register $(reg.name) with length $(length(reg))."

    #return getid(reg.bits[idx + 1])
    return reg.bits[idx + 1]
end


"Iterable method"
Base.iterate(reg::Register) = Base.iterate(reg.bits)
Base.iterate(reg::Register, state) = Base.iterate(reg.bits, state)

end  # module Registers
