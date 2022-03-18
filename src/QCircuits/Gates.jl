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

module Gates

using QuantumCircuits.QCircuits.QBase
using QuantumCircuits.QCircuits.Math
using QuantumCircuits.QCircuits.Registers

import Base: show, length, inv
import QuantumCircuits.QCircuits.QBase: tomatrix, setparameters!, simplify,
       standardGateError, decompose, bindparameters!

export X, Y, Z, S, Sd, T, Td, H, CX, U3, Rx, Ry, Rz, U, Sx,
      getqubits, getqubitsids, toU3,
      ChromosomeGate, DNAGate, Parameters, getArgs, decompose,
      ParameterT, getvalue


"Parameter type"
const ParamT = Number

mutable struct ParameterT
    value::Union{ParamT, Nothing}
end
ParameterT() = ParameterT(nothing)
ParameterT(param::ParameterT) = ParameterT(param.value)
Base.show(io::IO, p::ParameterT) = print(io, "$(p.value)")

const Parameter = Union{ParamT, ParameterT}

getvalue(p::ParamT) = p
function getvalue(p::ParameterT)
    @assert p.value != nothing "Unable get unset parameter."
    return p.value
end

#const ParamIdx = Int
# "The gates parameters structure, one plece for all parameters in cirquit."
# struct Parameters
#     params::Vector{ParamT}
#
#     Parameters() = new(ParamT[])
# end
# "Function return the next free parameter."
# function getnext(p::Parameters)
#     push!(p.params, zero(ParamT))
#
#     return length(p.params)
# end
# "Get/Set index"
# Base.getindex(p::Parameters, idx::ParamIdx) = p.params[idx]
# Base.setindex!(p::Parameters, v::ParamT, idx::ParamIdx) = p.params[idx] = v

"Gate which can be use in chromosome"
#abstract type ChromosomeGate <: QuantumGate end
#abstract type DNAGate <: ChromosomeGate end

#"Single Qubit Gate"
#abstract type SingleQubitGate <: QuantumGate end

"Quantum Universal Gate"
abstract type UniversalGate <: QuantumGate end


"Comparator"
Base.:(==)(c1::QuantumGate, c2::QuantumGate) = false

"Single qubit create macro"
macro singleQubitGate(name)
    eval(quote
        struct $name <: QuantumGate qubit::Qubit end

        Base.:(==)(g1::$name, g2::$name) = g1.qubit == g2.qubit
        Base.hash(g::$name, h::UInt) = hash(g.qubit, h)

        getArgs(g::$name) = (getid(g.qubit), )
    end)
end

# Single qubit quantum gates
@singleQubitGate(X)
@singleQubitGate(Y)
@singleQubitGate(Z)
@singleQubitGate(S)
@singleQubitGate(Sd)
@singleQubitGate(H)
@singleQubitGate(T)
@singleQubitGate(Td)
@singleQubitGate(Sx)

"CX Quantum Gate"
struct CX <: QuantumGate
    control::Qubit
    target::Qubit
end
Base.:(==)(g1::CX, g2::CX) = g1.control == g2.control && g1.target == g2.target
Base.hash(g::CX, h::UInt) = hash((g.control, g.target), h)
function Base.:(>)(x::CX, y::CX)
    if x.control == y.control
        return x.target > y.target
    else
        return x.control > y.control
    end
end
getArgs(g::CX) = (getid(g.control), getid(g.target))

"Rotation gates"
abstract type RotationGate <: QuantumGate end

"Rotation macro"
macro rotationGate(name)
    eval(quote
        mutable struct $name <: RotationGate
            qubit::Qubit
            θ::Parameter
        end
        $name(qubit::Integer) = $name(qubit, ParameterT(rand() * 2π))

        Base.:(==)(g1::$name, g2::$name) = g1.qubit == g2.qubit && g1.θ == g2.θ
        Base.hash(g::$name, h::UInt) = hash((g.qubit, g.θ), h)

        $name(qubis::Vector{T}, θ::ParamT) where T <: Integer = [$name(q, θ) for q in qubis]
        getArgs(g::$name) = (getid(g.qubit), g.θ)
    end)
end

@rotationGate(Rx)
@rotationGate(Ry)
@rotationGate(Rz)

appendparams!(param, gate::RotationGate) = appendparams!(param, gate.θ)
appendRandParams!(param, gate::RotationGate) = appendRandParams!(param, gate.θ)
Base.length(gate::RotationGate) = getlength(gate.θ)
function setparameters!(gate::RotationGate, params)
    @assert length(params) == length(gate) "RotationGate has $(length(gate)) parameters but $(length(param)) was provided."
    setparameter!(gate.θ, params, 1)
end


"The standard U3 gate"
mutable struct U3 <: UniversalGate
    qubit::Qubit
    θ::Parameter
    ϕ::Parameter
    λ::Parameter
end
U3(qubit::Qubit) = U3(qubit, ([ParameterT(p) for p in rand(3) * 2π])...)
getArgs(g::U3) = (getid(g.qubit), g.θ, g.ϕ, g.λ)

Base.show(io::IO, g::U3) = print(io, "U3($(g.qubit), $(g.θ), $(g.ϕ), $(g.λ))")

Base.:(==)(g1::U3, g2::U3) = g1.qubit == g2.qubit && g1.θ == g2.θ && g1.ϕ == g2.ϕ && g1.λ == g2.λ
Base.hash(g::U3, h::UInt) = hash((g.qubit, g.θ, g.ϕ, g.λ), h)

bindparameters!(gate::QuantumGate) = nothing
bindparameters!(gate::RotationGate) = gate.θ = getvalue(gate.θ)
function bindparameters!(gate::U3)
    gate.θ = getvalue(gate.θ)
    gate.ϕ = getvalue(gate.ϕ)
    gate.λ = getvalue(gate.λ)
end

# "The universar single qubit unitary gate macro"
# macro universalGate(name)
#     eval(quote
#         mutable struct $name <: UniversalGate
#             qubit::Qubit
#             θ::ParamT
#             ϕ::ParamT
#             λ::ParamT
#         end
#         $name(qubit::Integer) = $name(qubit, (rand(3) * 2π)...)
#     end)
# end

"The universar single qubit unitary gate"
mutable struct U <: UniversalGate
    qubit::Qubit
    β::Union{ParamT, Nothing}
    γ::Union{ParamT, Nothing}
    δ::Union{ParamT, Nothing}
end
U(qubit::Integer) = U(qubit, nothing, nothing, nothing)
getArgs(g::U) = (getid(g.qubit), g.β, g.γ, g.δ)


function setparameter!(param::ParameterT, params, i)
    param.value = params[i]
    return i + 1
end
setparameter!(param::ParamT, params, i) = i

function setparameters!(gate::U3, params)
    @assert length(params) == length(gate) "U3 has $(length(gate)) parameters but $(length(params)) was provided."
        i = setparameter!(gate.θ, params, 1)
        i = setparameter!(gate.ϕ, params, i)
        i = setparameter!(gate.λ, params, i)
end

function setparameters!(gate::U, params)
    @assert length(params) == 3 "U has 3 parameters but $(length(params)) was provided."
    gate.β = params[1]
    gate.γ = params[2]
    gate.δ = params[3]
end


"Return the qubits on which operate the gate"
getqubits(gate::QuantumGate) = (gate.qubit,) # All Single Qubit Gates has qubit field
getqubits(gate::CX) = (gate.control, gate.target)

"Return the qubits ids on which operate the gate"
getqubitsids(gate::QuantumGate) = (getid(gate.qubit),) # All Single Qubit Gates has qubit field
getqubitsids(gate::CX) = (getid(gate.control), getid(gate.target))

"Show method"
Base.show(io::IO, gate::QuantumGate) = print(io, "$(typeof(gate))($(gate.qubit))")
Base.show(io::IO, gate::CX) = print(io, "CX($(gate.control), $(gate.target))")



"Convert gate to U3 gate"
toU3(gate::X) = U3(gate.qubit, π, 0, π)


"Return the gate matrix"
u3matrix(θ, ϕ, λ) = [cos(θ/2)+0im              -exp(1im * λ) * sin(θ/2);
                     exp(1im * ϕ) * sin(θ/2)   exp(1im * λ + 1im * ϕ) * cos(θ/2)]
u2matrix(ϕ, λ) = u3matrix(π/2, ϕ, λ)
u1matrix(λ) = u3matrix(0, 0, λ)

# u3matrix(θ, -π/2, π/2)
Rxmatrix(θ) = [cos(θ/2)+0im -1im * sin(θ/2); -1im * sin(θ/2) cos(θ/2)]
# u3matrix(θ, 0, 0)
Rymatrix(θ) = [cos(θ/2)+0im -sin(θ/2); sin(θ/2) cos(θ/2)]
# exp(-1im * θ/2) * u1matrix(θ)
Rzmatrix(θ) = [exp(-1im * θ/2) 0; 0 exp(1im * θ/2)]
# Universary matrix
umatrix(β, γ, δ) = Rzmatrix(β) * Rymatrix(γ) * Rzmatrix(δ)

# Initial stata
const _0 = [1+0im; 0]
const _1 = [0im; 1]

const I0X0I = _0 * dagger(_0)
const I1X1I = _1 * dagger(_1)

# Gate matrixs
const Xmatrix = round.(u3matrix(π, 0, π), digits=2)
const Ymatrix = round.(u3matrix(π, π/2, π/2), digits=2)
const Zmatrix = round.(u1matrix(π), digits=2)
const Imatrix = [1.0 + 0im 0; 0 1]
const Hmatrix = u2matrix(0, π)
const Smatrix = round.(u1matrix(π/2), digits=2)
const Sdmatrix = round.(u1matrix(-π/2), digits=2)
const Tmatrix = u1matrix(π/4)
const Tdmatrix = u1matrix(-π/4)
const CX0matrix = [1.0 + 0im 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0]
const CX1matrix = [1.0 + 0im 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]

tomatrix(gate::X) = Xmatrix
tomatrix(gate::Y) = Ymatrix
tomatrix(gate::Z) = Zmatrix
tomatrix(gate::S) = Smatrix
tomatrix(gate::Sd) = Sdmatrix
tomatrix(gate::T) = Tmatrix
tomatrix(gate::Td) = Tdmatrix
tomatrix(gate::H) = Hmatrix
tomatrix(gate::U3) = u3matrix(getvalue(gate.θ), getvalue(gate.ϕ), getvalue(gate.λ))
function tomatrix(gate::U3, param)
    @assert length(param) == 3 "U3 has 3 parameters but $(length(param)) was provided."
    u3matrix(param[1], param[2], param[3])
end

tomatrix(gate::Rx) = Rxmatrix(getvalue(gate.θ))
tomatrix(gate::Ry) = Rymatrix(getvalue(gate.θ))
tomatrix(gate::Rz) = Rzmatrix(getvalue(gate.θ))
function tomatrix(gate::Rx, param)
    @assert length(param) == 1 "Rx has 1 parameters but $(length(param)) was provided."
    Rxmatrix(param[1])
end
function tomatrix(gate::Ry, param)
    @assert length(param) == 1 "Ry has 1 parameters but $(length(param)) was provided."
    Rymatrix(param[1])
end
function tomatrix(gate::Rz, param)
    @assert length(param) == 1 "Rz has 1 parameters but $(length(param)) was provided."
    Rzmatrix(param[1])
end

tomatrix(gate::U) = umatrix(gate.β, gate.γ, gate.δ)
function tomatrix(gate::U, param)
    @assert length(param) == 3 "U has 3 parameters but $(length(param)) was provided."
    umatrix(param[1], param[2], param[3])
end


function tomatrix(gate::CX)
    control = getid(gate.control)
    target = getid(gate.target)
    if abs(control - target) == 1
        if control < target
            return CX0matrix
        else
            return CX1matrix
        end
    else
        # https://quantumcomputing.stackexchange.com/questions/4252/how-to-derive-the-cnot-matrix-for-a-3-qubit-system-where-the-control-target-qu
        if control < target
            s0 = kron(Imatrix, I0X0I)
            s1 = kron(Imatrix, I1X1I)
            e0 = Imatrix
            e1 = Xmatrix
        else
            s0 = kron(Imatrix, Imatrix)
            s1 = kron(Imatrix, Xmatrix)
            e0 = I0X0I
            e1 = I1X1I
        end

        additional_update = abs(target - control) - 2
        for i in 1:additional_update
            s0 = kron(Imatrix, s0)
            s1 = kron(Imatrix, s1)
        end

        return kron(e0, s0) + kron(e1, s1)
    end
end


function tomatrix(qubits::Integer, gate::QuantumGate, params=nothing)
    # get gate unitary
    if params != nothing
        gate_unitary = tomatrix(gate, params)
    else
        gate_unitary = tomatrix(gate)
    end

    # create kroneker product
    if qubits == 1
        return gate_unitary
    end

    # Get rest
    idxs = getqubitsids(gate)
    minidx = minimum(idxs)
    maxidx = qubits - maximum(idxs) - 1
    #println("idxs: $idxs, minidx: $minidx, maxidx: $maxidx")

    return kron(eye(2^maxidx), gate_unitary, eye(2^minidx))
end

"Method return the error of standard gate"
standardGateError(gate::QuantumGate, params=nothing) = 0.0

function standardGateError(gate::UniversalGate, params=nothing)
    # get gate unitary
    if params != nothing
        m = tomatrix(gate, params)
    else
        m = tomatrix(gate)
    end

    # err = 1.0
    # err *= unitary_error(m, Xmatrix)
    # err *= unitary_error(m, Ymatrix)
    # err *= unitary_error(m, Zmatrix)
    # err *= unitary_error(m, Imatrix)
    # err *= unitary_error(m, Hmatrix)
    # err *= unitary_error(m, Smatrix)
    # err *= unitary_error(m, Sdmatrix)
    # err *= unitary_error(m, Tmatrix)
    # err *= unitary_error(m, Tdmatrix)
    #
    # return log(err + 1.0)

    return minimum([unitary_error(m, Xmatrix),
                    unitary_error(m, Ymatrix),
                    unitary_error(m, Zmatrix),
                    unitary_error(m, Imatrix),
                    unitary_error(m, Hmatrix),
                    unitary_error(m, Smatrix),
                    unitary_error(m, Sdmatrix),
                    unitary_error(m, Tmatrix),
                    unitary_error(m, Tdmatrix)])
end


appendparams!(param, gate) = nothing
appendparams!(param, value::ParameterT) = push!(param, value.value)
appendparams!(param, value::ParamT) = nothing
function appendparams!(param, gate::U3)
    appendparams!(param, gate.θ)
    appendparams!(param, gate.ϕ)
    appendparams!(param, gate.λ)
    nothing
end
function appendparams!(param, gate::U)
    push!(param, gate.β)
    push!(param, gate.γ)
    push!(param, gate.δ)
    nothing
end

appendRandParams!(param, value::ParameterT) = append!(param, rand(1) * 2π)
appendRandParams!(param, value::ParamT) = nothing
appendRandParams!(param, gate) = nothing
function appendRandParams!(param, gate::U3)
    appendRandParams!(param, gate.θ)
    appendRandParams!(param, gate.ϕ)
    appendRandParams!(param, gate.λ)
    nothing
end
appendRandParams!(param, gate::U) = append!(param, rand(3) * 2π)


getlength(value::ParameterT) = 1
getlength(value::ParamT) = 0
Base.length(gate) = 0
Base.length(gate::U3) = getlength(gate.θ) + getlength(gate.ϕ) + getlength(gate.λ)
Base.length(gate::U) = 3

################################################################################
const maxerr = 5e-5

"The library of single qubit gates"
const singlequbitsgates = Dict(Xmatrix => X,
                               Ymatrix => Y,
                               Zmatrix => Z,
                               Smatrix => S,
                               Sdmatrix => Sd,
                               Tmatrix => T,
                               Tdmatrix => Td,
                               Hmatrix => H)
# TODO  + I !!!!!

decompose(gate::QuantumGate) = [gate]

"The simplify the gates."
simplify(gate::QuantumGate) = [gate]

function simplify(gate::UniversalGate)
    exp = tomatrix(gate)

    # check for identify
    if unitary_error(Imatrix, exp) < maxerr
        return []
    end

    # single qubits gate
    for (matrix, mgate) in singlequbitsgates
        if unitary_error(matrix, exp) < maxerr
            # println("found: $(mgate(gate.qubit))")
            return [mgate(gate.qubit)]
        end
    end

    # two qubits gate
    for (matrix, mgate) in singlequbitsgates
        for (matrix2, mgate2) in singlequbitsgates
            # println("$mgate -> err: $(unitary_error(matrix, exp)) == $mgate")
            # println(matrix)
            # println(exp)
            if unitary_error(matrix2 * matrix, exp) < maxerr
                # println("found: $(mgate(gate.qubit))")
                return [mgate(gate.qubit), mgate2(gate.qubit)]
            end
        end
    end

    return [gate]
end


inv(gate::X) = gate
inv(gate::H) = gate
inv(gate::Z) = gate
inv(gate::Y) = gate
inv(gate::Rx) = Rx(gate.qubit, neg(gate.θ))
inv(gate::Ry) = Ry(gate.qubit, neg(gate.θ))
inv(gate::Rz) = Rz(gate.qubit, neg(gate.θ))
inv(gate::CX) = gate
inv(gate::T) = Td(gate.qubit)
inv(gate::S) = Sd(gate.qubit)
inv(gate::Td) = T(gate.qubit)
inv(gate::Sd) = S(gate.qubit)
inv(gate::U3) = U3(gate.qubit, neg(gate.θ), neg(gate.λ), neg(gate.ϕ))

function neg(p::ParameterT)
    @assert p.value != nothing "Unable inv unset parameter."
    return ParameterT(-p.value)
end
neg(value::ParamT) = -value

end  # module Gates
