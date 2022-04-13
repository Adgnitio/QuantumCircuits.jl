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


module Optimization

export gradientDescent, OptimizationFunction

"Standard gradient descent"
struct StandardGD{T <: AbstractFloat}
    α::T
end
function (m::StandardGD)(x, df)
    for i in 1:length(x)
        x[i] = x[i] - m.α * df[i]
    end
end

"Momentum optimization algoritm"
mutable struct Momentum{T <: AbstractFloat}
    α::T
    β::T
    V::Vector{T}

    Momentum{T}(alfa::T, beta::T, len) where {T <: AbstractFloat} = new(alfa, beta, zeros(T, len))
end
function (m::Momentum)(x, df)
    for i in 1:length(x)
        m.V[i] = m.β*m.V[i] + (1.0 - m.β)*df[i]
        x[i] = x[i] - m.α * m.V[i]
    end
end

"RMSprop optimization algoritm"
mutable struct RMSprop{T <: AbstractFloat}
    α::T
    β::T
    S::Vector{T}

    RMSprop{T}(alfa::T, beta::T, len) where {T <: AbstractFloat} = new(alfa, beta, zeros(T, len))
end
function (m::RMSprop)(x, df)
    for i in 1:length(x)
        m.S[i] = m.β*m.S[i] + (1.0 - m.β)*df[i]^2
        x[i] = x[i] - m.α * df[i] / (sqrt(m.S[i]) + 0.00000001)
    end
end

"Adam optimization algoritm"
mutable struct Adam{T <: AbstractFloat}
    α::T
    β1::T
    β2::T
    V::Vector{T}
    S::Vector{T}
    t::Int

    Adam{T}(alfa::T, beta1::T, beta2::T, len) where {T <: AbstractFloat} = new(alfa, beta1, beta2, zeros(T, len), zeros(T, len), 1)
end
function (a::Adam)(x, df)
    for i in 1:length(x)
        a.V[i] = a.β1*a.V[i] + (1.0 - a.β1)*df[i]
        a.S[i] = a.β2*a.S[i] + (1.0 - a.β2)*df[i]^2

        V = a.V[i] / (1.0 - a.β1^a.t)
        S = a.S[i] / (1.0 - a.β2^a.t)

        x[i] = x[i] - a.α * V / (sqrt(S) + 0.00000001)
    end
    a.t += 1
end

const Period = 2*π
const PeriodThreshold = 1.5 * Period
const FinalPeriodThreshold = 1.05 * Period


"Evolutionary Alpha optimization algoritm"
mutable struct Eva{T <: AbstractFloat}
#mutable struct Eva{T}
    α::Vector{T}
    df::Vector{T}
    x::Vector{T}
    val::T
    doUpdate::Bool
    lastSingularCount::Int

    Eva{T}(alfa::T, len) where {T <: AbstractFloat} = new(fill(alfa, len), zeros(T, len), zeros(T, len), 0, false, 0)
    #Eva{T}(alfa::T, len) where {T} = new(fill(alfa, len), zeros(T, len), zeros(T, len), false)
end
function (e::Eva)(x, val, df; argsArePeriodic=false, debug=false, useBigValInc=false)
    ϵ = 1e-20
    ϵ2 = 1e-10
    BigValIncreas = 1.1

    #println("==========================")
    #println("Val: $val, e.val: $(e.val)")
    maxH = 0.0
    maxDiffX = 0.0
    maxα = 0.0
    debug && val > e.val * BigValIncreas && e.doUpdate && useBigValInc && println("Big value increase: val: $val, e.val: $(e.val), max(α): $(maximum(e.α))")
    for i in 1:length(x)
        notsingular = true
        if e.doUpdate
            #i == 1 && @show i, df[i], e.df[i], x[i], e.x[i], e.α[i]

            diffx = x[i] - e.x[i]
            maxDiffX = max(maxDiffX, abs(diffx))
            if e.df[i] * df[i] < 0 && val > e.val  && abs(e.df[i]) > ϵ2 && abs(df[i]) > ϵ2
                # we are in singular point (there is no partial derivate)
                notsingular = false
                α = e.α[i]
            elseif e.lastSingularCount > 0
                α = e.α[i]
                e.lastSingularCount -= 1
            elseif abs(diffx) > ϵ
                H = (df[i] - e.df[i]) / diffx
                #@show H
                α = 1.0 / (abs(H) + ϵ)
                e.α[i] = α

                maxH = max(maxH, abs(H))
            else
                # the d'' = 0, so we are in sedle point, try to pertube some
                α = 1.2 * e.α[i]
                #α = e.α[i]
            end
        else
            α = e.α[i]
        end
        #@show α, x[i], df[i]
        #debug && !notsingular && println("$i singular: x[i]: $(x[i]), e.x[i]: $(e.x[i]), singular: df[i]: $(df[i]), e.df[i]: $(e.df[i])")

        # There was big increas of value
        if val > e.val * BigValIncreas && e.doUpdate && useBigValInc
            #debug && println("$i big diff: x[i]: $(x[i]), e.x[i]: $(e.x[i]), singular: df[i]: $(df[i]), e.df[i]: $(e.df[i])")
            e.α[i] = min(e.α[i] / 2, 0.1)
            x[i] = e.x[i] - e.α[i] * e.df[i]
            val = e.val
        else
            e.df[i] = df[i]
            if notsingular
                maxα = max(maxα, α)
                if argsArePeriodic
                    if x[i] > PeriodThreshold || x[i] < -PeriodThreshold
                        x[i] = x[i] % Period
                    end
                end
                e.x[i] = x[i]
                x[i] = x[i] - α * df[i]
            else
                ix = x[i]
                x[i] = (x[i] + e.x[i]) / 2
                e.x[i] = ix
                e.lastSingularCount = 10
                e.α[i] = min(0.1, α/2)
            end
        end
    end
    e.val = val
    if e.doUpdate
        return (maxH, maxα, maxDiffX)
    else
        e.doUpdate = true
        return (nothing, maxα, maxDiffX)
    end
end


function adjustPeriodicArgs(x)
    x[x .> FinalPeriodThreshold] .-= Period
    x[x .< -FinalPeriodThreshold] .+= Period

    return x
end

struct OptimizationFunction
    initial_check::Bool
    f_and_df::Function
    f::Function
end

function gradientDescent(f, df, x; α=0.001, maxItr=nothing, ϵ=1e-8,
                         argsArePeriodic=false, isExpectedZero=false,
                         debug=false, useBigValInc=false)
    of = OptimizationFunction(true, (x) -> (f(x), df(x)), f)
    return gradientDescent(of, x; α=α, maxItr=maxItr, ϵ=ϵ,
                             argsArePeriodic=argsArePeriodic, isExpectedZero=isExpectedZero,
                             debug=debug, useBigValInc=useBigValInc)
end


function gradientDescent(of::OptimizationFunction, x; α=0.001, maxItr=nothing, ϵ=1e-8,
                         argsArePeriodic=false, isExpectedZero=false,
                         debug=false, useBigValInc=false, checkFn=nothing)
    process = true

    #opt = StandardGD(α)
    #opt = Momentum{Float64}(α, 0.9, length(x))
    #opt = RMSprop{Float64}(α, 0.999, length(x))
    #opt = Adam{Float64}(α, 0.9, 0.999, length(x))
    opt = Eva{Float64}(α, length(x))
    #opt = Eva{ComplexF64}(α, length(x))

    if of.initial_check
        val = of.f(x)
        if isExpectedZero && val < ϵ
            return val, x, 0
        end
    else
        val = 10e10
    end

    # The best values
    best_x = x[:]
    best_val = val

    itr = 1
    while process
        # calc f and derivative
        val, dval = of.f_and_df(x)

        # save the best
        if best_val > val
            best_val = val
            best_x = x[:]
            debug && println("====================================")
            debug && println(join(best_x, ", "))
            debug && println("====================================")
        end

        #println("val: $val, x[1]: $(x[1]), dval[1]: $(dval[1])")
        # do update
        debug && println("Start iteration $itr, value: $val, best: $best_val, dval: $(maximum(abs.(dval)))")
        debug && checkFn != nothing && println("Real value: $(checkFn(x))")
        H, α, diffx = opt(x, val, dval, argsArePeriodic=argsArePeriodic, debug=debug, useBigValInc=useBigValInc)
        # for i in 1:length(x)
        #     x[i] = x[i] - α * dval[i]
        # end
        #@show x, val, dval
        # some debug
        itr % 1000 == 0 && println("After iteration $itr, value: $val, best: $best_val, dval: $(maximum(abs.(dval))), H: $H, α: $α, diffx: $diffx")
        #debug && println("After iteration $itr, H: $H, α: $α, diffx: $diffx")

        # Check if we have solution
        # check val if the zero is expected
        # @show all(abs.(dval) .< ϵ), any([isnan(p) for p in x]), (isExpectedZero && val < ϵ)
        if all(abs.(dval) .< ϵ) || any(isnan(p) for p in x) || (isExpectedZero && val < ϵ)
            process = false
        end
        if diffx < 2*ϵ && itr > 2000
            println("ERROR: We are going to slow, diffx is equal to: $diffx")
            process = false
        end

        # Check the iterate count
        itr += 1
        if maxItr != nothing && itr > maxItr
            process = false
        end
    end

    if argsArePeriodic
        x = adjustPeriodicArgs(x)
    end

    val = of.f(x)
    # save the best
    if best_val > val
        best_val = val
        best_x = x[:]
        debug && println("====================================")
        debug && println(join(best_x, ", "))
        debug && println("====================================")
    end
    itr -= 1
    #@show val, x, itr

    debug && println("Find solution in iteration $itr, value: $best_val")
    return best_val, best_x, itr
end

end  # module Optimization
