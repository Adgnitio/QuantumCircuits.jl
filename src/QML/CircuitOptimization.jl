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


module CircuitOptimization

using Zygote

using QuantumCircuits.QML.Optimization
using QuantumCircuits.QCircuits.QBase
using QuantumCircuits.QCircuits.Circuit
using QuantumCircuits.QCircuits.Math
using QuantumCircuits.QCircuits.Gates: ParamT

export findparam

const ϵ = 1e-7
const accept_ϵ = 1e-3

function multipleGradientDescent(qc, loss, addGlobalPhase; uϵ=ϵ, N=100, debug=false)
    maxItr=20000
    process = true

    # Get parameters
    params = getparameters(qc)
    if addGlobalPhase
        pushfirst!(params, rand() * 2π)
    end

    # best solution
    bestVal = 1e+5
    bestParams = []

    debug && println("Calculate the derivate.")
    dloss(params) = real(loss'(params))
    debug && dloss(params)

    i = 1
    while process
        debug && println("Calculate gradientDescent for iteration $i.")
        val, xparams, itr = gradientDescent(loss, dloss, params, α=1.0, maxItr=maxItr,
                                      argsArePeriodic=true, isExpectedZero=true, ϵ=uϵ)
        debug && println("Solution with error $val found in $itr iteration.")

        # save best solution
        if val < bestVal
            @assert abs(val - loss(xparams)) <= 1e-4 "Unexpected error in multipleGradientDescent: val=$val, loss=$(loss(xparams)), err=$(abs(val - loss(xparams)))"
            bestVal = val
            bestParams = xparams
        end

        if val < uϵ || i > N
            process = false
        else
            if i % 2 == 0
                params = 2π .- params
            else
                params = getRandParameters(qc)
                if addGlobalPhase
                    pushfirst!(params, rand() * 2π)
                end
            end
        end

        # Next iteration
        i += 1
    end

    if addGlobalPhase
        popfirst!(bestParams)
    end

    return bestVal, bestParams
end

function findparam(expmat, qc; N=100, addGlobalPhase=false, debug=false, trystandard=true)
    #println("getparameters")
    #params = getparameters(qc)

    lossPhase(params) = real(unitary_error(expmat, exp(1im * params[1]) * tomatrix(qc, params[2:end])))
    lossNoPhase(params) = real(unitary_error(expmat, tomatrix(qc, params)))

    if addGlobalPhase
        debug && println("Add global phase to unitary.")
        loss = lossPhase
    else
        loss = lossNoPhase
    end

    val, params = multipleGradientDescent(qc, loss, addGlobalPhase, N=N, debug=debug)
    sge = standardGateError(qc, params)
    params2 = nothing
    val2 = nothing

    if trystandard && val < ϵ && sge > ϵ
        debug && println("====================================================")
        debug && println("The optimal solution uses not standard gate.")

        #params2 = getparameters(qc)
        #debug && println("loss2")
        loss2(params) = real(unitary_error(expmat, tomatrix(qc, params))) + standardGateError(qc, params)
        #debug && loss2(params)

        val2, params2 = multipleGradientDescent(qc, loss2, false, N=N, debug=debug)
        val2std = loss(params2)
        sge2 = standardGateError(qc, params2)

        if (val2std - val) <= accept_ϵ && sge2 < sge
            params = params2
        else
            #@show val2std, val, ϵ, sge2, sge
            println("ERROR: The solution with constraint is worse than without them.")
        end
    end

    debug && println("setparameters")
    setparameters!(qc, params)
    debug && println("end")

    return params, params2, val, val2
end

end  # module CircuitOptimization
