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


using Test

using QuantumCircuits.QML.Optimization


testFun1(x) = 1.0*x[1]^2 + 100.0*x[2]^2
dtestFun1(x) = [2.0*x[1], 200.0*x[2]]

const ϵ = 0.000001

val, x, itr = gradientDescent(testFun1, dtestFun1, [10.0, 10.0], α=1.0, maxItr=10000000)
#val, x, itr = gradientDescent(testFun1, dtestFun1, [10.0, 10.0], α=1.0, maxItr=10000000, debug=true)
# Standard Unable
# Momentum Unable
# RMSprop      22
# Adam        384
# Moj           4
@test abs(val) < ϵ
@test abs(transpose(x)*x) < ϵ
#@test itr == 3
@test itr < 20

val, x, itr = gradientDescent(testFun1, dtestFun1, [10.0, 10.0], α=0.1, maxItr=10000000)
# Standard Unable
# Momentum    404
# RMSprop      94
# Adam        532
# Moj           4
@test abs(val) < ϵ
@test abs(transpose(x)*x) < ϵ
#@test itr == 3
@test itr < 20

val, x, itr = gradientDescent(testFun1, dtestFun1, [10.0, 10.0], α=0.01, maxItr=10000000)
# Standard Unable
# Momentum    720
# RMSprop    2896
# Adam       4580
# Moj           4
@test abs(val) < ϵ
@test abs(transpose(x)*x) < ϵ
@test itr == 3

val, x, itr = gradientDescent(testFun1, dtestFun1, [10.0, 10.0], α=0.001, maxItr=10000000)
# Standard  9549
# Momentum  9383
# RMSprop  15849
# Adam     17101
# Moj          3
@test abs(val) < ϵ
@test abs(transpose(x)*x) < ϵ
@test itr == 3

val, x, itr = gradientDescent(testFun1, dtestFun1, [10.0, 10.0], α=0.0001, maxItr=10000000)
# Standard  95561
# Momentum  95398
# RMSprop  108080
# Adam     109290
# Moj           3
@test abs(val) < ϵ
@test abs(transpose(x)*x) < ϵ
@test itr == 3
