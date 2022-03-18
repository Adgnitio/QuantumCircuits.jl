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

using QuantumCircuits
using QuantumCircuits.QML
using QuantumCircuits.Execute


################################################################################
################################################################################
function test_findparam(n, expqc, params=nothing)
    expmat = tomatrix(expqc)

    qc = QCircuit(n)
    for i in 0:(n-1)
        qc.u3(i)
    end
    if params != nothing
        setparameters!(qc, params)
    end
    params = getparameters(qc)

    @test expqc != qc
    @test expqc != simplify(qc)

    # fit the qc to expected
    findparam(expmat, qc)

    expqc != simplify(qc) && println(simplify(qc))
    expqc != simplify(qc) && println(params)
    @test expqc == simplify(qc)
end

################################################################################
n = 2
expqc = QCircuit(n)
expqc.x(0)

test_findparam(n, expqc)

# hard case

test_findparam(n, expqc, [3.266936951898874, 5.70062968255718, 4.413502891655148,
                    6.065806193650715, 4.586495424027681, 1.8946506980151154])

test_findparam(n, expqc, [3.146417882036646, 5.838659944632545, 0.7763087422764434,
                    3.890182287359087, 3.9144393001843962, 0.20934945580756525])

test_findparam(n, expqc, [6.027933728465509, 3.4939251107512956, 0.6586845450453815,
                    3.4399018540096193, 2.2051222189127144, 4.478774379878294])


################################################################################
n = 2
expqc = QCircuit(n)
expqc.x(0)
expqc.x(1)

test_findparam(n, expqc)

# hard case
test_findparam(n, expqc, [0.9904270500489234, 2.9438062464029695, 4.8575515480040545,
                    3.8392956982807447, 2.811189470141951, 5.58973178939339])

# hard case
test_findparam(n, expqc, [5.779528542008802, 4.8766335097778395, 3.6466973636363167,
                    2.9460221010727357, 2.924955576746307, 5.518914714482275])

################################################################################
n = 3
expqc = QCircuit(n)
expqc.z(0)
expqc.y(1)
expqc.x(2)

test_findparam(n, expqc)

################################################################################

n = 2
expqc = QCircuit(n)
expqc.h(0)
expqc.x(1)
expqc.cx(0, 1)
expqc.z(0)

expqc2 = QCircuit(n)
expqc2.h(0)
expqc2.cx(0, 1)
expqc2.z(0)
expqc2.x(1)

qc = QCircuit(n)
for i in 0:(n-1)
    qc.u3(i)
end
qc.cx(0, 1)
for i in 0:(n-1)
    qc.u3(i)
end

# test 1
expmat = tomatrix(expqc)
findparam(expmat, qc)
sqc = simplify(qc)
@test sqc == expqc || sqc == expqc2

# test 2
setparameters!(qc, getRandParameters(qc))
expmat = tomatrix(expqc2)
findparam(expmat, qc)
sqc = simplify(qc)
@test sqc == expqc || sqc == expqc2
