module Gates


using QuantumCircuits
using QuantumCircuits.QML
# using QGen.Quantum
# using QGen.Quantum.Gates
# using QGen.Quantum.QBase
# using QGen.Quantum.Circuit
# using QGen.Quantum.Registers
# using QGen.Quantum.Math
# using QGen.Genetic.CircOpt
# using QGen.Quantum.Execute
# using QGen.Quantum.Graph
# using QGen.Genetic.CircOpt

export ZZ, XX, YY, trotter, trotter2, findU4paramsZZYYXX, trotterU4,
       findU4paramsZZYYXXx2, trotter2U4, trotter2a

function ZZ(qc, q0, q1, t, usePulse=false)
    if usePulse
        qc.h(q1)
        qc.rzx(q0, q1, t)
        qc.x(q0)
        qc.rzx(q0, q1, -t)
        qc.x(q0)
        qc.h(q1)
    else
        qc.cx(q0, q1)
        qc.rz(q1, 2*t)
        qc.cx(q0, q1)
    end
end

function YY(qc, q0, q1, t, usePulse=false)
    if usePulse
        qc.sdg([q0, q1])
        qc.h(q0)
        qc.rzx(q0, q1, t)
        qc.x(q0)
        qc.rzx(q0, q1, -t)
        qc.x(q0)
        qc.h(q0)
        qc.s([q0, q1])
    else
        qc.rx([q0, q1], π/2)
        qc.cx(q0, q1)
        qc.rz(q1, 2*t)
        qc.cx(q0, q1)
        qc.rx([q0, q1], -π/2)
    end
end

function XX(qc, q0, q1, t, usePulse=false)
    if usePulse
        qc.h(q0)
        qc.rzx(q0, q1, t)
        qc.x(q0)
        qc.rzx(q0, q1, -t)
        qc.x(q0)
        qc.h(q0)
    else
        qc.ry([q0, q1], π/2)
        qc.cx(q0, q1)
        qc.rz(q1, 2*t)
        qc.cx(q0, q1)
        qc.ry([q0, q1], -π/2)
    end
end


function trotter(qc, qubits, t)
    for i in 1:(length(qubits)-1)
        ZZ(qc, qubits[i], qubits[i+1], t)
        YY(qc, qubits[i], qubits[i+1], t)
        XX(qc, qubits[i], qubits[i+1], t)
    end
end

function trotter2(qc, qubits, t)
    for i in 1:(length(qubits)-2)
        ZZ(qc, qubits[i], qubits[i+1], t/2)
        YY(qc, qubits[i], qubits[i+1], t/2)
        XX(qc, qubits[i], qubits[i+1], t/2)
    end

    i = length(qubits) - 1
    ZZ(qc, qubits[i], qubits[i+1], t)
    YY(qc, qubits[i], qubits[i+1], t)
    XX(qc, qubits[i], qubits[i+1], t)

    for i in 1:(length(qubits)-2)
        ZZ(qc, qubits[i], qubits[i+1], t/2)
        YY(qc, qubits[i], qubits[i+1], t/2)
        XX(qc, qubits[i], qubits[i+1], t/2)
    end
end
function trotter2a(qc, qubits, t, isFirst, isLast)
    if isFirst
        for i in 1:(length(qubits)-2)
            ZZ(qc, qubits[i], qubits[i+1], t/2)
            YY(qc, qubits[i], qubits[i+1], t/2)
            XX(qc, qubits[i], qubits[i+1], t/2)
        end
    end

    i = length(qubits) - 1
    ZZ(qc, qubits[i], qubits[i+1], t)
    YY(qc, qubits[i], qubits[i+1], t)
    XX(qc, qubits[i], qubits[i+1], t)

    if isLast
        for i in 1:(length(qubits)-2)
            ZZ(qc, qubits[i], qubits[i+1], t/2)
            YY(qc, qubits[i], qubits[i+1], t/2)
            XX(qc, qubits[i], qubits[i+1], t/2)
        end
    else
        for i in 1:(length(qubits)-2)
            ZZ(qc, qubits[i], qubits[i+1], t)
            YY(qc, qubits[i], qubits[i+1], t)
            XX(qc, qubits[i], qubits[i+1], t)
        end
    end
end
# function trotter2a(qc, qubits, t)
#     for i in 1:(length(qubits)-1)
#         ZZ(qc, qubits[i], qubits[i+1], t/2)
#         YY(qc, qubits[i], qubits[i+1], t/2)
#         XX(qc, qubits[i], qubits[i+1], t/2)
#     end
#     for i in (length(qubits)-1):-1:1
#         ZZ(qc, qubits[i], qubits[i+1], t/2)
#         YY(qc, qubits[i], qubits[i+1], t/2)
#         XX(qc, qubits[i], qubits[i+1], t/2)
#     end
# end

function trotterU4(qc, qubits, params)
    for i in 1:(length(qubits)-1)
        qc.u4(qubits[i], qubits[i+1], params)
    end
end


function trotter2U4(qc, qubits, t, isFirst, isLast, params, params2, params3)
    if isFirst
        for i in 1:(length(qubits)-2)
            qc.u4(qubits[i], qubits[i+1], params2)
        end
    end

    i = length(qubits) - 1
    qc.u4(qubits[i], qubits[i+1], params)

    if isLast
        for i in 1:(length(qubits)-2)
            qc.u4(qubits[i], qubits[i+1], params2)
        end
    else
        for i in 1:(length(qubits)-2)
            qc.u4(qubits[i], qubits[i+1], params3)
        end
    end
end
# function trotter2U4(qc, qubits, t, isFirst, isLast, params, params2, params3)
#     if isFirst
#         for i in 1:(length(qubits)-2)
#             qc.u4(qubits[i], qubits[i+1], params2)
#         end
#     end
#
#     i = length(qubits) - 1
#     qc.u4(qubits[i], qubits[i+1], params2)
#     qc.u4(qubits[i], qubits[i+1], params2)
#
#     if isLast
#         for i in 1:(length(qubits)-2)
#             qc.u4(qubits[i], qubits[i+1], params2)
#         end
#     else
#         for i in 1:(length(qubits)-2)
#             qc.u4(qubits[i], qubits[i+1], params3)
#         end
#     end
# end


function findU4paramsZZYYXX(t; debug=false)
    qc = QCircuit(2)
    ZZ(qc, 0, 1, t)
    YY(qc, 0, 1, t)
    XX(qc, 0, 1, t)
    expmat = tomatrix(qc)


    qr = QuantumRegister(2)
    qc = QCircuit(qr)
    qc.u4(qr[0], qr[1])

    params = getRandParameters(qc)
    setparameters!(qc, params)
    qc = decompose(qc)

    params, _, err, _  = findparam(expmat, qc, debug=debug, trystandard=false)

    @assert err < 1e-5 "The error of U gate should be small but it is $err."
    #println(round.(expmat, digits=2))
    #println(round.(tomatrix(qc), digits=2))
    #println(qc)

    return params
end


function findU4paramsZZYYXXx2(t; debug=false)
    qc = QCircuit(2)
    ZZ(qc, 0, 1, t)
    YY(qc, 0, 1, t)
    XX(qc, 0, 1, t)
    ZZ(qc, 0, 1, t)
    YY(qc, 0, 1, t)
    XX(qc, 0, 1, t)
    expmat = tomatrix(qc)


    qr = QuantumRegister(2)
    qc = QCircuit(qr)
    qc.u4(qr[0], qr[1])

    params = getRandParameters(qc)
    setparameters!(qc, params)
    qc = decompose(qc)

    params, _, err, _  = findparam(expmat, qc, debug=debug, trystandard=false)
    @assert err < 1e-5 "The error of U gate should be small but it is $err."

    #println(round.(expmat, digits=2))
    #println(round.(tomatrix(qc), digits=2))
    #println(qc)

    return params
end

end  # module Gate
