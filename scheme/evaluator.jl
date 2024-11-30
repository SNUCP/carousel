abstract type Evaluator end

"""
`LWEEvaluator` is a struct for LWE key-switching.
"""
struct LWEEvaluator <: Evaluator
    buffdec::Vector{UInt32}
    decpar::DecompParams{UInt32}

    LWEEvaluator(decpar::DecompParams{UInt32}) =
        new(Vector{UInt32}(undef, decpar.len), decpar)
end

"""
`keyswitch(x, y, eval)` computes the key-switching of an LWE ciphertext `x`. 
"""
keyswitch(x::LWE, y::Vector{LEV}, eval::LWEEvaluator) = begin
    res = LWE(x.n - length(y))
    keyswitchto!(res, x, y, eval)
    res
end

"""
`keyswitchto!(res, x, y, eval)` computes the key-switching of the LWE ciphertext `x` and stores it in `res`.
"""
function keyswitchto!(res::LWE, x::LWE, y::Vector{LEV}, eval::LWEEvaluator)
    # LMSS key-switching
    @assert x.n - res.n == length(y) "The dimension of input and output ciphertexts do not match."
    if length(y) ≠ 0
        @assert res.n == y[1].stack[1].n "The dimension of output ciphertext and key-switching key do not match."
        @assert y[1].len == eval.decpar.len "The size of key-switching key and the evaluator do not match."
    end

    decpar, buff = eval.decpar, eval.buffdec
    n = res.n

    @inbounds @simd for i = 1 : n
        res.a[i] = x.a[i]
    end
    res.b[] = x.b[]

    @inbounds for i = eachindex(y)
        decompto!(buff, x.a[n+i], decpar)
        for j = 1 : decpar.len
            muladdto!(res, y[i].stack[j], buff[j])
        end
    end
end

"""
`RingEvaluator` is a struct for external product and key-switching for RLWE ciphertexts.
"""
struct RingEvaluator <: Evaluator
    ffter::Transformer
    transb::TransPoly; transa::TransPoly
    bdec::Vector{RingPoly}; adec::Vector{RingPoly}
    transbdec::Vector{TransPoly}; transadec::Vector{TransPoly}
    decpar::DecompParams{UInt64}

    function RingEvaluator(ffter::Transformer, decpar::DecompParams{UInt64})
        N = ffter.N

        transb, transa = TransPoly(N), TransPoly(N)
        len = max(decpar.len, 3)
        bdec, adec = [RingPoly(N) for _ = 1 : len], [RingPoly(N) for _ = 1 : len]
        transbdec, transadec = [TransPoly(N) for _ = 1 : len], [TransPoly(N) for _ = 1 : len]
        
        new(ffter, transb, transa, bdec, adec, transbdec, transadec, decpar)
    end
end

"""
`gadgetprod(x, y, eval)` computes and returns the gadget product of `x` and `y`.
"""
gadgetprod(x::RingPoly, y::TransRLEV, eval::RingEvaluator) = begin
    res = TransRLWE(x.N)
    gadgetprodto!(res, x, y, eval)
    res
end

"""
`gadgetprodto!(x, y, eval)` computes the gadget product of `x` and `y` and stores it in `res`.
"""
function gadgetprodto!(res::TransRLWE, x::RingPoly, y::TransRLEV, eval::RingEvaluator)
    @assert x.N == y.stack[1].N "The dimension of input polynomial and ciphertext do not match."
    @assert res.N == x.N "The dimension of input and output ciphertexts do not match."
    @assert y.len == eval.decpar.len "The parameter of RLEV ciphertext and evaluator do not match."

    adec, transadec, decpar, ffter = eval.adec, eval.transadec, eval.decpar, eval.ffter

    decompto!(adec, x, decpar)

    initialise!(res)
    
    @inbounds for i = 1 : y.len
        fftto!(transadec[i], adec[i], ffter)
        muladdto!(res, transadec[i], y.stack[i])
    end
end

"""
`keyswitch(x, y, eval)` computes and returns the key-switch of `x` with key-switching key `y`.
"""
keyswitch(x::RLWE, y::TransRLEV, eval::RingEvaluator) = begin
    res = TransRLWE(x.N)
    keyswitchto!(res, x, y, eval)
    res
end

"""
`keyswitchto!(res, x, y, eval)` computes the key-switch of `x` with key-switching key `y` and stores it to `res`.
"""
function keyswitchto!(res::TransRLWE, x::RLWE, y::TransRLEV, eval::RingEvaluator)
    @assert x.N == y.stack[1].N "The dimension of input polynomial and ciphertext do not match."
    @assert res.N == x.N "The dimension of input and output ciphertexts do not match."
    @assert y.len == eval.decpar.len "The parameter of RLEV ciphertext and evaluator do not match."

    adec, transadec, decpar, ffter = eval.adec, eval.transadec, eval.decpar, eval.ffter

    initialise!(res.a)
    fftto!(res.b, x.b, ffter)
    
    decompto!(adec, x.a, decpar)
    @inbounds for i = 1 : decpar.len
        fftto!(transadec[i], adec[i], ffter)
        muladdto!(res, transadec[i], y.stack[i])
    end
end

"""
`decomptobuff_a!(x, eval)` decompose the mask `a` of input RLWE ciphertext `x` to the buffer of the evaluator.  
"""
function decomptobuff_a!(x::RLWE, eval::RingEvaluator)
    ffter, tmpb, adec, transadec, decpar = eval.ffter, eval.transb, eval.adec, eval.transadec, eval.decpar

    fftto!(tmpb, x.b, ffter)

    decompto!(adec, x.a, decpar)

    @inbounds for i = 1 : decpar.len
        fftto!(transadec[i], adec[i], ffter)
    end
end

"""
`hoistedkeyswitchto!(res, y, eval)` performs hoisted key-switching from the buffer.
"""
function hoistedkeyswitchto!(res::TransRLWE, y::TransRLEV, eval::RingEvaluator)
    @assert y.len == eval.decpar.len "The parameter of RLEV ciphertext and evaluator do not match."

    copy!(res.b, eval.transb)
    initialise!(res.a)
    
    @inbounds for i = 1 : eval.decpar.len
        muladdto!(res, eval.transadec[i], y.stack[i])
    end
end

"""
`rotate(x, idx, ksk, eval)` computes and returns the `idx`-th rotation of an RLWE ciphertext `x`.
"""
rotate(x::RLWE, idx::Int64, ksk::TransRLEV, eval::RingEvaluator) = begin
    res = TransRLWE(x.N)
    rotateto!(res, x, idx, ksk, eval)
    res
end

"""
`rotateto!(res, x, idx, ksk, eval)` computes the `idx`-th rotation of an RLWE ciphertext `x` and stores it in `res`.
"""
rotateto!(res::TransRLWE, x::RLWE, idx::Int64, ksk::TransRLEV, eval::RingEvaluator) = begin
    keyswitchto!(res, x, ksk, eval)
    rotate!(res.a, idx); rotate!(res.b, idx)
end

"""
`hoistedrotate!(res, idx, ksk, eval)` computes the hoisted `idx`-th rotation and stores it in `res`.
"""
function hoistedrotateto!(res::TransRLWE, idx::Int64, ksk::TransRLEV, eval::RingEvaluator)
    hoistedkeyswitchto!(res, ksk, eval)
    rotate!(res.a, idx); rotate!(res.b, idx)
end

"""
`externalprod(x, y, eval)` computes and returns the external product of an RLWE ciphertext `x` and an RGSW ciphertext `y`. 
"""
externalprod(x::RLWE, y::TransRGSW, eval::RingEvaluator) = begin
    res = TransRLWE(x.N)
    externalprodto!(res, x, y, eval)
    res
end

"""
`externalprodto!(res, x, y, eval)` computes the external product of an RLWE ciphertext `x` and an RGSW ciphertext `y` and stores it in `res`.
"""
function externalprodto!(res::TransRLWE, x::RLWE, y::TransRGSW, eval::RingEvaluator)
    @assert x.N == y.basketb.stack[1].N "The dimension of input ciphertexts do not match."
    @assert res.N == x.N "The dimension of input and output ciphertexts do not match."
    @assert y.basketb.len == eval.decpar.len "The parameter of RLEV ciphertext and evaluator do not match."

    ffter, decpar = eval.ffter, eval.decpar
    bdec, adec, transbdec, transadec = eval.bdec, eval.adec, eval.transbdec, eval.transadec

    decompto!(bdec, x.b, decpar)
    decompto!(adec, x.a, decpar)

    initialise!(res)

    @inbounds for i = 1 : decpar.len
        fftto!(transbdec[i], bdec[i], ffter)
        fftto!(transadec[i], adec[i], ffter)
        muladdto!(res, transbdec[i], y.basketb.stack[i])
        muladdto!(res, transadec[i], y.basketa.stack[i])
    end
end

"""
`decomptobuff!(x, eval)` decompose the input RLWE ciphertext `x` to the buffer of the evaluator.  
"""
function decomptobuff!(x::RLWE, eval::RingEvaluator)
    ffter, decpar = eval.ffter, eval.decpar
    bdec, adec, transbdec, transadec = eval.bdec, eval.adec, eval.transbdec, eval.transadec

    decompto!(bdec, x.b, decpar)
    decompto!(adec, x.a, decpar)

    @inbounds for i = 1 : decpar.len
        fftto!(transbdec[i], bdec[i], ffter)
        fftto!(transadec[i], adec[i], ffter)
    end
end

"""
`hoistedkeyswitchto!(res, y, eval)` performs hoisted external product from the buffer.
"""
function hoistedexternalprodto!(res::TransRLWE, y::TransRGSW, eval::RingEvaluator)
    @assert res.N == y.basketb.stack[1].N "The dimension of input and output ciphertexts do not match."
    @assert y.basketb.len == eval.decpar.len "The parameter of RLEV ciphertext and evaluator do not match."

    transbdec, transadec = eval.transbdec, eval.transadec

    initialise!(res)
    
    @inbounds for i = 1 : eval.decpar.len
        muladdto!(res, transbdec[i], y.basketb.stack[i])
        muladdto!(res, transadec[i], y.basketa.stack[i])
    end
end

function fftto!(res::TransRLWE, x::RLWE, eval::RingEvaluator)
    fftto!(res.b, x.b, eval.ffter)
    fftto!(res.a, x.a, eval.ffter)
end

function ifftto!(res::RLWE, x::TransRLWE, eval::RingEvaluator)
    ifftto!(res.b, x.b, eval.ffter)
    ifftto!(res.a, x.a, eval.ffter)
end