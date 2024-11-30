"""
`LWEEncryptor` is a struct for encryption and decryption of LWE-based ciphertexts.
"""
struct LWEEncryptor
    rng::ChaChaStream
    key::LWEkey
    σ::Float64
    LWEEncryptor(rng::ChaChaStream, key::LWEkey, σ::Real) = 
        new(rng, key, σ)
end

"""
`LWEsample(entor)` returns an LWE sample.
"""
LWEsample(entor::LWEEncryptor) = begin
    e = gaussian(entor.rng, UInt32, entor.σ)
    a = rand(entor.rng, UInt32, entor.key.n)
    b = e - reduce(+, a .* entor.key.key)
    LWE(b, a)
end

"""
`LWE_encrypt(m, entor)` returns an LWE encryption of message `m`.
"""
LWE_encrypt(m::UInt32, entor::LWEEncryptor) = begin
    res = LWEsample(entor)
    res.b[] += m
    res
end

"""
`phase(ct, entor)` returns the phase of an LWE encryption `ct`.
"""
phase(ct::LWE, entor::LWEEncryptor) = begin
    reduce(+, ct.a .* entor.key.key) + ct.b[]
end

#=====================================================================================================#

"""
`LEV_encrypt(m, entor, params)` returns an LEV encryption of message `m`.
"""
function LEV_encrypt(m::UInt32, entor::LWEEncryptor, params::DecompParams{UInt32})
    stack = Vector{LWE}(undef, params.len)
    @inbounds for i = eachindex(stack)
        stack[i] = LWE_encrypt(m * params.gvec[i], entor)
    end
    LEV(stack)
end

#=====================================================================================================#

"""
`RLWEEncryptor` is a struct for encryption and decryption of RLWE-based ciphertexts.
"""
struct RLWEEncryptor
    rng::ChaChaStream
    key::TransRLWEkey
    σ::Float64
    ffter::Transformer
    decparams::DecompParams
    ringbuff::Vector{RingPoly}
    transbuff::TransPoly

    RLWEEncryptor(rng::ChaChaStream, key::RLWEkey, σ::Real, ffter::Transformer) = begin
        new(rng, fft(key, ffter), σ, ffter, DecompParams{UInt64}(3, 22), [RingPoly(key.N) for _ = 1 : 3], TransPoly(key.N))
    end
end

"""
`RLWE_sample(entor)` returns an RLWE sample. 
"""
function RLWE_sample(entor::RLWEEncryptor)
    res = RLWE(entor.ffter.N)
    RLWE_sample_to!(res, entor)
    res
end

"""
`RLWE_sample_to!(res, entor)` samples an RLWE sample and stores it in res.
"""
function RLWE_sample_to!(res::RLWE, entor::RLWEEncryptor)
    rng, key, σ, ffter = entor.rng, entor.key, entor.σ, entor.ffter
    decpar, ringbuff, transbuff = entor.decparams, entor.ringbuff, entor.transbuff
    N = ffter.N

    randringpoly!(res.a, rng)
    initialise!(res.b)

    decompto!(ringbuff, res.a, decpar)

    @inbounds for i = 1 : 3
        fftto!(transbuff, ringbuff[i], ffter)
        multo!(transbuff, transbuff, key)
        ifftto!(ringbuff[i], transbuff, ffter)
        @. res.b.coeffs <<= 22
        subto!(res.b, res.b, ringbuff[i])
    end
    
    @inbounds @simd for i = 1 : N
        res.b.coeffs[i] += gaussian(rng, UInt64, σ)
    end
end
 
"""
`RLWE_encrypt(m, entor)` returns an RLWE encryption of `m`.
"""
RLWE_encrypt(m::UInt64, entor::RLWEEncryptor) = begin
    res = RLWE_sample(entor)
    # Note that the η-basis representation of a constant c is identical to (-c, -c, …, -c). 
    @. res.b.coeffs -= m
    res
end

"""
`RLWE_encrypt_a(m, entor)` returns an RLWE encryption of `m ⋅ s`.
"""
RLWE_encrypt_a(m::UInt64, entor::RLWEEncryptor) = begin
    res = RLWE_sample(entor)
    # Note that the η-basis representation of a constant c is identical to (-c, -c, …, -c). 
    @. res.a.coeffs -= m
    res
end

RLWE_encrypt(m::RingPoly, entor::RLWEEncryptor) = begin
    res = RLWE_sample(entor)
    addto!(res.b, res.b, m)
    res
end

RLWE_encrypt_a(m::RingPoly, entor::RLWEEncryptor) = begin
    res = RLWE_sample(entor)
    addto!(res.a, res.a, m)
    res
end

"""
`phase(ct, entor)` returns the phase of an RLWE encryption `ct`.
"""
function phase(ct::RLWE, entor::RLWEEncryptor)
    key, ffter = entor.key, entor.ffter
    decpar, ringbuff, transbuff = entor.decparams, entor.ringbuff, entor.transbuff
    N = ffter.N

    res = RingPoly(N)

    # We utilise the decomposition trick from tfhe-go, to reduce the FFT noise.
    decompto!(ringbuff, ct.a, decpar)

    @inbounds for i = 1 : 3
        fftto!(transbuff, ringbuff[i], ffter)
        multo!(transbuff, transbuff, key)
        ifftto!(ringbuff[i], transbuff, ffter)
        @. res.coeffs <<= 22
        addto!(res, res, ringbuff[i])
    end
    
    addto!(res, res, ct.b)

    res
end

#=====================================================================================================#

"""
`RLEV_encrypt(m, entor)` returns an RLEV encryption of `m`.
"""
function RLEV_encrypt(m::Union{Integer, RingPoly}, entor::RLWEEncryptor, params::DecompParams{UInt64})
    res = RLEV([RLWE(entor.ffter.N) for _ = 1 : params.len])
    RLEV_encrypt_to!(res, m, entor, params)
    res
end

"""
`RLEV_encrypt_to!(res, m, entor)` computes an RLEV encryption of `m` and stores it in `res`.
"""
function RLEV_encrypt_to!(res::RLEV, m::Integer, entor::RLWEEncryptor, params::DecompParams{UInt64})
    m = m % UInt64

    @inbounds for i = 1 : params.len
        RLWE_sample_to!(res.stack[i], entor)
        @. res.stack[i].b.coeffs -= m << params.gveclog[i]
    end
end

function RLEV_encrypt_to!(res::RLEV, m::RingPoly, entor::RLWEEncryptor, params::DecompParams{UInt64})
    @inbounds for i = 1 : params.len
        RLWE_sample_to!(res.stack[i], entor)
        @. res.stack[i].b.coeffs += m.coeffs << params.gveclog[i]
    end
end

"""
`RLEV_encrypt_a(m, entor)` returns an RLEV encryption of `m ⋅ s`.
"""
function RLEV_encrypt_a(m::Union{Integer, RingPoly}, entor::RLWEEncryptor, params::DecompParams{UInt64})
    res = RLEV([RLWE(entor.ffter.N) for _ = 1 : params.len])
    RLEV_encrypt_a_to!(res, m, entor, params)
    res
end

function RLEV_encrypt_a_to!(res::RLEV, m::Integer, entor::RLWEEncryptor, params::DecompParams{UInt64})
    @inbounds for i = 1 : params.len
        RLWE_sample_to!(res.stack[i], entor)
        @. res.stack[i].a.coeffs -= m << params.gveclog[i]
    end
end

function RLEV_encrypt_a_to!(res::RLEV, m::RingPoly, entor::RLWEEncryptor, params::DecompParams{UInt64})
    @inbounds for i = 1 : params.len
        RLWE_sample_to!(res.stack[i], entor)
        @. res.stack[i].a.coeffs += m.coeffs << params.gveclog[i]
    end
end

#=====================================================================================================#

"""
`RGSW_encrypt(m, entor)` returns an RGSW encryption of `m`.
"""
function RGSW_encrypt(m::Union{Integer, RingPoly}, entor::RLWEEncryptor, params::DecompParams{UInt64})
    basketb = RLEV_encrypt(m, entor, params)
    basketa = RLEV_encrypt_a(m, entor, params)
    RGSW(basketb, basketa)
end

function RGSW_encrypt_to!(res::RGSW, m::Union{Integer, RingPoly}, entor::RLWEEncryptor, params::DecompParams{UInt64})
    RLEV_encrypt_to!(res.basketb, m, entor, params)
    RLEV_encrypt_a_to!(res.basketa, m, entor, params)
end