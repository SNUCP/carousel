"""
LWEkey is a struct for LWE secret key s ∈ Zⁿ.
"""
struct LWEkey
    n::Int64
    key::Vector{UInt32}
end

binary_lwekey(rng::ChaChaStream, n::Int64, hw::Int64 = 0) = 
    LWEkey(n, uniform_binary(rng, n, hw) .% UInt32)

ternary_lwekey(rng::ChaChaStream, n::Int64, hw::Int64 = 0) = 
    LWEkey(n, uniform_ternary(rng, n, hw) .% UInt32)

block_binary_lwekey(rng::ChaChaStream, d::Int64, ℓ::Int64) =
    LWEkey(d * ℓ, block_binary(rng, d, ℓ) .% UInt32)

"""
RLWEkey is a struct for the RLWE secret key s ∈ R.
"""
RLWEkey = RingPoly

"""
TransRLWEkey is a struct for Fourier transformed RLWE secret key s ∈ R.
"""
TransRLWEkey = TransPoly

"""
Outputs secret key for RLWE.
If hw = 0, the secret is sampled from a uniform binary distribution.
Otherwise, it outputs a secret key with hamming weight at most hw. 
"""
binary_ringkey(rng::ChaChaStream, N::Int64, hw::Int64 = 0) = RLWEkey(uniform_binary(rng, N, hw) .% UInt64, N)

"""
Outputs secret key for RLWE.
If hw = 0, the secret is sampled from a uniform ternary distribution.
Otherwise, it outputs a secret key with hamming weight at most hw. 
"""
ternary_ringkey(rng::ChaChaStream, N::Int64, hw::Int64 = 0) = RLWEkey(uniform_ternary(rng, N, hw) .% UInt64, N)

"""
Outputs secret key for RLWE.
It takes LWE key as an input, and fills the rest with uniform ternary key.
"""
function partial_ringkey(rng::ChaChaStream, N::Int64, lwekey::LWEkey)
    n = lwekey.n

    if n ≥ N
        RLWEkey(signed.(lwekey.key[1:N]) .% UInt64, N)
    else
        RLWEkey(vcat(lwekey.key .% UInt64, uniform_ternary(rng, N - n) .% UInt64), N)
    end
end