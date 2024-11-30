"""
uniform_binary(N, hw = h) outputs a uniform binary Int64 vector with Hamming weight h.
If hw is not given or set to zero, it outputs a truly uniform binary vector.
"""
uniform_binary(rng::ChaChaStream, N::Int64, hw::Int64 = 0) = begin
    if hw == 0
        rand(rng, [0, 1], N)
    else
        res = zeros(Int64, N)
        @inbounds @simd for _ = 1 : hw
            res[ceil(Int64, rand(rng) * N)] = 1
        end
        res
    end
end

"""
uniform_ternary(N, hw = h) outputs a uniform ternary Int64 vector with Hamming weight h.
If hw is not given or set to zero, it outputs a truly uniform ternary vector.
"""
uniform_ternary(rng::ChaChaStream, N::Int64, hw::Int64 = 0) = begin
    if hw == 0
        rand(rng, [-1, 0, 1], N)
    else
        res = zeros(Int64, N)
        tmp = [-1, 1]
        @inbounds @simd for _ = 1 : hw
            res[ceil(Int64, rand(rng) * N)] = rand(rng, tmp)
        end
        res
    end
end

"""
block_binary(d, ℓ) outputs an Int64 vector sampled from the block binary distribution,
with block length ℓ and block number d.
"""
function block_binary(rng::ChaChaStream, d::Int64, ℓ::Int64)
    vec = zeros(Int64, d * ℓ)
    @inbounds @simd for i = 0 : d - 1
        idx = rand(rng, 0 : ℓ)
        if idx != 0
            vec[i * ℓ + idx] = 1
        end
    end
    vec
end

"""
gaussian outputs a random sample from the Gaussian distribution.
"""
gaussian(rng::ChaChaStream, T::Type, σ::Float64) = 
    round(Int64, σ * randn(rng)) % T