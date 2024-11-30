"""
`divbits(a, bit)` computes ⌊a/2ᵇⁱᵗ⌉. 
"""
@inline divbits(a::T, bit::Int64) where T <: Unsigned = begin
    mask = one(T) << (bit - 1)
    a >>> bit + (a & mask) >>> (bit - 1)
end

"""
`native(x::Float64)` computes ⌊x⌉ (mod 2⁶⁴)
"""
native(x::Float64) = begin
    x -= round(x * 5.421010862427522e-20) * 1.8446744073709552e19
    x == 9.223372036854776e18 ? 0x8000000000000000 : round(Int64, x) % UInt64
end

"""
`ord(p, m)` returns order of p modulo m.
"""
function ord(p::Integer, m::Integer)
    res = 1
    acc = p % m
    while acc != 1
        acc = (acc * p) % m
        res += 1
    end

    res
end

@inline zeropadto(a::Vector{T}, n::Int64) where T = begin
    @assert n ≥ length(a)
    vcat(a, zeros(T, n - length(a)))
end

function primitive_root_finder(Q::Int64)
    test = (Q-1) .÷ collect(keys(factor(Dict, Q-1)))

    g = 2
    while true
        if all(powermod.(g, test, Q) .!= 1)
            break
        end
        g += 1

        if g > Q
            error("$Q may not be a prime number.")
        end
    end

    g
end
