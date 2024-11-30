"""
RingPoly is a struct for the modulo polynomial arithmetic over a polynomial ring.
"""
struct RingPoly
    coeffs::Vector{UInt64}
    N::Int64

    RingPoly(coeffs::Vector{UInt64}) = new(coeffs, length(coeffs))
    RingPoly(coeffs::Vector{UInt64}, N::Int64) = new(coeffs, N)
    RingPoly(N::Int64) = new(zeros(UInt64, N), N)
end

"""
TransPoly is a struct for the modulo polynomial arithmetic over a polynomial ring.
"""
struct TransPoly
    coeffs::Vector{Float64}
    N::Int64

    TransPoly(coeffs::Vector{Float64}) = new(coeffs, length(coeffs))
    TransPoly(coeffs::Vector{Float64}, N::Int64) = new(coeffs, N)
    TransPoly(N::Int64) = new(zeros(Float64, N), N)
end

const PolyType = Union{RingPoly, TransPoly}

randringpoly!(res::RingPoly, rng::ChaChaStream) = begin
    @inbounds @simd for i = eachindex(res.coeffs)
        res.coeffs[i] = rand(rng, UInt64)
    end 
end

Base.copy!(dst::T, src::T) where {T<:PolyType} = begin
    @assert dst.N == src.N "The ring dimension of the destination polynomial should be same to the input polynomial."

    copy!(dst.coeffs, src.coeffs)
end

neg(x::T) where {T<:PolyType} = begin
    res = T(x.N)
    negto!(res, x)
end

negto!(res::T, x::T) where {T<:PolyType} = begin
    @assert res.N == x.N "Output polynomial should have the same ring degree to the input polynomial."
    
    @inbounds @simd for i = eachindex(res.coeffs)
        res.coeffs[i] = -x.coeffs[i]
    end
end

add(x::T, y::T) where {T<:PolyType} = begin
    res = T(x.N)
    addto!(res, x, y)
end

addto!(res::T, x::T, y::T) where {T<:PolyType} = begin
    @assert x.N == y.N "Input polynomials should have the same ring degree."
    @assert res.N == x.N "Output polynomial should have the same ring degree to the input polynomial."
    
    @inbounds @simd for i = eachindex(res.coeffs)
        res.coeffs[i] = x.coeffs[i] + y.coeffs[i]
    end 
end

sub(x::T, y::T) where {T<:PolyType} = begin    
    res = T(x.N)
    subto!(res, x, y)
end

subto!(res::T, x::T, y::T) where {T<:PolyType} = begin
    @assert x.N == y.N "Input polynomials should have the same ring degree."
    @assert res.N == x.N "Output polynomial should have the same ring degree to the input polynomial."
    
    @inbounds @simd for i = eachindex(res.coeffs)
        res.coeffs[i] = x.coeffs[i] - y.coeffs[i]
    end 
end

rotate(p::RingPoly, idx::Int64) = RingPoly(circshift(p.coeffs, idx))
rotate!(p::RingPoly, idx::Int64) = circshift!(p.coeffs, idx)

rotateto!(res::RingPoly, p::RingPoly, idx::Int64) = begin
    @. res.coeffs = p.coeffs
    circshift!(res.coeffs, idx)
end

initialise!(p::RingPoly) = begin
    @inbounds @simd for i = eachindex(p.coeffs)
        p.coeffs[i] = zero(UInt64)
    end
end

rotate(p::TransPoly, idx::Int64) = TransPoly(circshift(p.coeffs, -idx))

rotate!(p::TransPoly, idx::Int64) = circshift!(p.coeffs, -idx)

rotateto!(res::TransPoly, p::TransPoly, idx::Int64) = begin
    @. res.coeffs = p.coeffs
    circshift!(res.coeffs, -idx)
end

initialise!(p::TransPoly) = begin
    @inbounds @simd for i = eachindex(p.coeffs)
        p.coeffs[i] = .0
    end
end


mul(x::TransPoly, y::TransPoly) = begin
    res = deepcopy(x)
    multo!(res, x, y)
end

multo!(res::TransPoly, x::TransPoly, y::TransPoly) = begin
    @assert x.N == y.N "Input polynomials should have the same length."
    @assert x.N == res.N "Input polynomial and output polynomial should have the same length."

    @inbounds @simd for i = eachindex(res.coeffs)
        res.coeffs[i] = x.coeffs[i] * y.coeffs[i]
    end
end

muladdto!(res::TransPoly, x::TransPoly, y::TransPoly) = begin
    @assert x.N == y.N "Input Polynomials should have the same length."
    @assert x.N == res.N "Input polynomial and output polynomial should have the same length."

    @inbounds @simd for i = eachindex(res.coeffs)
        res.coeffs[i] += x.coeffs[i] * y.coeffs[i]
    end
end

mulsubto!(res::TransPoly, x::TransPoly, y::TransPoly) = begin
    @assert x.N == y.N "Input Polynomials should have the same length."
    @assert x.N == res.N "Input polynomial and output polynomial should have the same length."
   
    @inbounds @simd for i = eachindex(res.coeffs)
        res.coeffs[i] -= x.coeffs[i] * y.coeffs[i]
    end
end

function fft(p::RingPoly, ffter::SubringTransformer)
    @assert p.N == ffter.N "Polynomial length does not match the FFT parameter."
    res = TransPoly(Vector{Float64}(undef, p.N))
    fftto!(res.coeffs, p.coeffs, ffter)
    res
end

function fftto!(res::TransPoly, p::RingPoly, ffter::SubringTransformer)
    @assert res.N == p.N "The input and output polynomial lengths do not match."
    @assert p.N == ffter.N "Polynomial length does not match the FFT parameter."
    fftto!(res.coeffs, p.coeffs, ffter)
end

function ifft(p::TransPoly, ffter::SubringTransformer)
    @assert p.N == ffter.N "Polynomial length does not match the FFT parameter."
    res = RingPoly(Vector{UInt64}(undef, p.N))
    ifftto!(res.coeffs, p.coeffs, ffter)
    res
end

function ifftto!(res::RingPoly, p::TransPoly, ffter::SubringTransformer) 
    @assert res.N == p.N "The input and output polynomial lengths do not match."
    @assert p.N == ffter.N "Polynomial length does not match the FFT parameter."
    ifftto!(res.coeffs, p.coeffs, ffter)
end