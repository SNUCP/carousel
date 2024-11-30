"""
LWE is a struct for the LWE cryptosystem. Modulus is fixed by 2³².
"""
struct LWE
    n::Int64
    b::RefUInt32
    a::Vector{UInt32}

    LWE(b::UInt32, a::Vector{UInt32}) = new(length(a), Ref(b), a)
    LWE(n::Int64) = new(n, Ref(zero(UInt32)), Vector{UInt32}(undef, n))    
end

add(x::LWE, y::LWE) = 
    LWE(x.b[] + y.b[], x.a + y.a)

addto!(res::LWE, x::LWE, y::LWE) = begin
    res.b[] = x.b[] + y.b[]
    @. res.a = x.a + y.a
end

sub(x::LWE, y::LWE) = 
    LWE(x.b[] - y.b[], x.a - y.a)

subto!(res::LWE, x::LWE, y::LWE) = begin
    res.b[] = x.b[] - y.b[]
    @. res.a = x.a - y.a
end

muladdto!(res::LWE, x::LWE, y::UInt32) = begin
    res.b[] += x.b[] * y
    @. res.a += x.a * y
end

initialise!(x::LWE) = begin
    x.b[] = 0
    @. x.a = 0
end

###########################################################################################

"""
LEV is a struct for the LEV (Gadget Encryption) cryptosystem. Modulus is fixed by 2³².
"""
struct LEV
    len::Int64
    stack::Vector{LWE}

    LEV(stack::Vector{LWE}) = new(length(stack), stack)
end

add(x::LEV, y::LEV) = 
    LEV((@. add(x.stack, y.stack)))

addto!(res::LEV, x::LEV, y::LEV) = begin
    @inbounds @simd for i = 1 : res.len
        addto!(res.stack[i], x.stack[i], y.stack[i])
    end
end

sub(x::LEV, y::LEV) = 
    LEV((@. sub(x.stack, y.stack)))

subto!(res::LEV, x::LEV, y::LEV) = begin
    @inbounds @simd for i = 1 : res.len
        subto!(res.stack[i], x.stack[i], y.stack[i])
    end
end

initialise!(x::LEV) = begin
    @inbounds @simd for i = 1 : x.l
        initialise!(x.stack[i])
    end
end

###########################################################################################

"""
RLWE is a struct for the RLWE cryptosystem. Modulus is fixed by 2⁶⁴.
"""
struct RLWE
    N::Int64
    b::RingPoly
    a::RingPoly

    RLWE(b::RingPoly, a::RingPoly) = begin
        @assert b.N == a.N "The polynomial degree of `b` and `a` are different."
        new(b.N, b, a)
    end

    RLWE(N::Int64) = new(N, RingPoly(N), RingPoly(N))
    RLWE(b::RingPoly) = new(b.N, deepcopy(b), RingPoly(b.N))
end

"""
TransRLWE is a struct for the Fourier transformed RLWE cryptosystem.
"""
struct TransRLWE
    N::Int64
    b::TransPoly
    a::TransPoly

    TransRLWE(b::TransPoly, a::TransPoly) = begin
        @assert b.N == a.N "The polynomial degree of `b` and `a` are different."
        new(b.N, b, a)
    end
    
    TransRLWE(N::Int64) = new(N, TransPoly(N), TransPoly(N))
    TransRLWE(b::TransPoly) = new(b.N, deepcopy(b), TransPoly(b.N))
end

const RLWEType = Union{RLWE, TransRLWE}

Base.copy!(dst::T, src::T) where {T<:RLWEType} = begin
    @assert dst.N == src.N "The size of the ciphertexts `dst` and `src` should match. "
    copy!(dst.b, src.b)
    copy!(dst.a, src.a)
end

add(x::T, y::T) where {T<:RLWEType} =
    T(add(x.b, y.b), add(x.a, y.a))

addto!(res::T, x::T, y::T) where {T<:RLWEType} = begin
    addto!(res.b, x.b, y.b)
    addto!(res.a, x.a, y.a)
end

sub(x::T, y::T) where {T<:RLWEType} = 
    T(sub(x.b, y.b), sub(x.a, y.a))

subto!(res::T, x::T, y::T) where {T<:RLWEType} = begin
    subto!(res.b, x.b, y.b)
    subto!(res.a, x.a, y.a)
end

mul(x::Integer, y::T) where {T<:RLWEType} = 
    T(x * y.b, x * y.a)

multo!(res::T, x::Integer, y::T) where {T<:RLWEType} = begin
    @. res.b = x * y. b
    @. res.a = x * y. a
end

muladdto!(res::T, x::Integer, y::T) where {T<:RLWEType} = begin
    @. res.b += x * y. b
    @. res.a += x * y. a
end

mulsubto!(res::T, x::Integer, y::T) where {T<:RLWEType} = begin
    @. res.b -= x * y. b
    @. res.a -= x * y. a
end

initialise!(x::RLWEType) = begin
    initialise!(x.b)
    initialise!(x.a)
end

mul(x::TransPoly, y::TransRLWE) = 
    TransRLWE(mul(x, y.b), mul(x, y.a))

multo!(res::TransRLWE, x::TransPoly, y::TransRLWE) = begin
    multo!(res.b, x, y.b)
    multo!(res.a, x, y.a)
end

muladdto!(res::TransRLWE, x::TransPoly, y::TransRLWE) = begin
    muladdto!(res.b, x, y.b)
    muladdto!(res.a, x, y.a)
end

mulsubto!(res::TransRLWE, x::TransPoly, y::TransRLWE) = begin
    mulsubto!(res.b, x, y.b)
    mulsubto!(res.a, x, y.a)
end

fft(ct::RLWE, ffter::Transformer) =
    TransRLWE(fft(ct.b, ffter), fft(ct.a, ffter))

fftto!(res::TransRLWE, ct::RLWE, ffter::Transformer) = begin
    fftto!(res.b, ct.b, ffter)
    fftto!(res.a, ct.a, ffter)
end

ifft(ct::TransRLWE, ffter::Transformer) = 
    RLWE(ifft(ct.b, ffter), ifft(ct.a, ffter))

ifftto!(res::RLWE, ct::TransRLWE, ffter::Transformer) = begin
    ifftto!(res.b, ct.b, ffter)
    ifftto!(res.a, ct.a, ffter)
end

###########################################################################################

"""
RLEV is a struct for the RLEV (Gadget Encryption) cryptosystem. Modulus is fixed by 2⁶⁴.
"""
struct RLEV
    len::Int64
    stack::Vector{RLWE}

    function RLEV(stack::Vector{RLWE}) 
        new(length(stack), stack)
    end
end

"""
TransRLEV is a struct for the Fourier transformed RLEV (Gadget Encryption) cryptosystem. 
"""
struct TransRLEV
    len::Int64
    stack::Vector{TransRLWE}

    function TransRLEV(stack::Vector{TransRLWE}) 
        new(length(stack), stack)
    end
end

const RLEVType = Union{RLEV, TransRLEV}

Base.copy!(dst::T, src::T) where {T<:RLEVType} = begin
    @assert dst.len == src.len "The size of the ciphertexts `dst` and `src` should match. "
    @inbounds for i = 1 : src.len
        copy!(dst.stack[i], src.stack[i])
    end
end

add(x::T, y::T) where {T<:RLEVType} = 
    T((@. add(x.stack, y.stack)))

addto!(res::T, x::T, y::T) where {T<:RLEVType} = begin
    @assert x.len == y.len "The size of the input ciphertexts should match."
    @assert x.len == res.len "The size of the input and output ciphertext should match."
    @inbounds for i = 1 : res.len
        addto!(res.stack[i], x.stack[i], y.stack[i])
    end
end

sub(x::T, y::T) where {T<:RLEVType} = 
    T((@. sub(x.stack, y.stack)))

subto!(res::T, x::T, y::T) where {T<:RLEVType} = begin
    @assert x.len == y.len "The size of the input ciphertexts should match."
    @assert x.len == res.len "The size of the input and output ciphertext should match."
    @inbounds @simd for i = 1 : res.len
        subto!(res.stack[i], x.stack[i], y.stack[i])
    end
end

rotate!(x::RLEVType, idx::Int64) = begin
    @inbounds @simd for i = 1 : x.len
        rotate!(x.stack[i].a, idx)
        rotate!(x.stack[i].b, idx)
    end
end

initialise!(x::RLEVType) = begin
    @inbounds @simd for i = 1 : x.len
        initialise!(x.stack[i])
    end
end

fft(ct::RLEV, ffter::Transformer) =
    TransRLEV((fft.(ct.stack, Ref(ffter))))

fftto!(res::TransRLEV, ct::RLEV, ffter::Transformer) = begin
    @inbounds @simd for i = 1 : res.len
        fftto!(res.stack[i], ct.stack[i], ffter)
    end
end

ifft(ct::TransRLEV, ffter::Transformer) = begin
    RLEV((ifft.(ct.stack, Ref(ffter))))
end

ifftto!(res::RLEV, ct::TransRLEV, ffter::Transformer) = begin
    @inbounds @simd for i = 1 : res.len
        ifftto!(res.stack[i], ct.stack[i], ffter)
    end
end

###########################################################################################

"""
RGSW is a struct for the RGSW cryptosystem. Modulus is fixed by 2⁶⁴.
"""
struct RGSW
    basketb::RLEV
    basketa::RLEV

    function RGSW(basketb::RLEV, basketa::RLEV)
        new(basketb, basketa)
    end
end

"""
TransRGSW is a struct for the Fourier transformed RGSW cryptosystem. 
"""
struct TransRGSW
    basketb::TransRLEV
    basketa::TransRLEV

    function TransRGSW(basketb::TransRLEV, basketa::TransRLEV)
        new(basketb, basketa)
    end
end

const RGSWType = Union{RGSW, TransRGSW}

add(x::T, y::T) where {T<:RGSWType} = 
    T(add(x.basketb, y.basketb), add(x.basketa, y.basketa))

addto!(res::T, x::T, y::T) where {T<:RGSWType} = begin
    addto!(res.basketb, x.basketb, y.basketb)
    addto!(res.basketa, x.basketa, y.basketa)
end

sub(x::T, y::T) where {T<:RGSWType} = 
    T(sub(x.basketb, y.basketb), sub(x.basketa, y.basketa))

subto!(res::T, x::T, y::T) where {T<:RGSWType} = begin
    subto!(res.basketb, x.basketb, y.basketb)
    subto!(res.basketa, x.basketa, y.basketa)
end

initialise!(x::RGSWType) = begin
    initialise!(x.basketb)
    initialise!(x.basketa)
end

fft(ct::RGSW, ffter::Transformer) =
    TransRGSW(fft(ct.basketb, ffter), fft(ct.basketa, ffter))

fftto!(res::TransRGSW, ct::RGSW, ffter::Transformer) = begin
    fftto!(res.basketb, ct.basketb, ffter)
    fftto!(res.basketa, ct.basketa, ffter)
end

ifft(ct::TransRGSW, ffter::Transformer) =
    RGSW(ifft(ct.basketb, ffter), ifft(ct.basketa, ffter))

ifftto!(res::RGSW, ct::TransRGSW, ffter::Transformer) = begin
    ifftto!(res.basketb, ct.basketb, ffter)
    ifftto!(res.basketa, ct.basketa, ffter)
end