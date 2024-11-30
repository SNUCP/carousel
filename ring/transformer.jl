"""
`CyclicTransformer` is a struct which support the fast Fourier transform (FFT) over a ring R[X]/(Xᴺ-1).
You can generate a CyclicTransformer over a ring R[X]/(Xᴺ-1) by `CyclicTransformer(N)`.
Currently, we only cover the case where N is even.
"""
mutable struct CyclicTransformer
    const N::Int64
    const plan::Plan; const iplan::Plan

    function CyclicTransformer(N::Int64)
        # Define plans.
        plan  = plan_rfft(Vector{Float64}(undef, N))
        iplan = plan_brfft(Vector{ComplexF64}(undef, N>>1+1), N)

        new(N, plan, iplan)
    end
end

"""
`real_fft(a, ffter)` outputs discrete Fourier transform of a real vector `a`.
"""
function real_fft(a::Vector{Float64}, ffter::CyclicTransformer)
    res = Vector{ComplexF64}(undef, ffter.N >> 1 + 1)
    real_fftto!(res, a, ffter)
    res
end

"""
`real_fftto!(res, a, ffter)` computes discrete Fourier transform of a real vector `a` and stores it in the array `res`.
"""
real_fftto!(res::Vector{ComplexF64}, a::Vector{Float64}, ffter::CyclicTransformer) =
    matmul!(res, ffter.plan, a)

"""
`real_ifft(a, ffter)` outputs the coefficient form of the discrete Fourier transform `a`.
"""
function real_ifft(a::Vector{ComplexF64}, ffter::CyclicTransformer)
    res = Vector{Float64}(undef, ffter.N)
    real_ifftto!(res, a, ffter)
    res
end

"""
`real_ifftto!(res, a, ffter)` computes the coefficient form of the discrete Fourier transform `a` and stores it in the array `res`.
"""
real_ifftto!(res::Vector{Float64}, a::Vector{ComplexF64}, ffter::CyclicTransformer) =
    matmul!(res, ffter.iplan, a)

"""
`SubringTransformer` is a struct which support the fast Fourier transform (FFT) over the subring.
You can generate a SubringTransformer over a ring Z[X]/Φₘ(X) invariant by X→Xᵖ by `SubringTransformer(m, p)`.
Currently, we only cover the case where d = ord(p, m) is even. 
By doing so, we can make the result of Subring FFT real.
"""
mutable struct SubringTransformer
    const m::Int64; const N::Int64
    const Ψ::Vector{ComplexF64}; const Ψinv::Vector{ComplexF64}
    const buff1::Vector{ComplexF64}; const buff2::Vector{Float64}
    const ffter::CyclicTransformer

    function SubringTransformer(m::Int64, p::Int64)
        @assert isprime(m) "$m is not a prime number."

        g, d = primitive_root_finder(m), ord(p, m); N = (m-1) ÷ d

        @assert iseven(d) "d = ord(p, m) is not even."

        Ñ = factor(Vector, N) ⊆ [2, 3, 5] ? N : nextprod((2, 3, 5), 2N)

        # Define buffers. 
        buff1 = Vector{ComplexF64}(undef, Ñ>>1+1)
        buff2 = Vector{Float64}(undef, Ñ)

        # Generate the resolution of unity.
        Ψ = zeros(Float64, Ñ)
        mmod = Modulus(m)
        idx = g
        @inbounds for _ = 0 : d-1, i = 1 : N
            Ψ[i] += cos(-2*π*idx/m)
            idx = barrett(idx * g, mmod)
        end

        # Inverse resolution of unity is the same, 
        # but we scale them to make FFT more efficient.
        Ψinv = zeros(Float64, Ñ)
        @. Ψinv[1:N] = ((@view Ψ[1:N]) - d) / m

        # Define CyclicTransformer.
        ffter = CyclicTransformer(Ñ)

        # Pre-compute the FFT form of resolution of unities.
        Ψ, Ψinv = real_fft(Ψ, ffter) / Ñ, real_fft(Ψinv, ffter) / Ñ
        
        new(m, N, Ψ, Ψinv, buff1, buff2, ffter)
    end
end

"""
`fft(a, ffter)` returns the subring FFT of the UInt64 vector `a`.
"""
function fft(a::Vector{UInt64}, ffter::SubringTransformer)
    res = Vector{Float64}(undef, ffter.N)
    fftto!(res, a, ffter)
    res
end

"""
`fftto!(res, a, ffter)` computes the subring FFT of the UInt64 vector `a` and stores it in `res`.
"""
@views function fftto!(res::Vector{Float64}, a::Vector{UInt64}, ffter::SubringTransformer)
    N, Ψ = ffter.N, ffter.Ψ
    buff1, buff2 = ffter.buff1, ffter.buff2
    
    if N == length(buff2)
        buff2[1] = signed(a[1])
        @. buff2[2:end] = signed(a[end:-1:2])
        
        real_fftto!(buff1, buff2, ffter.ffter)
        @. buff1 *= Ψ
        real_ifftto!(res, buff1, ffter.ffter)
    else
        buff2[1] = signed(a[1])
        @. buff2[2:N] = signed(a[end:-1:2])
        @. buff2[N+1:end] = 0
    
        real_fftto!(buff1, buff2, ffter.ffter)
        @. buff1 *= Ψ
        real_ifftto!(buff2, buff1, ffter.ffter)
        @. res = buff2[1:N] + buff2[N+1:2N]
    end
end

"""
`ifft(a, ffter)` returns the subring iFFT of the Float64 vector `a`.
"""
function ifft(a::Vector{Float64}, ffter::SubringTransformer)
    res = Vector{UInt64}(undef, ffter.N)
    ifftto!(res, a, ffter)
    res
end

"""
`ifftto!(res, a, ffter)` computes the subring iFFT of the real vector `a` and stores it in `res`.
"""
@views function ifftto!(res::Vector{UInt64}, a::Vector{Float64}, ffter::SubringTransformer)
    N, Ψinv = ffter.N, ffter.Ψinv
    buff1, buff2 = ffter.buff1, ffter.buff2
    
    if N == length(buff2)
        buff2[1] = a[1]
        @. buff2[2:end] = a[end:-1:2]

        real_fftto!(buff1, buff2, ffter.ffter)
        @. buff1 *= Ψinv
        real_ifftto!(buff2, buff1, ffter.ffter)
        @. res = native(buff2)
    else
        buff2[1] = a[1]
        @. buff2[2:N] = a[end:-1:2]
        @. buff2[N+1:end] = 0

        real_fftto!(buff1, buff2, ffter.ffter)
        @. buff1 *= Ψinv
        real_ifftto!(buff2, buff1, ffter.ffter)
        @. res = native(buff2[1:N] + buff2[N+1:2N])
    end
end

const Transformer = Union{CyclicTransformer, SubringTransformer}