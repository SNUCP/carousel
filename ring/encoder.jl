"""
resolution(m, p, r) returns resolutions of cyclotomic polynomial of degree m over Z_pʳ.
"""
function resolution(m::Int64, p::Int64, r::Int64)
    d = ord(p, m)
    N = (m - 1) ÷ d

    Kp, _ = finite_field(p, d, " ")
    ZZy, y = ZZ["y"]
    ZZpr, _ = residue_ring(ZZ, p^r)
    ZZprz, z = polynomial_ring(ZZpr, "z")

    poly = ZZprz(lift(ZZy, defining_polynomial(Kp)))
    Kpr, _ = residue_ring(ZZprz, poly)
    Kprx, x = Kpr["x"]

    # Order of the multiplicative group of the multiplicative ring.
    order = big(p)^((r - 1) * d) * (big(p)^d - 1)

    # Find the m-th root of unity over the finite ring Kpr.
    #TODO Make it deterministic instead of random.
    ξ = rand(Kpr)
    while true
        if ξ^order == Kpr(1)
            ζ = ξ^(order ÷ m)
            ζ ≠ Kpr(1) && break
        end
        ξ = rand(Kpr)
    end
    ζ = ξ^(order ÷ m)

    # Compute the first factor.
    factx = Kprx(1)
    for j = 0:d-1
        factx *= x - ζ^(powermod(p, j, m))
    end

    # Embed into the integer coefficient polynomials.
    factz = ZZprz(0)
    for j = reverse(0:d)
        factz = factz * z + coeff(coeff(factx, j).data, 0).data
    end

    # Compute its compliment.
    cyclo = ZZprz(cyclotomic(m, y))
    compz = cyclo ÷ factz

    # Compute the resolution of unity.
    R, _ = residue_ring(ZZprz, factz)
    τ = (R(compz)^(order - 1)).data * compz

    # Extract the meaningful values.
    g = primitive_root_finder(m)
    resol = Vector{Int64}(undef, N)
    pr, idx = p^r, 1
    tmp = pr - coeff(τ, 0).data
    for i = 1:N
        resol[i] = coeff(τ, idx).data + tmp
        resol[i] ≥ pr && (resol[i] -= pr)
        idx = idx * g % m
    end

    # Use an adequate resolution to make the extraction part easier.
    idx = findfirst(x -> gcd(x, p) == 1, resol)
    circshift!(resol, -idx + 1)

    resol
end

find_and_save_resolution(m::Int64, p::Int64, r::Int64) = begin
    resol = resolution(m, p, r)
    save_resolution("m$(m)p$(p)r$(r)", resol)
    resol
end

save_resolution(paramname::String, resol::Vector{Int64}) = begin
    open(String(@__DIR__) * "/resolutions.jl", "a") do file
        println(file, paramname, " = ", resol)
    end
end

function load_resolution(m::Int64, p::Int64, r::Int64)
    resol = open(String(@__DIR__) * "/resolutions.jl") do file
        for str in eachline(file)
            if occursin("m$(m)p$(p)r$(r)", str)
                idx1 = last(findfirst("[", str))
                idx2 = first(findfirst("]", str))
                resol = parse.(Int64, split(str[idx1+1:idx2-1], ","))
                return resol
            end
        end
    end

    if isnothing(resol)
        @info "No resolution of unity for given parameters found. Resolution will be found and saved in resolution.jl."
        find_and_save_resolution(m, p, r)
    else
        resol
    end
end

"""
`SubringEncoder` is a struct which support the conversion between the η-basis and ξ-basis.
You can generate a SubringEncoder over a ring Z_pʳ[X]/Φₘ(X) invariant by X→Xᵖ by `SubringEncoder(m, p, r)`.
Currently, we only cover the case where d = ord(p, m) and N = (m-1)/d are even. 
"""
mutable struct SubringEncoder
    p::Int64
    r::Int64
    m::Int64
    N::Int64
    pr::Modulus
    ffter::CyclicTransformer
    buff1::Vector{ComplexF64}
    buff2::Vector{Float64}
    maskconst::UInt64
    resol::Vector{ComplexF64}
    invresol::Vector{ComplexF64}

    function SubringEncoder(m::Int64, p::Int64, r::Int64)
        d = ord(p, m)
        N = (m - 1) ÷ d

        # Define CyclicTransformer and Modulus.
        ffter = CyclicTransformer(N)
        pr = Modulus(p^r)

        # Load resolution of unity and compute inverse of it.
        resol = UInt64.(load_resolution(m, p, r))
        invresol = barrett.(m * resol .+ d, Ref(pr))

        # Compute mask constant, which is used during the sample extraction.
        maskconst = invmod(resol[1], pr.Q)

        # Compute the FFT transform of the resolution and inverse resolution.
        resol, invresol = real_fft(Float64.(resol), ffter) / N, real_fft(Float64.(invresol), ffter) / N

        # Define the buffers.
        buff1, buff2 = Vector{ComplexF64}(undef, 1 + N >> 1), Vector{Float64}(undef, N)

        new(p, r, m, N, pr, ffter, buff1, buff2, maskconst, resol, invresol)
    end
end

"""
`encode!(a, encoder)` computes an in-place encoding of `a`.
"""
encode!(a::RingPoly, encoder::SubringEncoder) = encodeto!(a.coeffs, a.coeffs, encoder)

"""
`encodeto!(res, a, encoder)` computes the encoding of `a` and stores it in `res`.
"""
@views function encodeto!(res::Vector{UInt64}, a::Vector{UInt64}, encoder::SubringEncoder)
    N, pr = encoder.N, encoder.pr
    ffter, buff1, buff2, resol = encoder.ffter, encoder.buff1, encoder.buff2, encoder.resol

    @assert length(a) == N "The length of input vector does not match the parameter."

    buff2[1] = barrett(signed(a[1]), pr)
    buff2[2:end] .= barrett.(signed.(a[end:-1:2]), Ref(pr))

    real_fftto!(buff1, buff2, ffter)
    @. buff1 *= resol
    real_ifftto!(buff2, buff1, ffter)

    res .= barrett.(native.(buff2), Ref(pr))
end

"""
`decode!(a, encoder)` does an in-place decoding of `a`.
"""
decode!(a::Vector{UInt64}, encoder::SubringEncoder) = decodeto!(a, a, encoder)

"""
`decodeto!(res, a, encoder)` computes the decoding of `a` and stores it in `res`.
"""
@views function decodeto!(res::Vector{UInt64}, a::Vector{UInt64}, encoder::SubringEncoder)
    N, pr = encoder.N, encoder.pr
    ffter, buff1, buff2, invresol = encoder.ffter, encoder.buff1, encoder.buff2, encoder.invresol

    @assert length(a) == N "The length of input vector does not match the parameter."

    buff2[1] = barrett(signed(a[1]), pr)
    buff2[2:end] .= barrett.(signed.(a[end:-1:2]), Ref(pr))

    real_fftto!(buff1, buff2, ffter)
    @. buff1 *= invresol
    real_ifftto!(buff2, buff1, ffter)

    res .= barrett.(native.(buff2), Ref(pr))
end