"""
DecompParams is a struct for digit decomposition parameters.
"""
struct DecompParams{T}
    len::Int64; mask::T
    logB::Int64; halfB::T
    gveclog::Vector{Int64}; gvec::Vector{T}

    function DecompParams{T}(len::Int64, logB::Int64) where {T<:Unsigned}
        B = one(T) << logB
        gveclog = trailing_ones(typemax(T)) .- collect(1:len) * logB
        gveclog[end] < 0 && (gveclog .-= gveclog[end])
        gvec = one(T) .<< gveclog
        mask = B - 1
        new(len, mask, logB, B >> 1, gveclog, gvec)
    end
end

"""
`decomp(a, params)` returns a gadget decomposition of the 32-bit integer `a`.
"""
decomp(a::UInt32, params::DecompParams{UInt32}) = begin
    res = Vector{UInt32}(undef, params.len)
    decompto!(res, a, params)
    res
end

"""
`decompto!(avec, a, params)` computes a gadget decomposition of the 32-bit integer `a` and stores it in `avec`.
"""
decompto!(avec::Vector{UInt32}, a::UInt32, params::DecompParams{UInt32}) = begin
    ai = divbits(a, params.gveclog[end])
    @inbounds for i = params.len : -1 : 2
        avec[i] = ai & params.mask
        ai >>= params.logB
        ai += avec[i] >> (params.logB - 1)
        avec[i] -= (avec[i] & params.halfB) << 1
    end
    avec[1] = ai & params.mask
    avec[1] -= (avec[1] & params.halfB) << 1
end

"""
`decomp(a, params)` returns a gadget decomposition of the integer polynomial `a`.
"""
decomp(a::RingPoly, params::DecompParams{UInt64}) = begin
    res = [RingPoly(a.N) for _ = 1 : params.len]
    decompto!(res, a, params)
    res
end

"""
`decompto!(avec, a, params)` computes a gadget decomposition of the integer polynomial `a` and stores it in `avec`.
"""
decompto!(avec::Vector{RingPoly}, a::RingPoly, params::DecompParams{UInt64}) = begin
    @. avec[1].coeffs = a.coeffs >> params.gveclog[end]
    @inbounds for j = params.len : -1 : 2
        @. avec[j].coeffs = avec[1].coeffs & params.mask
        @. avec[1].coeffs >>= params.logB
        @. avec[1].coeffs += avec[j].coeffs >> (params.logB - 1)
        @. avec[j].coeffs -= (avec[j].coeffs & params.halfB) << 1
    end
    @. avec[1].coeffs &= params.mask
    @. avec[1].coeffs -= (avec[1].coeffs & params.halfB) << 1
end