"""
An Extractor type for the subring.
"""
mutable struct SubringExtractor
    const m::Int64
    const d::UInt32
    const idx::Vector{Int64}
    const buff::Vector{UInt32}

    function SubringExtractor(ffter::SubringTransformer)
        m, N = ffter.m, ffter.N
        d, g = (m - 1) ÷ N, primitive_root_finder(m)

        idx = Vector{Int64}(undef, m - 1)
        @inbounds for i = 0:m-1
            idx[powermod(g, i, m)] = (i % N) + 1
        end

        buff = Vector{UInt32}(undef, N)

        new(m, d, idx, buff)
    end
end

"""
extract! performs sample extract for input polynomial a, and stores the value in res.
"""
function extractto!(res::Vector{UInt32}, a::RingPoly, extor::SubringExtractor)
    m, d, idx, buff = extor.m, extor.d, extor.idx, extor.buff

    @inbounds @simd for i = eachindex(a.coeffs)
        buff[i] = divbits(a.coeffs[i], 32) % UInt32
    end

    @. res = -buff * d
    @inbounds for i = 2:m-1
        res[idx[i]] += buff[idx[m-i+1]]
    end
end