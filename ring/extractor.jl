"""
An Extractor type for the subring.
"""
mutable struct SubringExtractor
    const M::Array{UInt32, 2}
    const buff::Vector{UInt32}

    function SubringExtractor(ffter::SubringTransformer)
        N = ffter.N

        ffta = fft(vcat(one(UInt64), zeros(UInt64, N-1)), ffter)
        fftb = similar(ffta)
        buff = Vector{UInt64}(undef, N)

        M = zeros(UInt32, N, N)
        @inbounds for i = 1 : N
            @. fftb = ffta
            circshift!(fftb, 1-i)
            @. fftb *= ffta
            ifftto!(buff, fftb, ffter)
            @simd for j = eachindex(buff)
                buff[j] != 0 && (M[j==1 ? 1 : N-j+2, (N-j+i)%N+1] = buff[j] % UInt32)
            end
        end

        new(M, buff)
    end
end

"""
extract! performs sample extract for input polynomial a, and stores the value in res.
"""
function extractto!(res::Vector{UInt32}, a::RingPoly, extor::SubringExtractor)
    buff, M = extor.buff, extor.M

    @inbounds @simd for i = eachindex(a.coeffs)
        buff[i] = divbits(a.coeffs[i], 32) % UInt32
    end
    
    matmul!(res, M, buff)
end