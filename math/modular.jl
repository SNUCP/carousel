"""
Modulus supports a modular arithmetic given modulus Q.
"""
struct Modulus
    Q::UInt64
    r::UInt128
    
    function Modulus(Q::Integer)
        logQ = ceil(Int64, log2(Q))
        @assert logQ ≤ 31 "Modulus should be smaller than 2³¹."

        new(Q, floor(UInt128, (big(1)<<64)/Q))
    end
end

"""
barrett(x, Q) returns x % Q using barret reduction.
"""
barrett(x::UInt64, Q::Modulus) = begin
    t = (Q.r * x) >> 64
    res = (x - Q.Q * t) % UInt64
    res ≥ Q.Q && (res -= Q.Q)
    res
end

barrett(x::Int64, Q::Modulus) = begin
    if x ≥ 0
        res = barrett(UInt64(x), Q)
    else
        res = barrett(UInt64(-x), Q)
        res = res ≠ zero(UInt64) ? Q.Q - res : res
    end 
    
    res
end