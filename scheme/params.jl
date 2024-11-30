struct Parameters
    p::Int64; r::Int64              # plaintext modulus = pʳ
    d::Int64; ℓ::Int64              # LWE dimension = d ⋅ ℓ
    α::Real                         # LWE noise std

    m::Int64                        # base cyclotomic ring degree
    expansion_rate::Int64           # expansion rate for EBS.
    β::Real                         # RLWE noise std

    brk_len::Int64; brk_logB::Int64 # Gadget parameter for blind rotation keys
    ksk_len::Int64; ksk_logB::Int64 # Gadget parameter for key-switching keys
end

params4 = Parameters(
    2, 2, 315, 2,
    249095.878,
    
    87211,
    1, 17121749.034,

    4, 8, 6, 2
)

params4_v2 = Parameters(
    2, 2, 315, 2,
    249095.878,
    
    65537,
    1, 6405.772,

    3, 10, 6, 2
)

params8 = Parameters(
    2, 3, 340, 2,
    100147.879,
    
    87211,
    1, 17121749.034,

    4, 8, 7, 2
)

params8_v2 = Parameters(
    2, 3, 340, 2,
    100147.879,
    
    65537,
    1, 6405.772,

    3, 10, 7, 2
)

params16 = Parameters(
    2, 3, 375, 2,
    27965.788,
    
    174763,
    1, 3.2,

    5, 7, 7, 2
)

params3 = Parameters(
    3, 1, 300, 2,
    430329.792,

    176419,
    1, 3.2,

    4, 8, 6, 2
)

params9 = Parameters(
    3, 2, 350, 2,
    69558.762,

    176419,
    1, 3.2,

    4, 8, 7, 2
)

params5 = Parameters(
    5, 1, 325, 2,
    173012.160,

    38923,
    1, 147047830.338,
    6, 5, 6, 2
)

params5_v2 = Parameters(
    5, 1, 325, 2,
    173012.160,

    221401,
    1, 3.2,
    3, 10, 6, 2
)

params7 = Parameters(
    7, 1, 340, 2,
    100147.879,

    137089,
    1, 88360.643,
    4, 8, 6, 2
)

params11 = Parameters(
    11, 1, 360, 2,
    48312.769,

    83791,
    1, 366.434,
    4, 8, 7, 2
)