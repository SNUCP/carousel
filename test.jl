include("Carousel.jl")

function main()
    # Parameter Selection. You can choose from the parameters in scheme/params.jl
    param = params4

    _, entor, boot = setup(param)
    pr = boot.encoder.pr.Q

    Δ = round(UInt64, 2.0^64 / pr)
    oneoverΔ = pr / 2.0^64

    message = rand(0:pr-1) % UInt64
    input = RLWE_encrypt(message * Δ, entor)

    for _ = 1:10
        @time bootstrapto!(input, input, boot)
        res = phase(input, entor)

        for i = eachindex(res.coeffs)
            res.coeffs[i] = round(UInt64, res.coeffs[i] * oneoverΔ) % pr
        end

        decode!(res.coeffs, boot.encoder)
        @assert (Int.(res.coeffs[1]) == message)
    end
end

main()