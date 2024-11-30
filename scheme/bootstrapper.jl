mutable struct Bootstrapper
    const lweeval::LWEEvaluator; const ringeval::RingEvaluator
    const encoder::SubringEncoder; const extor::SubringExtractor
    const brk::Vector{TransRGSW}; const rtk::Vector{TransRLEV}; const ksk::Vector{LEV}
    const exp_rate::Int64; const mask::TransPoly; const decpar::DecompParams{UInt64}
    const acc::RLWE; const tacc::TransRLWE; 
    const RLWEtmp::Vector{RLWE}; const TransRLWEtmp::Vector{TransRLWE}
    const LWEtmp::LWE; const LWEks::LWE; 
end

function setup(params::Parameters)
    m, p, r = params.m, params.p, params.r

    # Define DecompParams.
    decparLWE = DecompParams{UInt32}(params.ksk_len, params.ksk_logB)
    decparRing = DecompParams{UInt64}(params.brk_len, params.brk_logB)
    decpar = DecompParams{UInt64}(3, 22)

    # Define structs.
    ffter = SubringTransformer(m, p)
    encoder = SubringEncoder(m, p, r)
    extor = SubringExtractor(ffter)
    lweeval = LWEEvaluator(decparLWE)
    ringeval = RingEvaluator(ffter, decparRing)

    # Compute mask.
    tmp = ringeval.adec[1]
    tmp.coeffs[1] = encoder.maskconst
    encode!(tmp, encoder)
    @inbounds @simd for i = eachindex(tmp.coeffs)
        (tmp.coeffs[i] ≥ encoder.pr.Q÷2) && (tmp.coeffs[i] -= encoder.pr.Q)
    end
    mask = fft(tmp, ffter)
    
    # Sample the secret keys.
    rng = ChaCha12Stream()
    lwekey = block_binary_lwekey(rng, params.d, params.ℓ)
    ringkey = partial_ringkey(rng, ffter.N, lwekey)

    # Define Encryptors.
    entorLWE = LWEEncryptor(rng, lwekey, params.α)
    entorRing = RLWEEncryptor(rng, ringkey, params.β, ffter)

    N, n, exp_rate = ffter.N, lwekey.n, params.expansion_rate

    @assert N ≥ n "Choose LWE dimension smaller than ring dimension."

    # Generate Blind Rotation keys.
    brk = Vector{TransRGSW}(undef, n)
    buffRLEV = RLEV([RLWE(N) for _ = 1 : decparRing.len])
    buffRGSW = RGSW(buffRLEV, deepcopy(buffRLEV))
    @inbounds for i = 1 : n
        RGSW_encrypt_to!(buffRGSW, lwekey.key[i] % UInt64, entorRing, decparRing)
        brk[i] = fft(buffRGSW, ffter)
    end

    # Generate key-switching keys.
    rtk = Vector{TransRLEV}(undef, N-1)
    buffkey = deepcopy(ringkey)
    @inbounds for i = 1 : N - 1
        rotate!(buffkey, 1)
        RLEV_encrypt_to!(buffRLEV, buffkey, entorRing, decparRing)
        rtk[i] = fft(buffRLEV, ffter)
        rotate!(rtk[i], -i)
    end

    ksk = [LEV_encrypt(ringkey.coeffs[i] % UInt32, entorLWE, decparLWE) for i = n+1 : N]

    ℓ = params.ℓ
    RLWEtmp, TransRLWEtmp, LWEtmp, LWEks, acc, tacc = 
            [RLWE(N) for _ = 1 : ℓ], [TransRLWE(N) for _ = 1 : ℓ], LWE(N), LWE(n), RLWE(N), TransRLWE(N)

    entorLWE, entorRing, Bootstrapper(lweeval, ringeval, encoder, extor, brk, rtk, ksk, exp_rate, mask, decpar, acc, tacc, RLWEtmp, TransRLWEtmp, LWEtmp, LWEks)
end

function mask!(in::RLWE, boot::Bootstrapper)
    extor, mask, LWEtmp = boot.extor, boot.mask, boot.LWEtmp
    acc, tacc, decpar, decvecb, decveca, buffb, buffa, ffter = 
            boot.acc, boot.tacc, boot.decpar, boot.ringeval.bdec, boot.ringeval.adec, boot.RLWEtmp[1].b, boot.RLWEtmp[1].a, boot.ringeval.ffter
    
    initialise!(acc)
    decompto!(decvecb, in.b, decpar); decompto!(decveca, in.a, decpar)
    @inbounds for i = 1 : 3
        fftto!(tacc.b, decvecb[i], ffter); fftto!(tacc.a, decveca[i], ffter)
        multo!(tacc.b, mask, tacc.b); multo!(tacc.a, mask, tacc.a)
        ifftto!(buffb, tacc.b, ffter); ifftto!(buffa, tacc.a, ffter)
        @. acc.b.coeffs <<= 22; @. acc.a.coeffs <<= 22
        addto!(acc.b, acc.b, buffb); addto!(acc.a, acc.a, buffa)
    end
    copy!(in, acc)
end

function bootstrapto!(res::RLWE, in::RLWE, boot::Bootstrapper)
    # Sample Extraction
    extor, mask, LWEtmp = boot.extor, boot.mask, boot.LWEtmp
    acc, tacc, decpar, decvecb, decveca, buffb, buffa, ffter = 
            boot.acc, boot.tacc, boot.decpar, boot.ringeval.bdec, boot.ringeval.adec, boot.RLWEtmp[1].b, boot.RLWEtmp[1].a, boot.ringeval.ffter

    initialise!(acc)
    decompto!(decvecb, in.b, decpar); decompto!(decveca, in.a, decpar)
    @inbounds for i = 1 : 3
        fftto!(tacc.b, decvecb[i], ffter); fftto!(tacc.a, decveca[i], ffter)
        multo!(tacc.b, mask, tacc.b); multo!(tacc.a, mask, tacc.a)
        ifftto!(buffb, tacc.b, ffter); ifftto!(buffa, tacc.a, ffter)
        @. acc.b.coeffs <<= 22; @. acc.a.coeffs <<= 22
        addto!(acc.b, acc.b, buffb); addto!(acc.a, acc.a, buffa)
    end
    extract!(LWEtmp, acc, extor)

    # Key switching
    lweeval, LWEks, ksk = boot.lweeval, boot.LWEks, boot.ksk
    keyswitchto!(LWEks, LWEtmp, ksk, lweeval)

    # Blind Rotation
    ringeval, RLWEtmp, TransRLWEtmp, brk, rtk = boot.ringeval, boot.RLWEtmp, boot.TransRLWEtmp, boot.brk, boot.rtk
    setacc!(acc, boot.encoder)
    Δ = round(UInt64, 1.8446744073709552e19 / boot.encoder.pr.Q)
    @. acc.b.coeffs *= Δ; @. acc.a.coeffs *= Δ
    blind_rotate!(LWEks, acc, tacc, RLWEtmp, TransRLWEtmp, brk, rtk, ringeval)

    copy!(res.a, acc.a); copy!(res.b, acc.b)
end

function extract!(res::LWE, x::RLWE, extor::SubringExtractor)
    res.b[] = divbits(x.b.coeffs[1], 32) % UInt32
    extractto!(res.a, x.a, extor)
end

function setacc!(acc::RLWE, encoder::SubringEncoder; f::Function = x::UInt64 -> x::UInt64)
    initialise!(acc.a)

    N, pr = acc.b.N, encoder.pr
    @inbounds @simd for i = eachindex(acc.b.coeffs)
        acc.b.coeffs[i] = f(barrett(round(UInt64, pr.Q * i / N), pr))
    end

    encode!(acc.b, encoder)
end

function blind_rotate!(in::LWE, acc::RLWE, tacc::TransRLWE, tmp::Vector{RLWE}, tmpfft::Vector{TransRLWE}, brk::Vector{TransRGSW}, rtk::Vector{TransRLEV}, eval::RingEvaluator)
    N = acc.N
    ℓ = length(tmp); d = in.n ÷ ℓ

    scale = N / (1 << 32)
 
    rotate!(acc.b, round(Int64, in.b[] * scale) % N)
    @inbounds for i = 0 : d - 1
        decomptobuff!(acc, eval)

        # Multiply sᵢ
        initialise!(tacc)
        for j = 1 : ℓ
            idx = round(Int64, in.a[i*ℓ+j] * scale) % N
            if idx > 0
                hoistedexternalprodto!(tmpfft[j], brk[i*ℓ+j], eval)
                subto!(tacc, tacc, tmpfft[j])
            end
        end

        # Rotation
        for j = 1 : ℓ
            idx = round(Int64, in.a[i*ℓ+j] * scale) % N
            if idx > 0
                ifftto!(tmp[j], tmpfft[j], eval)
                rotateto!(tmpfft[j], tmp[j], idx, rtk[idx], eval)
                addto!(tacc, tacc, tmpfft[j])
            end
        end
        
        ifftto!(eval.bdec[1], tacc.b, eval.ffter); ifftto!(eval.adec[1], tacc.a, eval.ffter)
        addto!(acc.b, acc.b, eval.bdec[1]); addto!(acc.a, acc.a, eval.adec[1])
    end
end