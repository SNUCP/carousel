using ChaChaCiphers: ChaChaStream, ChaCha12Stream
using FFTW: plan_rfft, plan_brfft, Plan
using Primes: totient, isprime, prevprime, nextprime, factor
using Nemo: finite_field, polynomial_ring, defining_polynomial, ZZ, coeff, lift, zzModPolyRingElem, residue_ring, cyclotomic
using Base.Threads: @threads, @spawn, @sync
using LinearAlgebra: mul! as matmul!

RefBool = Base.RefValue{Bool}
RefUInt32 = Base.RefValue{UInt32}

include("math/modular.jl")
include("math/arithmetic.jl")
include("math/sampler.jl")

include("ring/transformer.jl")
include("ring/ringpoly.jl")
include("ring/encoder.jl")
include("ring/extractor.jl")

include("ciphertext/key.jl")
include("ciphertext/decomposition.jl")
include("ciphertext/ciphertext.jl")
include("ciphertext/encryptor.jl")

include("scheme/evaluator.jl")
include("scheme/params.jl")
include("scheme/bootstrapper.jl")