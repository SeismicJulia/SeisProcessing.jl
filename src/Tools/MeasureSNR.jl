"""
    MeasureSNR(signal, noisy; db=false)

Measure the signal-to-noise ratio between the clean input `signal` and the
contaminated input `noisy`.

# Arguments
- `signal::Array{Real, N}`: N-dimensional clean signal. `N` must be <= 5.
- `noisy::Array{Real, N}`: N-dimensional noisy signal of same size as `signal`.
- `db::Bool=false`: Flag is false if the signal-to-noise ratio is measured by
amplitude. FLag is true if snr is in dB.

# Example
```julia
julia> d = SeisHypEvents(); dnoisy = SeisAddNoise(d, 2);
MeasureSNR(d, dnoisy)
```
"""
function MeasureSNR(signal::Array{Ts, N},
            noisy::Array{Tn, N}; db::Bool=false) where {Ts <: Real, Tn <: Real, N}

    if db==false
        snr = norm(signal)/norm(signal-noisy)
    elseif db==true
        snr = 20*log10(norm(signal)/norm(signal-noisy))
    end
    return snr

end
