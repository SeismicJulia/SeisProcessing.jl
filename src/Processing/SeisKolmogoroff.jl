"""
    SeisKolmogoroff(in)

Transform a wavelet into its minimum phase equivalent.

# Arguments
- `in::Array{Real,1}`: input wavelet.

# Example
```julia
julia> using PyPlot
julia> w = Ricker()
julia> wmin = SeisKolmogoroff(w)
julia> plot(w); plot(wmin)
```

# Reference
* Claerbout, Jon F., 1976, Fundamentals of geophysical data processing.
McGraw-Hill Inc.
"""
function SeisKolmogoroff(w::Array{T,1}) where T<:Real

    nw = length(w)
    nf = 8*nextpow(2,nw)
    W = fft(cat(w,zeros(nf-nw),dims=1))
    A = log.(abs.(W) .+ 1.e-8)
    a = 2*ifft(A)
    n2 = floor(Int, nf/2)
    a[n2+2:nf] .= 0.0
    a[1] = a[1]/2
    A = exp.(fft(a))
    a = real(ifft(A))
    wmin = real(a[1:nw])

end
