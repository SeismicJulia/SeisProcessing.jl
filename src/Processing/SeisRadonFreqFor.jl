"""
    SeisRadonFreqFor(in, nt; <keyword arguments>)

Transform a tau-p gather to a time-offset gather using a frequency domain
forward parabolic or linear Radon operator.

# Arguments
- `in::Array{Float64,2}`: 2D Radon panel, `in[1:ntau,1:np]`, where `ntau` is the
number of intercept times and `np` the number of curvatures or ray parameters.
- `nt::Int`: number of time samples.

# Keyword arguments
- `order="parab"`: `"parab"` for parabolic transform, `"linear"`
for linear transform.
- `dt=0.004`: time sampling interval in seconds.
- `h=collect(0.0:20.0:1000.0)`: offset vector; `h[1:nh]`.
- `href=0.0`: reference offset for parabolic Radon Transform. If the
defautl value `href=0.0` is given, `href` is set to `max(abs(h))`.
- `p=collect(-0.05:0.01:2.2)`: `p[1:np]`. If `order="parab"`, `p`
is a vector of residual moveout ("curvatures") at reference offset `href` in
seconds; if `order=linear`, `p` is a vector of ray parameters in s/m.
- `flow=0.0`: minimum frequency in the data in Hz.
- `fhigh=125.0`: maximum frequency in the data in Hz.

# Output
- `d`: data synthetized via forward Radon modeling, `d[1:nt, 1:nh]`.

# References
*  Hampson, D., 1986, Inverse velocity stacking for multiple elimination:
Canadian Journal of Exploration Geophysics, 22, 44-55.
*  Sacchi, M.D. and Ulrych, T.J., 1995, High-resolution velocity gathers and
offset space reconstruction: Geophysics, 60, 1169-1177.
"""
function SeisRadonFreqFor(m::Array{Tm,2}, nt::Int; order="parab",
                            dt=0.004, href=0.0,
                            h=collect(0.0:20.0:1000.0),
                            p=collect(-0.05:0.01:2.2),
                            flow=0.0, fhigh=125.0)where{Tm<:Real}

    if order=="parab"
        ind = 2
        href == 0 && (href = maximum(abs.(h)))
    elseif order=="linear"
        ind = 1
        href = 1.0
    else
        error("Order should be equal to \"parab\" or \"linear\"")
    end
    ntau = size(m, 1)
    np = length(p)
    np == size(m, 2) || error("lenght(p) must be equal to size(m, 2)")
    nh = length(h)
    nw = 2*nextpow(2,ntau)
    m = cat( m, zeros(nw-ntau, np),dims=1)
    M = fft(m, 1)
    iw_low = round(Int, flow*dt*nw+1)
    iw_high = round(Int, fhigh*dt*nw+1)
    D = zeros(ComplexF64, nw, nh)
    for iw = iw_low:iw_high
        w  = 2.0*pi*(iw-1)/(nw*dt)
        L = zeros(ComplexF64, nh, np)
        for ip = 1:np
            for ih = 1:nh
                phi = w*p[ip]*(h[ih]/href)^ind
      	        L[ih, ip] = exp(-im*phi)
            end
        end
        D[iw, :] = L*M[iw,:]
    end
    for iw = round(Int, nw/2)+2:nw
        D[iw, :] = conj(D[nw-iw+2, :])
    end
    d = real(ifft(D, 1))
    d = d[1:nt, :]
    return d
end
