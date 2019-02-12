"""
    SeisGain(d, t; <keyword arguments>)

Gain a group of traces. Input and output are 2D.

# Arguments
- `d::Array{Real,2}`: two dimensional data.

# Keyword arguments
- `dt::Real=0.004`: sampling interval in secs.
- `kind::AbstractString="time"`: if kind="time", gain = t.^a . * exp(-bt);
                         if kind="agc", automatic gain control is applied.
- `param::Vector{Real}=[2.0,0.0]`: if kind="time", param = [a,b];
                                   if kind="agc", param = [agc_gate]
- `norm::Int=0`: `norm=0` no normalization; `norm=1` normalize each trace by
                  amplitude; `norm=2` normalize each trace by rms value/


# Example
```julia
julia> using PyPlot
julia> d = SeisHypEvents();
       dout = SeisGain(d, kind="agc", param=[0.05]);
       SeisPlotTX([d dout]);
```

Credits: Juan I. Sabbione, Aaron Staton, Mauricio D. Sacchi, 2016

"""
function SeisGain(d::Array{Td,2}; dt::Real=0.004, kind::AbstractString="time",
                    param::Vector{Tp}=[2.0,0.0], norm::Int=0) where {Td<:Real,Tp<:Real}

    nt = size(d,1)
    nx = size(d[:,:],2)
    dout = zero(d)

    if kind == "time"   # Geometrical spreading-like gain

        a,b = param
        t = collect(0:1:nt-1)*dt

        tgain = [(t[i]^a)*exp(b*t[i]) for i in 1:nt]

        dout = [d[i,k]*tgain[i] for k in 1:nx, i in 1:nt]

    end

    if kind=="agc"   #AGC

        L = floor(Int,round(param[1]/(2dt)))
        h = triang(2L+1)

        for k = 1:nx
            e = [d[i,k]^2 for i in 1:nt]
            c = conv(e,h)[L+1:nt+L]
            nc = length(c)
            rms = [sqrt(abs(c[i])) for i in 1:nc]
            epsi = 1.e-10*maximum(rms)
            op = [rms[i]/(rms[i]^2 + epsi) for i in 1:nc]
            dout[:,k] = d[:,k] .* op
        end
    end

    if norm==1     # Normalize by amplitude

        for k = 1:nx
            amax = maximum([abs(d[i,k]) for i in 1:nt])
            dout[:,k] = dout[:,k] ./ amax
        end

    end

    if norm==2;    # Normalize by rms

        for k = 1:nx
            amax = [sqrt(sum(d[i,k]^2)/nt) for i in 1:nt];
            dout[:,k] = dout[:,k] ./ amax;
        end

    end

    return dout

end

function triang(n::Integer)
    convert(Array{Float32,1},[1 - abs((k - (n-1)/2))/(n/2) for k in 0:(n-1)])
end
