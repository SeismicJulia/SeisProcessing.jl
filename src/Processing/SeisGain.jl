"""
    SeisGain(d ; <keyword arguments>)

Gain a group of traces. Input and output are 2D.

# Arguments
- `d::Array{Real,2}`: two dimensional data.

# Keyword arguments
- `dt::Real=0.004`: sampling interval in secs.
- `kind::AbstractString="time"`: if kind="time", gain = t.^a . * exp(-bt);
                         if kind="agc", automatic gain control is applied.
- `coef::Vector{Real}=[2.0,0.0]`: if kind="time", coef = [a,b];
                                   if kind="agc", coef = [agc_gate]
- `normal::Int=0`: `normal=0` no normalization; `normal=1` normalize each trace by
                  amplitude; `normal=2` normalize each trace by rms value/


# Example
```julia
julia> using PyPlot
julia> d = SeisHypEvents();
       dout = SeisGain(d, kind="agc", coef=[0.05]);
       SeisPlotTX([d dout]);
```

Credits: Juan I. Sabbione, Aaron Staton, Mauricio D. Sacchi, 2016

"""
function SeisGain(d::Array{Td,2}; dt::Real=0.004, kind::AbstractString="time",
                    coef::Vector{Tp}=[2.0,0.0], normal::Int=0) where {Td<:Real,Tp<:Real}

    nt = size(d,1)
    nx = size(d,2)
    dout = zero(d)

    if kind == "time"   # Geometrical spreading-like gain

        a,b = coef
        t = collect(0:1:nt-1)*dt

        tgain = [(t[i]^a)*exp(b*t[i]) for i in 1:nt]
        dout = [d[i,k]*tgain[i] for i in 1:nt, k in 1:nx]

    end

    if kind=="agc"   #AGC

        L = floor(Int,round(coef[1]/(2dt)))
        h = triang(2L+1)

        for k = 1:nx
            e = [d[i,k]^2 for i in 1:nt]
            c = conv(promote(e,h)...)[L+1:nt+L]
            nc = length(c)
            rms = [sqrt(abs(c[i])) for i in 1:nc]
            epsi = 1.e-10*maximum(rms)
            op = [rms[i]/(rms[i]^2 + epsi) for i in 1:nc]
            dout[:,k] = d[:,k] .* op
        end
    end

    if normal==1     # Normalize by amplitude

        for k = 1:nx
            amax = maximum([abs(d[i,k]) for i in 1:nt])
            dout[:,k] = dout[:,k] ./ amax
        end

    end

    if normal==2;    # Normalize by rms

        for k = 1:nx
            amax = [sqrt(sum(d[i,k]^2)/nt) for i in 1:nt];
            dout[:,k] = dout[:,k] ./ amax;
        end

    end

    return dout

end

"""
    SeisGain(in,out,parameters; <keyword arguments>)

    Gain a group of traces. Input and output are file names.

# Arguments
- `in::String`: Input file - Seis format
- `out::String`: Output file - Seis format.
- `parameters` : list of the keyword arguments for the function SeisGain.

# Keyword arguments
- `group="gather"` : Options are all, some or gather
- `key=["imx","imy"]` : Defines type of gather
- `itrace=1` : Initial trace number
- `ntrace=10000` : Total number of traces to process at once

# Example
```
julia> param = Dict(:dt=>0.01,:kind=>"time",:coef=>[2.0,0.0],:normal=>0)
julia> SeisGain(filein,fileout, param,group="all")
```

"""
function SeisGain(in::String,out::String,parameters;group="gather",key=["imx","imy"],itrace=1,ntrace=10000)

	kind = get(parameters,:kind,"time")
	coef = get(parameters,:coef,[2.0,0.0])
	normal = get(parameters,:normal,0)

    headers = SeisMain.SeisReadHeaders(in);
    dt = headers[1].d1
    parameters = Dict(:dt=>dt,:kind=>kind,:coef=>coef,:normal=>normal)


	SeisProcessFile(in,out,[SeisGain],[parameters];group=group,key=key,itrace=itrace,ntrace=ntrace)
end



function triang(n::Integer)
    convert(Array{Float32,1},[1 - abs((k - (n-1)/2))/(n/2) for k in 0:(n-1)])
end
