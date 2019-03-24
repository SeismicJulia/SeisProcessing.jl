"""
    SeisBandPass(in; <keyword arguments>)

Apply a bandpass filter to seismic data. Input and output data is 2D in tx domain. Filter is applied in fx domain.

# Arguments
- `in`: Input 2D data array in tx domain. Time is first dimension.

# Keyword arguments
- `dt=0.004`: time sampling interval in secs.
- `fa=0`: lower frequency cut [Hz]
- `fb=0`: lower frequency banpass [Hz]
- `fc=60`: upper frequency bandpass [Hz]
- `fd=80`: upper frequency cut [Hz]


# Example
```
julia> d = SeisLinearEvents(); SeisPlotAmplitude(d,125,0.004);
julia> d_filter = SeisBandPass(d;dt=0.004,fa=2,fb=8,fc=12,fd=20); SeisPlotAmplitude(d_filter,125,0.004)
julia> SeisPlot([d d_filter],title="Data and Filtered data")

```
"""
function SeisBandPass(d;dt=0.004,fa=0,fb=0,fc=60,fd=80)

#println(size(d))
	v = size(d)
    nt = v[1]
    #println(nt)
	dn = reshape(d,nt,:)
    println(size(d))
	nx = size(dn,2)
    #println(nx," ",typeof(nx))
	nf = iseven(nt) ? nt : nt + 1
	df = 1/nf/dt
	nw = round(Int,nf/2) + 1

	if(fd*dt*nf < nw)
		iw_max = round(Int,floor(fd*dt*nf))
	else
		iw_max = round(Int,floor(0.5/dt))
	end

	dp = pad_first_axis(dn,nf)
	m = fft(dp,1)/sqrt(size(dp,1))
	if fa > 0.
		m[1,:] *= 0.
	end
	for iw=2:iw_max
		f = df*(iw-1)
		if (f<fa)
			m[iw,:] *= 0.
		elseif (f >= fa && f < fb)
			m[iw,:] *= (f-fa)/(fb-fa)
		elseif (f >= fb && f <= fc)
			m[iw,:] *= 1.
		elseif (f > fc && f <= fd)
			m[iw,:] *= 1. - (f-fc)/(fd-fc)
		else
			m[iw,:] *= 0.
		end
	end
	m[iw_max:end,:] .= 0.

	# symmetries
	for iw=nw+1:nf
		m[iw,:] = conj(m[nf-iw+2,:])
	end
	dn = real(bfft(m,1)/sqrt(size(m,1)))
    println(nx
	dout = dn[1:nt,1:nx]
    #println(size(dout), " ",size(d))
    prinltn(v)
	return reshape(dout,v);
end


"""
    SeisBandPass(in,out,parameters; <keyword arguments>)

# Arguments
- `in::String`: Input file - Seis format
- `out::String`: Output file - Seis format.
- `parameters` : list of the keyword arguments for the function SeisBandPass.

# Keyword arguments
- `group="gather"` : Options are all, some or gather
- `key=["imx","imy"]` : Defines type of gather
- `itrace=1` : Initial trace number
- `ntrace=10000` : Total number of traces to process at once

# Example
```
julia> param = Dict(:fa>=0,:fb=>0,:fc=>60,:fd=>80)
julia> SeisBandPass(filein,fileout, param,group="all")
```

"""
function SeisBandPass(in::String,out::String,parameters;group="gather",key=["imx","imy"],itrace=1,ntrace=10000)

    fa = get(parameters,:fa,0)
    fb = get(parameters,:fb,0)
    fc = get(parameters,:fc,60)
    fd = get(parameters,:fd,80)

    ext = SeisMain.ReadTextHeader(in);
    dt = ext.d1
    parameters = Dict(:dt=>dt,:fa=>fa,:fb=>fb,:fc=>fc,:fd=>fd )


	SeisProcessFile(in,out,[SeisBandPass],[parameters];group=group,key=key,itrace=itrace,ntrace=ntrace)
end



function pad_first_axis(a,N1)
	n1 = size(a,1)
	nx = size(a,2)
	b = zeros(N1,nx)
	for ix = 1 : nx
		b[1:n1,ix] = a[:,ix]
	end
	return b
end
