function SeisNMO(in;dt=0.001,offset=1000.,tnmo=[0.],vnmo=[1500.],max_stretch=1000)

	nt,nx = size(in)
	if length(offset) < size(in,2)
		offset = offset[1]*fill!(similar(in[1,:]),one(eltype(in[1,:])))
	end


	# interpolate tau,v pairs to match sampling of input data
	if (length(vnmo) == 1)
		tnmo = convert(Float64,tnmo[1]);
		vnmo = convert(Float64,vnmo[1]);
		ti = collect(0:1:nt-1)*dt
		vi = ones(1,nt)*vnmo

	else
		tnmo = (convert(Array{Float64,1},vec(tnmo)),)
		vnmo = convert(Array{Float64,1},vec(vnmo))
		ti = collect(0:1:nt-1)*dt
#g = InterpIrregular(tnmo, vnmo, BCnan, InterpLinear)
		g = interpolate(tnmo, vnmo, Gridded(Linear()))
		ge = extrapolate(g, Line())
		vi = ge(ti)
	end
	out = zeros(typeof(in[1,1]),size(in))
	M = zeros(nt,1)
	for it = 1:nt
		for ix = 1:nx
			time = sqrt(ti[it]^2 + (offset[ix]/vi[it]).^2)
			stretch = (time-ti[it])/(ti[it]+1e-10)
			if (stretch<max_stretch/100)
				its = round(Int,time/dt)+1
				it1 = round(Int,floor(time/dt))+1
				it2 = it1+1
				a = its-it1
				if (it2 <= nt)
					out[it,ix] = (1-a)*in[it1,ix]+a*in[it2,ix]
				end
			end
		end
	end
	return out
end

"""
    SeisNMO(in,out,parameters; <keyword arguments>)

# Arguments
- `in::String`: Input file - Seis format.
- `out::String`: Output file - Seis format.
- `parameters` : list of the keyword arguments for the function SeisMute.

# Keyword arguments
- `group="gather"` : Options are all, some or gather
- `key=["imx","imy"]` : Defines type of gather
- `itrace=1` : Initial trace number
- `ntrace=10000` : Total number of traces to process at once

# Example
```
julia> param = Dict(:tmute=>0.0, :vmute=>10000, :taper=>0.05,:dt=>0.01)
julia> SeisMute(filein,fileout, param,group="some")
```

"""
function SeisMute(in::String,out::String,parameters;group="gather",key=["imx","imy"],itrace=1,ntrace=10000)

	tnmo = get(parameters,:tnmo,[0.])
	vnmo = get(parameters,:vnmo,[1500.])
	max_stretch = get(parameters,:max_stretch,1000)

	if (group=="all")
		headers = SeisMain.SeisReadHeaders(in);
		offset = SeisMain.ExtractHeader(headers,"h")
		dt = headers[1].d1
		parameters = Dict(:offset=>offset,:tnmo=>tmute,:vnmo=>vmute,:max_stretch=>max_stretch,:dt=>dt)

		SeisProcessFile(in,out,[SeisNMO],[parameters];group=group)
	else
		itrace_in = 1
		itrace_out = 1
		nx = SeisMain.GetNumTraces(in)
		while itrace_in <= nx
			h = SeisMain.SeisReadHeaders(in,group=group,key=key,itrace=itrace_in,ntrace=ntrace)
			offset = SeisMain.ExtractHeader(h,"h")
			dt = h[1].d1
			num_traces = size(h,1)
			parameters = Dict(:offset=>offset,:tnmo=>tmute,:vnmo=>vmute,:max_stretch=>taper,:dt=>dt)
			SeisProcessFile(in,out,[SeisMute],[parameters];group=group,key=key,itrace=itrace_in,ntrace=ntrace)
			itrace_in += num_traces
			itrace_out += num_traces
		end
	end
end
