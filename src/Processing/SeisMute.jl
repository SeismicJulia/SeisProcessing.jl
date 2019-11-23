"""
    SeisMute(d; <keyword arguments>) -> Array{Real,2}

Mute trace noise before the first arrival.
The muting function depends on the offset between source and receiver and the velocity of the first layer as
t_mute = t_0 + offset/v1.


# Arguments
- `d`: Array of input traces.

# Keyword arguments
- `offset=[0.]`: Vector of distances between source and receiver.
- `tmute=0.`: initial muting time
- `vmute=1500.`: Velocity of the first layer.
- `taper=0.1`: Taper of the muting function
- `dt::Real=0.004`: sampling interval in secs.

"""
function SeisMute(in;offset=[0.],tmute=0.,vmute=1500.,taper=0.1,dt=0.004)

	vec= size(in)
	in = reshape(in,vec[1],:)
	nt,nx = size(in)
	out = copy(in)
	for it = 1:nt
		for ix = 1:nx
			t = sqrt(tmute^2 + (offset[ix]/vmute).^2)
			if ((it-1)*dt <= t )
					if ((it-1)*dt > t - taper)
						out[it,ix] *= 1 - (t - (it-1)*dt)/taper
					else
						out[it,ix] = 0.
					end
			end
		end
	end
	out = reshape(out,vec)
	return out
end


"""
    SeisMute(in,out,parameters; <keyword arguments>)

Mutes noise before first arrival of traces in a file. Saves the muted traces to an output file.
The muting function depends on the offset between source and receiver and the velocity of the first layer as
t_mute = t_0 + offset/v1.


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

	tmute = get(parameters,:tmute,0.)
	vmute = get(parameters,:vmute,1500.)
	taper = get(parameters,:taper,0.1)

	if (group=="all")
		headers = SeisMain.SeisReadHeaders(in);
		offset = SeisMain.ExtractHeader(headers,"h")
		dt = headers[1].d1
		parameters = Dict(:offset=>offset,:tmute=>tmute,:vmute=>vmute,:taper=>taper,:dt=>dt)

		SeisProcessFile(in,out,[SeisMute],[parameters];group=group)
	else
		itrace_in = itrace #1
		#itrace_out = 1
		nx = SeisMain.GetNumTraces(in)
		while itrace_in <= nx
			h = SeisMain.SeisReadHeaders(in,group=group,key=key,itrace=itrace_in,ntrace=ntrace)
			offset = SeisMain.ExtractHeader(h,"h")
			dt = h[1].d1
			#num_traces = size(h,1)
			parameters = Dict(:offset=>offset,:tmute=>tmute,:vmute=>vmute,:taper=>taper,:dt=>dt)
			#SeisProcessFile(in,out,[SeisMute],[parameters];group=group,key=key,itrace=itrace_in,ntrace=ntrace)
			d1,h1,e1 = SeisMain.SeisRead(in,group=group,key=key,itrace=itrace_in,ntrace=ntrace)

			d2 = SeisMute(d1;parameters...)
			d1 = copy(d2)

			SeisMain.SeisWrite(out,d1,h1,e1,itrace=itrace_in)

			#println(itrace_in," ",ntrace)

			itrace_in += ntrace

			#itrace_in += num_traces
			#itrace_out += num_traces
		end
	end
end
