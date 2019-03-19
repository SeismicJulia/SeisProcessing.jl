function SeisMute(in;offset=[0.],tmute=0.,vmute=1500.,taper=0.1,dt=0.001)

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

function SeisMute(in::String,out::String,parameters;group="gather",key=["imx","imy"],ntrace=10000)

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
		itrace_in = 1
		itrace_out = 1
		nx = SeisMain.GetNumTraces(in)
		while itrace_in <= nx
			h = SeisMain.SeisReadHeaders(in,group=group,key=key,itrace=itrace_in,ntrace=ntrace)
			offset = SeisMain.ExtractHeader(h,"h")
			dt = h[1].d1
			num_traces = size(h,1)
			parameters = Dict(:offset=>offset,:tmute=>tmute,:vmute=>vmute,:taper=>taper,:dt=>dt)
			SeisProcessFile(in,out,[SeisMute],[parameters];group=group,key=key,itrace=itrace_in,ntrace=ntrace)
			itrace_in += num_traces
			itrace_out += num_traces
		end
	end
end
