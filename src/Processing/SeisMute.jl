function SeisMute(in;offset=[0.],tmute=0.,vmute=1500.,taper=0.1,dt=0.001)

	vec= size(d)
	in = reshape(d,vec[1],:)
	nx = size(in,2)
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

	headers = SeisMain.SeisReadHeaders(in);
	offset = SeisMain.ExtractHeader(headers,"h")
	dt = headers[1].d1
	parameters = Dict(:offset=>offset,:tmute=>tmute,:vmute=>vmute,:taper=>taper,:dt=>dt)
	#if (adj==true)
		SeisProcess(d,m,[SeisMute],[parameters];group=group,key=key,ntrace=ntrace)
	#else
	#	SeisProcess(m,d,[SeisMute],[parameters];key=["imx"])
	#end

end
