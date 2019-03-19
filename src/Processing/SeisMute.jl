function SeisMute(in;offset=[0.],tmute=0.,vmute=1500.,taper=0.1,dt=0.001)

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

	return out
end



function SeisMute(m::ASCIIString,d::ASCIIString,adj;tmute=0.,vmute=1500.,taper=0.1,ntrace=100000,key=["sx"])

	param["group"] = "some"
	param["ntrace"] = ntrace
	parameters = Dict(:tmute=>tmute,:vmute=>vmute,:taper=>taper)
	if (adj==true)
		SeisProcess(d,m,param)
		SeisProcess(d,m,[SeisMute],[parameters];key=["imx"])
	else
		SeisProcess(m,d,param)
		SeisProcess(m,d,[SeisMute],[parameters];key=["imx"])
	end

end
