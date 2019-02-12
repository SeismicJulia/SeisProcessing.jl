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
	out = zeros(size(in))
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
