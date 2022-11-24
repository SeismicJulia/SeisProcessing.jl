function SeisFXDecon(d,param=Dict())
	#FX_DECON: SNR enhancement using fx-deconvolution.

	#  [DATA_f] = fx_decon(D1,D2,dt,lf,mu,flow,fhigh);
	# param should be defined, for example
	# param = Dict("dt"=>dt,"filter_length"=>lf,"mu"=>mu,"fmax"=>fhigh)

	#  OUT  de:     filtered data
	#
	#  Reference: Canales, 1984, Random noise reduction, 54.th. Ann. Internat.
	#             Mtg., Soc. Expl. Geophys., Expanded Abstracts, pp. 525-527
	#
	#  Note: Canales' method is modified to use non Toeplitz system of equations
	#        with backward and forward prediction filters
	#



	dt = get(param,"dt",0.002)
	filter_length = get(param,"filter_length",10)
	mu = get(param,"mu",0.01)
	nt = size(d,1)
	nx = size(d,2)
	nf = nt
	dw = 2. *pi/nf/dt
	nw = convert(Int,floor(nf/2)) + 1
	fmax = get(param,"fmax",convert(Int,floor(0.5/dt)))
	if(fmax*dt*nf < nw)
		iw_max = convert(Int,floor(fmax*dt*nf))
	else
		iw_max = convert(Int,floor(0.5/dt))
	end
	D = fft(d,1)
	Df = zeros(ComplexF64,size(D))
	Db = zeros(ComplexF64,size(D))
	for iw=1:iw_max
		y = vec(D[iw,:])
		yf = fxdecon(y,"f",filter_length,mu);
		yb = fxdecon(y,"b",filter_length,mu);
		Df[iw,:] = yf
		Db[iw,:] = yb
	end
	# symmetries
	for iw=nw+1:nf
		Df[iw,:] = conj(Df[nf-iw+2,:])
		Db[iw,:] = conj(Db[nf-iw+2,:])
	end
	dout = ( ifft(Df,1) .+ ifft(Db,1) )/ 2
	out = real(dout[1:nt,:])
        out[:,1:filter_length] = 2*out[:,1:filter_length]
        out[:,nx-filter_length+1:nx] = 2*out[:,nx-filter_length+1:nx]

	return out
end

function fxdecon(d,direction,filter_length,mu)
	nx = length(d);
	if (direction=="b")
		y  = d[1:nx-filter_length];
		C  = d[2:nx-filter_length+1];
		R  = d[nx-filter_length+1:nx];
		M = hankel(C,R);
	else
		y  = d[filter_length+1:nx,1];
		C  = d[filter_length:nx-1];
		R = reverse(d[1:filter_length]);
		M = toeplitz(C,R);
	end
	B = M'*M;
	beta = B[1,1]*mu/100
println(beta)
	ab = (B + beta*Matrix(I,filter_length,filter_length))\M'*y
	temp = M*ab
	if (direction=="b")
		d = [temp ; zeros(filter_length)]
	else
		d = [zeros(filter_length) ; temp]
	end

	return d
end

function hankel(c, r)

	nc = length(r)
	nr = length(c)
	if (nc > 1)
		c = [ c ; r[2:nc] ]
	end
	m = c[ones(Int,nr) * collect(Int,1:nc)' + collect(Int,0:nr-1)*ones(Int,nc)']
	m = reshape(m,nr,nc)

	return m
end

function toeplitz(c, r)

	nc = length(r)
	nr = length(c)
	m = zeros(ComplexF64,nr, nc);
	for i = 1:min(nc, nr)
		m[i:nr,i] .= c[1:nr-i+1]
	end
	for i = 1:min(nr, nc-1)
		m[i,i+1:nc] .= r[2:nc-i+1];
	end

	return m
end
