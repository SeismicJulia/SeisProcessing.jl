function SeisDiff(d;dt=0.001,pow=-2,rot=0)

	nt = size(d,1)
	D = fft(d,1)
	dw = 2. * pi/nt/dt
	nw = convert(Int,floor(nt/2)) + 1
	eps = pow < 0 ? 10.0 : 0.0
	for iw=1:nw
		D[iw,:] *= exp(rot*1im)*((iw*dw) + eps)^pow
	end
	# symmetries
	for iw=nw+1:nt
		D[iw,:] = conj(D[nt-iw+2,:])
	end
	d = real(ifft(D,1))
	return d

end
