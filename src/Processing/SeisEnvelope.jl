function SeisEnvelope(d)

	D = fft(d,1)
	D[1:Int(floor(size(d,1)/2)),:] .= 0.0
	return 2*abs.(ifft(D,1))

end
