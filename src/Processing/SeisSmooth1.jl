function SeisSmooth1(A;L=7)
	a = isodd(L) ? 0 : 1
	kernel = ones(eltype(A),L)/L;
	A = conv(A,kernel)
	A = A[floor(Int,L/2) + 1:end-floor(Int,L/2) + a]
	return A
end


