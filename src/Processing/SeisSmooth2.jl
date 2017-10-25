function SeisSmooth2(A;L1=7,L2=7,Nrepeat=1)
    a1 = isodd(L1) ? 0 : 1
    a2 = isodd(L2) ? 0 : 1
    kernel = ones(eltype(A),L1,L2)/(L1*L2);
    for irepeat = 1 : Nrepeat
	A = conv2(A,kernel)
	A = A[floor(Int,L1/2) + 1:end-floor(Int,L1/2) + a1,floor(Int,L2/2) + 1:end-floor(Int,L2/2) + a2]
	end
    return A
end


