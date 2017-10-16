function SeisSmoothGathers(in, out,adj;Nsmooth=3,Nrepeat=1)

	if (adj==false)
		ngather = size(in,3)
		for igather = 1 : ngather
			A = in[:,:,igather];
			A = smooth_angles(A;Nsmooth=Nsmooth,Nrepeat=Nrepeat)
			out[:,:,igather] = A
		end
	else
		ngather = size(out,3)
		for igather = 1 : ngather
			A = out[:,:,igather];
			A = smooth_angles(A;Nsmooth=Nsmooth,Nrepeat=Nrepeat)
			in[:,:,igather] = A
		end
	end
return in,out
end

function smooth_angles(A;Nsmooth=3,Nrepeat=1)

	for iz = 1 : size(A,1)
		a = A[iz,:]
		for irepeat = 1 : Nrepeat
			a = mean_filter(a,Nsmooth)
		end
		A[iz,:] = a
	end

	return A;
end

function mean_filter(a,nw)

	b = vec(zeros(length(a),1))
	for ix = 1 : length(a)
		sum  = 0.
		nsum = 0.
		for iw = 1 : nw
			index1 = ix - Int(floor(nw/2)) - 1 + iw
			if (index1>0 && index1<length(a))
				sum += a[index1]
				nsum += 1.
			end
		end
		b[ix] = sum/nsum
	end
	return b
end
