"""
    SeisFKFilter(in; <keyword arguments>)

Removes energy from seismic data by applying filters in the frequency-wavenumber domain


# Arguments
* `in`: input seismic array in TX domain

# Keyword arguments
* `dt=0.004`: time sampling interval in secs.
* `dx=10`: space sampling interval in meters.
* `va=-2000,vb=-3000,vc=3000,vd=2000`: apparent velocity range to filter

# Output
Filtered data in TX domain

# Example
```julia
julia> d = SeisLinearEvents(); df = SeisFKFilter(d,va=-4000,vb=-6000,vc=6000,vd=4000);
```
*Credits: Aaron Stanton,2017*

"""
function SeisFKFilter(d;dt=0.004,dx=10,va=-2000,vb=-3000,vc=3000,vd=2000)

    nt = size(d,1)
    nx = size(d,2)
    pa = 1. /va
    pb = 1. /vb
    pc = 1. /vc
    pd = 1. /vd
    m = fft(d,1)
    m = fft(m,2)
    nf = size(m,1)
    dw = 1. /nf/dt
    nw = isodd(nf) ? round(Int,floor(nf/2)) + 1 : round(Int,floor(nf/2))
    nk = size(m,2)
    dk = 1. /nk/dx
    for iw=1:nw
	w = dw*(iw-1)
	for ik=1:nk
	    k = ik < round(Int,floor(nk/2)) + 1 ? dk*(ik-1) : -(dk*nk - dk*(ik-1))
	    p = k/w;
	    if p < pa
		m[iw,ik] *= 0.
	    elseif p >= pa && p < pb
		m[iw,ik] *= (p-pa)/(pb-pa)
	    elseif p >= pb && p <= pc
				m[iw,ik] *= 1.
	    elseif p > pc && p <= pd
		m[iw,ik] *= 1. - (p-pc)/(pd-pc)
	    elseif p > pd
				m[iw,ik] *= 0.
	    end
	end
    end
    # symmetries
    m = ifft(m,2)
    for iw=nw+1:nf
	m[iw,:] = conj(m[nf-iw+2,:])
    end
    d = real(ifft(m,1))

    return d[1:nt,1:nx];
end
