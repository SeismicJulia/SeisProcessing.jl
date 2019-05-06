"""
    SeisDecimate(in; <keyword arguments>)

Apply random or regular decimation to a multidimensional array of data. Input and output have the same dimensions.

# Arguments
- `in`: input data as 2D, or 3D,4D,5D tensors. The first dimension is time.

# Keyword arguments
- `mode="random"`: decimation mode. Random is default, else decimates uniformly.
- `perc=50`: percentage of traces equal to zero (decimated). Only for random mode
- `incx1=1`: consecutive traces zeroed in first dimension. Only for regular decimation
- `incx2=1`: consecutive traces zeroed in second dimension.
- `incx3=1`: consecutive traces zeroed in third dimension.
- `incx4=1`: consecutive traces zeroed in fourth dimension.


# Example
```julia
julia> d = SeisLinearEvents(); deci = SeisDecimate(d);
```

*Credits: Aaron Stanton,2017*

"""
function SeisDecimate(in;mode="random",perc=50,incx1=1,incx2=1,incx3=1,incx4=1)

    if (mode=="random") # decimate data randomly
        out = copy(in)
        out = reshape(out,size(in,1),:)
	mask = rand(1,size(in,2)*size(in,3)*size(in,4)*size(in,5));
	mask[(LinearIndices(mask .< perc/100))[findall(mask .< perc/100)]] .= 0;
	mask[(LinearIndices(mask .>= perc/100))[findall(mask .>= perc/100)] ] .= 1;
	for it = 1 : size(in,1)
	    out[it,:] = out[it:it,:].*mask;
	end
    out = reshape(out,size(in));
    else # decimate data regularly with respect to 4 spatial dimensions
        out = zero(in)
	if (size(in,5)) > 1
	    out[:,1:incx1:end,1:incx2:end,1:incx3:end,1:incx4:end] = in[:,1:incx1:end,1:incx2:end,1:incx3:end,1:incx4:end]
	elseif (size(in,4)) > 1
            out[:,1:incx1:end,1:incx2:end,1:incx3:end] = in[:,1:incx1:end,1:incx2:end,1:incx3:end]
        elseif (size(in,3)) > 1
            out[:,1:incx1:end,1:incx2:end] = in[:,1:incx1:end,1:incx2:end]
        elseif (size(in,2)) > 1
            out[:,1:incx1:end] = in[:,1:incx1:end]
        end

    end

    return out
end
