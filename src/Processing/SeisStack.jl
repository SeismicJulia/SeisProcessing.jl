"""
    stacktraces(in; <keyword arguments>)

Stack an array of traces
# Arguments
- `in`: 2D seismic data
- `normalize=true`: Normalize result by number of traces

```
"""
function stacktraces(d;normalize=true)
	if (normalize == true)
		val = sum(d[:,:],dims=2)/size(d[:,:],2)
	else
		val = sum(d[:,:],dims=2)
	end
	return val
end
