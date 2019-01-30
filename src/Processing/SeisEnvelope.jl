"""
    SeisEnvelope(in; <keyword arguments>)

Calculate the envelope of input seismic array.

# Arguments
* `in`: input seismic traces

# Output
Array of envelopes

# Example
```julia
julia> using SeisPlot
julia> d = SeisLinearEvents(); SeisPlotTX(SeisEnvelope(d));
```
*Credits: Aaron Stanton,2017*

"""
function SeisEnvelope(d)

	D = fft(d,1)
	D[1:Int(floor(size(d,1)/2)),:] .= 0.0
	return 2*abs.(ifft(D,1))

end
