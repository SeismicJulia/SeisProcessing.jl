"""
      SeisProcessFile(in,out,operators,parameters;<keyword arguments>)

Run processing flows that read and write from disk

f is a function that has the following syntax: d2 = f(d1,param),
where
param is list of keyword arguments for the function.
Note that f can be a vector of functions. They will be executed
sequentially on the same group of traces.

# Arguments
- `in::String`: input filename
- `out::String`: output filenames
- `operator`
- `parameters`

# Keyword arguments
- `group="gather"` : Options are all, some or gather
- `key=["imx","imy"]` : Defines type of gather
- `itrace=1` : Initial traces
- `ntrace=10000` : Total number of traces to process at once

# Example
Apply a bandpass filter to a seismic cube sequentially, by shot gather.
Assume dt is equal to 0.002.
```
julia> operators = [SeisBandPass]
julia> param = [Dict(:dt=>0.002, :fa=>20,:fb=>30,:fc=>80,:fd=>90)]
julia> SeisProcessFile(filein,fileout,operators,param,key=["sx"])
```

"""
function SeisProcessFile(in::String,out::String,operators,parameters;group="gather",key=["imx","imy"],itrace=1,ntrace=0)

	if (group=="all")
		d1,h1,e1 = SeisMain.SeisRead(in,group=group,key=key,itrace=1,ntrace=ntrace)
		for j = 1 : length(operators)
			op = operators[j]
			d2 = op(d1;parameters[j]...)
			d1 = copy(d2)
		end
		SeisMain.SeisWrite(out,d1,h1,e1)
	else
		#itrace_in = 1
		#itrace_out = 1
        itrace_in = itrace
        itrace_out = itrace
        n_tr = SeisMain.GetNumTraces(in)
        if ntrace == 0
            nx = n_tr
        else
            nx = ntrace
        end
println("nx = ",nx," itrace+nx = ",itrace+nx)
        while itrace_out < itrace+nx && itrace_out <= n_tr
			d1,h1,e1 = SeisMain.SeisRead(in,group=group,key=key,itrace=itrace_in,ntrace=ntrace)
            nt=size(d1,1)
			num_traces_in = size(reshape(d1,nt,:),2)
			for j = 1 : length(operators)
				op = operators[j]
				d2 = op(d1;parameters[j]...)
				d1 = copy(d2)
			end
			num_traces_out = size(d1,2)
			SeisMain.SeisWrite(out,d1,h1,e1,itrace=itrace_out)
			itrace_in += num_traces_in
			itrace_out += num_traces_out
            println("itrace_out = ",itrace_out)

		end
	end

end

function SeisProcess(in::Array{String,1},out::Array{String,1},operators,parameters;group="gather",key=["imx","imy"],ntrace=10000)
	for j = 1 : length(in)
		SeisProcess(in[j],out[j],parameters;group=group,key=key,ntrace=ntrace)
	end

end

function SeisProcess(in1::String,in2::String,out::String,operators,parameters;group="gather",key=["imx","imy"],ntrace=10000)
	group = get(param,"group","some")
	key = get(param,"key",["imx","imy"])
	ntrace = get(param,"ntrace",100)

	if (group=="all")
		d1,h1 = SeisMain.SeisRead(in1,group=group,key=key,itrace=1,ntrace=ntrace)
		d2,h2 = SeisMain.SeisRead(in2,group=group,key=key,itrace=1,ntrace=ntrace)
		for ifunc = 1 : length(operators)
			func = operators[ifunc]
			d3,h3 = func(d1,d2,param)
			d1 = copy(d3)
		end
		SeisMain.SeisWrite(out,d1,h1)
	else
		itrace_in = 1
		itrace_out = 1
		nx = SeisMain.GetNumTraces(in1)
		while itrace_in <= nx
			d1,h1 = SeisMain.SeisRead(in1,group=group,key=key,itrace=itrace_in,ntrace=ntrace)
			d2,h2 = SeisMain.SeisRead(in2,group=group,key=key,itrace=itrace_in,ntrace=ntrace)
            nt = size(d1,1)
			num_traces_in = size(reshape(d1,nt,:),2)
			for ifunc = 1 : length(operators)
				func = operators[ifunc]
				d3,h3 = func(d1,d2,param)
				d1 = copy(d3)
			end
            nt_out = size(d1,1)
			num_traces_out = size(reshape(d1,nt_out,:),2)
			SeisMain.SeisWrite(out,d1,h1,itrace=itrace_out)
			itrace_in += num_traces_in
			itrace_out += num_traces_out
		end

	end

end
