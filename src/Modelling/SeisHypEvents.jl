"""
   SeisHypEvents(; <keyword arguments>)

Generate two dimensional data consisting of hyperbolic events.

# Arguments
- `ot::Real=0.0`: first sample for the time axis in secs.
- `dt::Real=0.004`: time sampling interval in secs.
- `nt::Int=301`: number of time samples.
- `ox::Real=-1000.0`: first sample for spatial dimension in meters.
- `dx::Real=20.0`: sample interval for the spatial dimension in meters.
- `nx::Int=101`: number of samples for the spatial dimension.
- `tau::Vector{Real}=[0.2, 0.6, 0.9]`: intercept traveltimes for each event.
- `vel::Vector{Real}=[1500.0, 2000.0, 3000.0]`: rms velocities in m/s
- `apex::Vector{Real}=[0.0, 0.0, 0.0]`: apex-shifts in meters.
- `amp::Vector{Real}=[1.0, -1.0, 1.0]`: amplitudes for each event.
- `wavelet::AbstractString="ricker"`: wavelet used to model the events.
- `f0::Vector{Real}=[20.0]`: central frequency of wavelet for each event.
# Output
- `d::Array{Real, 2}`: two dimensional data consisting of hyperbolic events.
# Examples
```julia
julia> using SeisPlot
julia> d = SeisHypEvents(); SeisPlotTX(d);
julia> d = SeisHypEvents(apex=[100, 200, -300], f0=[30, 20, 15]);
SeisPlotTX(d);
```
"""
function SeisHypEvents(; ot=0.0, dt=0.004, nt=301,
                        ox1=-1000.0, dx1=20.0, nx1=101,
                        ox2=-1000.0, dx2=20.0, nx2=101,
                        ox3=-1000.0, dx3=20.0, nx3=101,
                        ox4=-1000.0, dx4=20.0, nx4=101,
                        tau::Vector{T1}=[0.2, 0.6, 0.9],
                        vel::Vector{T2}=[1500.0, 2000.0, 3000.0],
                        apex::Vector{T3}=[0.0, 0.0, 0.0],
                        amp::Vector{T4}=[1.0, -1.0, 1.0],
                        wavelet="ricker",
                        f0::Vector{T5}=[20.0]) where {T1<:Real, T2<:Real, T3<:Real, T4<:Real, T5<:Real}

    x1 = ox1 .+ collect(0:1:nx1-1)*dx1;
    x2 = ox2 .+ collect(0:1:nx2-1)*dx2;
    x3 = ox3 .+ collect(0:1:nx3-1)*dx3;
    x4 = ox4 .+ collect(0:1:nx4-1)*dx4;

    if length(f0) != length(tau)
	f0 = f0[1] .+ 0.0*tau
    end

    nf = 4*nextpow(2,nt)
    dw = 2.0*pi/(nf*dt)
    nw = floor(Int, floor(nf/2)) + 1
    d = zeros(Float64, nt, nx1, nx2, nx3, nx4);
    D = zeros(Complex{Float64}, nf, nx1,nx2,nx3,nx4)

    nevents = length(tau)

    for ievent = 1:nevents
        if wavelet=="ricker"
            w = Ricker(f0=f0[ievent], dt=dt)
            nw = length(w)
            w = cat(w, zeros(nf-nw),dims=1)
            W = fft(w)
            delay = dt*(floor(Int, nw/2))
            shift = sqrt.(tau[ievent]^2 .+ ((x1 .- apex[ievent])/vel[ievent]).^2)'
        else
            error("Only the Ricker wavelet is implemented at the moment")
        end
        for ix4 = 1:nx4
         for ix3 = 1:nx3
          for ix2 = 1:nx2
           for ix1 = 1:nx1
            for iw = 1:nw
                wrs = (iw - 1)*dw
                shift = sqrt.(tau[ievent]^2 .+ ((x1[ix1]-apex[ievent])/vel[ievent]).^2
                                             + ((x2[ix2]-apex[ievent])/vel[ievent]).^2
                                             + ((x3[ix3]-apex[ievent])/vel[ievent]).^2
                                             + ((x4[ix4]-apex[ievent])/vel[ievent]).^2)'

                D[iw:iw, :] += amp[ievent]*W[iw]*exp.(-1im*wrs*(shift .- delay))
            end
          end
        end
     end
    end
    for iw = nw+1:nf
        D[iw, :] = conj(D[nf-iw+2, :])
    end
    d = ifft(D, 1)
    d = real(d[1:nt, :,:,:,:])
    # extent = Extent(Int32(nt), Int32(nx), Int32(1), Int32(1), Int32(1),
    #                 Float32(ot), Float32(ox), Float32(0), Float32(0),
    #                 Float32(0), Float32(dt), Float32(dx), Float32(1),
    #                 Float32(1), Float32(1), "Time", "Offset", "ix2", "ix3",
    #                 "ix4", "s", "m", "index", "index", "index",
    #                 "Hyperbolic events")
    # return d, extent
    return d
end
