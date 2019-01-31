using SeisProcessing
using Test

nt=250;
dt=0.004;
nx=30;
dx=5;

param = Dict(:nt=>nt, :dt=>dt,:nx1=>nx,:dx1=>dx,
              :p1=>[0.0003, -0.0005, 0.000], :p2=>[0,0,0],:p3=>[0,0,0],:p4=>[0,0,0],
              :tau=>[0.3, 0.4,0.5], :amp=>[ -1.0, 0.7,1],:f0=>20);


d=SeisLinearEvents(;param...);

h = collect(0:nx-1)*dx;


# 2- Define linear moveouts p
pmin = -2*dt/dx;
pmax = 2*dt/dx;
np = 60
dp = (pmax-pmin)/(np-1)
p = collect(pmin:dp:pmax);

# 3- Transform to Radon domain

m = SeisRadonFreqInv(d, order="linear", dt=dt, h=h, p=p, flow=2, fhigh=80,
                     mu=0.00001 )

# 4- Transform back to data domain
d2 = SeisRadonFreqFor(m, nt, order="linear", dt=dt, h=h, p=p, flow=2, fhigh=80)

# 5- Calculate and print modeling error (%)
e2 = 100*(sum((d-d2).^2))/sum(d.^2);
println("Relative error = ", e2, " %")
@test e2 < 0.1
