<a name="logo"/>
<div align="center">
<a href="http://saig.physics.ualberta.ca/" target="_blank">
<img src="https://saig.physics.ualberta.ca/lib/tpl/dokuwiki/images/logo.png" alt="SAIG Logo" width="240" height="106"></img>
</a>
</div>

# SeisProcessing.jl

[![Build Status](https://travis-ci.com/SeismicJulia/SeisProcessing.jl.svg?branch=master)](https://travis-ci.com/SeismicJulia/SeisProcessing.jl)

This package contains Processing tools for SeismicJulia project.

At the moment, it is updated and tested against Julia v1.

## Installation

To use this package you must first install the [Julia](http://julialang.org/downloads/) programming language.
Then, run the Julia application and add the package
```
julia>using Pkg
julia>Pkg.add(PackageSpec(url="https://github.com/SeismicJulia/SeisProcessing.jl.git"))
```
or
```
julia>]
pkg> add ""https://github.com/SeismicJulia/SeisProcessing.jl.git"
```

Finally, in the julia prompt type
```
julia>using SeisProcessing
```

## Contents
1. Modelling
SeisLinearEvents, SeisParabEvents, SeisHypEvents, SeisAddNoise, Berlage, Ormsby, Ricker, Hamming

2. Processing
SeisBandPass, SeisDecimate, SeisDelay, SeisDiff, SeisEnvelope, SeisFKFilter, SeisFXDecon, SeisGain, SeisKolmogoroff, SeisMute, SeisNMO, SeisPWD, SeisRadonFreqFor, SeisRadonFreqInv, SeisSemblance, SeisSincInterp1D, SeisSmooth1, SeisSmooth2, SeisSmoothGathers, SeisStack, MeasureSNR, TriangleFilter


## Basic usage
Please, direct to the folder examples. 
Also, the following code produces the figure below.

```Julia
using PyPlot, SeisProcessing
dt = 0.002;
w1 = Ricker(dt=dt);
nw = length(w1);
nc = floor(Int, nw/2);
t1 = dt*collect(-nc:1:nc);

w2 = Ormsby(dt=dt, f=[5.0, 10.0, 30.0, 55.0]);
nw = length(w2);
t2 = dt*collect(0:1:nw-1);

w3 = Berlage(dt=dt);
nw = length(w3);
t3 = dt*collect(0:1:nw-1);

subplot(131)
plot(t1, w1)
axis("tight")
xlabel("Time (s)")
title("Ricker wavelet")

subplot(132)
plot(t2, w2)
axis("tight")
xlabel("Time (s)")
title("Ormsby wavelet")

subplot(133)
plot(t3, w3)
axis("tight")
xlabel("Time (s)")
title("Berlage wavelet")

tight_layout()

```

## For developers: contributing to the project

* New at GitHub? These [basic commands](http://seismic.physics.ualberta.ca/docs/git_basic_commands.pdf)
and this [dictionary](http://seismic.physics.ualberta.ca/docs/git_dictionary.pdf) might help.
* This [tutorial](http://seismic.physics.ualberta.ca/docs/develop_SeismicJulia.pdf) provides the basics
steps you need to follow in order to fork the main repository, change the source code in your forked
repository, commit the changes, and make pull requests using GitHub.
* For contributions to the package, please follow the general guidelines given here:
[Modifications.md](https://github.com/SeismicJulia/Seismic.jl/blob/master/Modifications.md).


## Credits


If you use the SeismicJulia project, please cite the following paper
```
@article{stanton2016efficient,
  title={Efficient geophysical research in Julia},
  author={Stanton, Aaron and Sacchi, Mauricio D},
  journal={CSEG GeoConvention 2016},
  pages={1--3},
  year={2016}
}
```