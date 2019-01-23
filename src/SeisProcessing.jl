module SeisProcessing
    using Interpolations
    using Requires,
    FFTW,
    LinearAlgebra,
    DSP
    include("Processing/Processing.jl")
    include("Modelling/Modelling.jl")
    include("Tools/Tools.jl")
end
