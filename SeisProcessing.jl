module SeisProcessing
    using Interpolations,Requires,Compat
    include("Processing/Processing.jl")
    include("Reconstruction/Reconstruction.jl")
end
