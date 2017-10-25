module SeisProcessing
    using Interpolations,Requires,Compat
    include("Processing/Processing.jl")
    include("Tools/Tools.jl")
    include("Windows/Windows.jl")
end
