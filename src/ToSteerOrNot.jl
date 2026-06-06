module ToSteerOrNot

using LinearAlgebra, Statistics, Random  # std lib.
export mean

using Optim, Plots
export Optim, LBFGS, Adam, ConjugateGradient

import ForwardDiff, DifferentiableEigen

import IterativeSolvers, LinearMaps

include("optim.jl")
include("utils.jl")

include("extra.jl")  # misc funcitons in files I didn't include here!

include("two_whole.jl")
include("better.jl")
include("trajectories.jl")

# include("threedims.jl")  # my simple 3D code
# include("three_wigner.jl")  # new fancy 3D code

end # module ToSteerOrNot
