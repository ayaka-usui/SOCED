using Arpack, SparseArrays, LinearAlgebra
using JLD
# using FLoops

# define functions used here
# include("ades.jl")
# include("acre.jl")
include("pascaltriangle.jl")
# include("in2b.jl")
# include("b2in.jl")
# include("epsilon2.jl")
include("cutMsizeEne2.jl")
include("Hsocfunccutoff_spinless.jl")
include("Hintfunccutoff2_spinless.jl")
include("correctionint.jl")

function main2_spinless(g0::Float64, Msize::Int64, Np::Int64, specnum::Int64)

    # define cutoff of energy in Fock states
    Enecutoff = Msize - 1 + Np/2
    matp = zeros(Int64,Msize+1,Np+1)
    pascaltriangle!(Msize,Np,matp) # the size is Msize+1 times Np+1
    # note the indices are m+1 and n+1 for N^m_n
    indvec = cutMsizeEne2(Msize,Np,matp,Enecutoff) #Ene0minumhalf
    maxmatpcut = length(indvec)

    # Hamiltonian
    # construct the single-particle Hamiltonian
    mat0 = spzeros(Float64,maxmatpcut,maxmatpcut)
    Hsocfunccutoff_spinless!(indvec,Msize,Np,matp,mat0)

    # construct the interaction Hamiltonian
    mat1 = spzeros(Float64,maxmatpcut,maxmatpcut)
    Hintfunccutoff2_spinless!(indvec,Msize,Np,matp,mat1)

    # diagonalisation
    # lambda, phi = eigs(mat0+g0*mat1,nev=specnum,which=:SR)
    # lambda0 = eigvals(Array(mat0+g0*mat1))

    # return lambda
    return mat0+g0*mat1

end
