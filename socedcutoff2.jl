using Arpack, SparseArrays, LinearAlgebra
using JLD
# using FLoops

# define functions used here
include("ades.jl")
include("acre.jl")
include("pascaltriangle.jl")
include("in2b.jl")
# include("b2in.jl")
include("epsilon.jl")
include("cutMsizeEne.jl")
include("Hsocfunccutoff.jl")
include("Hsocfunccutoffk1W1.jl")
include("Hintfunccutoff2.jl")
include("correctionint.jl")

function main2(gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega::Float64, Msize0::Int64, Np::Int64, specnum::Int64)

    # define cutoff of energy in Fock states
    Msize = Msize0*2
    Enecutoff = Msize0 - 1 + Int64(Np/2)
    matp = zeros(Int64,Msize+1,Np+1)
    pascaltriangle!(Msize,Np,matp) # the size is Msize+1 times Np+1
    # note the indices are m+1 and n+1 for N^m_n
    indvec = cutMsizeEne(Msize0,Np,matp,Enecutoff) #Ene0minumhalf
    maxmatpcut = length(indvec)

    # Hamiltonian
    # construct the single-particle Hamiltonian
    mat0 = spzeros(ComplexF64,maxmatpcut,maxmatpcut)
    Hsocfunccutoff!(indvec,Msize0,Np,matp,ksoc,Omega,mat0)

    # construct the interaction Hamiltonian
    mat1 = spzeros(Float64,maxmatpcut,maxmatpcut)
    mat2 = spzeros(Float64,maxmatpcut,maxmatpcut)
    mat3 = spzeros(Float64,maxmatpcut,maxmatpcut)
    Hintfunccutoff2!(indvec,Msize0,Np,matp,mat1,mat2,mat3)

    # diagonalisation
    lambda, phi = eigs(mat0+gdown*mat1+gup*mat2+gdu*mat3,nev=specnum,which=:SR)

    # apply Abel approach (correction of interaction strength)
    if isapprox(ksoc,0)
       vecg = [gdown, gup, gdu]
       vecgind = sortperm(vecg)
       vecg .= vecg[vecgind]
       gcorrected = zeros(Float64,3)
       gcorrected[1] = correctionint(Msize0,vecg[1],real(lambda[1]))
       gcorrected[2] = correctionint(Msize0,vecg[2],real(lambda[2]))
       gcorrected[3] = correctionint(Msize0,vecg[3],real(lambda[3]))
       gcorrected .= gcorrected[vecgind]
       println("ksoc=0")
       return gcorrected, lambda
   else
       println("ksoc~=0")
       return lambda
   end

end

function createHtotal(Msize0::Int64, Np::Int64)
# gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega::Float64,
# specnum::Int64

    # define cutoff of energy in Fock states
    Msize = Msize0*2
    Enecutoff = Msize0 - 1 + Int64(Np/2)
    matp = zeros(Int64,Msize+1,Np+1)
    pascaltriangle!(Msize,Np,matp) # the size is Msize+1 times Np+1
    # note the indices are m+1 and n+1 for N^m_n
    indvec = cutMsizeEne(Msize0,Np,matp,Enecutoff) #Ene0minumhalf
    maxmatpcut = length(indvec)

    # Hamiltonian
    # construct the single-particle Hamiltonian
    # mat0 = spzeros(ComplexF64,maxmatpcut,maxmatpcut)
    matho = spzeros(ComplexF64,maxmatpcut,maxmatpcut)
    matsoc = spzeros(ComplexF64,maxmatpcut,maxmatpcut)
    matW = spzeros(ComplexF64,maxmatpcut,maxmatpcut)
    # Hsocfunccutoff!(indvec,Msize0,Np,matp,ksoc,Omega,mat0)
    Hsocfunccutoffk1W1!(indvec,Msize0,Np,matp,matho,matsoc,matW)

    # construct the interaction Hamiltonian
    mat1 = spzeros(Float64,maxmatpcut,maxmatpcut)
    mat2 = spzeros(Float64,maxmatpcut,maxmatpcut)
    mat3 = spzeros(Float64,maxmatpcut,maxmatpcut)
    Hintfunccutoff2!(indvec,Msize0,Np,matp,mat1,mat2,mat3)

    # diagonalisation
    # lambda, phi = eigs(mat0+gdown*mat1+gup*mat2+gdu*mat3,nev=specnum,which=:SR)

   #  # apply Abel approach (correction of interaction strength)
   #  if isapprox(ksoc,0)
   #     vecg = [gdown, gup, gdu]
   #     vecgind = sortperm(vecg)
   #     vecg .= vecg[vecgind]
   #     gcorrected = zeros(Float64,3)
   #     gcorrected[1] = correctionint(Msize0,vecg[1],real(lambda[1]))
   #     gcorrected[2] = correctionint(Msize0,vecg[2],real(lambda[2]))
   #     gcorrected[3] = correctionint(Msize0,vecg[3],real(lambda[3]))
   #     gcorrected .= gcorrected[vecgind]
   #     println("ksoc=0")
   #     return gcorrected, lambda
   # else
   #     println("ksoc~=0")
   #     return lambda
   # end

   return matho,matsoc,matW,mat1,mat2,mat3

end

function diagonaliseHtotsingle(Msize0::Int64, Np::Int64, gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega::Float64, specnum::Int64)

    matho,matsoc,matW,mat1,mat2,mat3 = createHtotal(Msize0,Np)
    lambda, ~ = eigs(matho+ksoc*matsoc+Omega*matW+gdown*mat1+gup*mat2+gdu*mat3,nev=specnum,which=:SR)

    return lambda

end

function diagonaliseHtotW(Msize0::Int64, Np::Int64, gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega0::Float64, Omega1::Float64, NOmega::Int64, specnum::Int64)

    matho, matsoc, matW, mat1, mat2, mat3 = createHtotal(Msize0,Np)
    save("data_Htot.jld", "matho", matho, "matsoc", matsoc, "matW", matW, "mat1", mat1, "mat2", mat2, "mat3", mat3)

    arrayOmega = LinRange(Omega0,Omega1,NOmega)
    arraylambda = zeros(ComplexF64,NOmega,specnum)

    for jj = 1:NOmega
        arraylambda[jj,:], ~ = eigs(matho+ksoc*matsoc+arrayOmega[jj]*matW+gdown*mat1+gup*mat2+gdu*mat3,nev=specnum,which=:SR)
    end

    # return arraylambda
    save("data_arraylambda.jld", "arrayOmega", arrayOmega, "arraylambda", arraylambda)

    # arraylambda = load("data_arraylambda.jld")["data"]

end
