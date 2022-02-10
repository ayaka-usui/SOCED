using Arpack, SparseArrays, LinearAlgebra
using JLD
# using FLoops

# define functions used here
include("pascaltriangle.jl")
include("cutMsizeEnespinless.jl")
include("cutMsizeEnespinmixed.jl")
include("Hsocfunccutoffk1W1.jl")
include("Hintfunccutoff2.jl")

function createHtotal(Msize0::Int64, Np::Int64)

    # create Fock basis
    # for down down and up up
    Enecutoff = Msize0 - 1 + Np/2
    matp = zeros(Int64,Msize0+1,Np+1)
    pascaltriangle!(Msize0,Np,matp) # note the indices are m+1 and n+1 for N^m_n
    indvec = cutMsizeEnespinless(Msize0,Np,matp,Enecutoff)
    maxmatpcut = length(indvec)

    # for down up
    # Enecutoff = Msize0 - 1 + Np/2
    matp2 = zeros(Int64,Msize0+1,Np)
    pascaltriangle!(Msize0,Np-1,matp2) # note the indices are m+1 and n+1 for N^m_n
    indvec2 = cutMsizeEnespinmixed(Msize0,Np,matp2,Enecutoff)
    maxmatpcut2 = length(indvec2)

    # construct single particle Hamiltonian
    matho = spzeros(Float64,maxmatpcut+maxmatpcut2+maxmatpcut,maxmatpcut+maxmatpcut2+maxmatpcut)
    matsoc = spzeros(ComplexF64,maxmatpcut+maxmatpcut2+maxmatpcut,maxmatpcut+maxmatpcut2+maxmatpcut)
    matW = spzeros(Float64,maxmatpcut+maxmatpcut2+maxmatpcut,maxmatpcut+maxmatpcut2+maxmatpcut)
    Hsocfunccutoffk1W1!(indvec,indvec2,Msize0,Np,matp,matp2,matho,matsoc,matW)

    # construct interaction Hamiltonian
    matdowndown = spzeros(Float64,maxmatpcut+maxmatpcut2+maxmatpcut,maxmatpcut+maxmatpcut2+maxmatpcut)
    matupup = spzeros(Float64,maxmatpcut+maxmatpcut2+maxmatpcut,maxmatpcut+maxmatpcut2+maxmatpcut)
    matdownup = spzeros(Float64,maxmatpcut+maxmatpcut2+maxmatpcut,maxmatpcut+maxmatpcut2+maxmatpcut)
    Hintfunccutoff2!(indvec,indvec2,Msize0,Np,matp,matp2,matdowndown,matupup,matdownup)

    return matho, matdowndown, matupup, matdownup, matsoc, matW

end

function diagonaliseHtotsingle(Msize0::Int64, Np::Int64, gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega::Float64, specnum::Int64)

    matho, matdowndown, matupup, matdownup, matsoc, matW = createHtotal(Msize0,Np)
    lambda, phi = eigs(matho + gdown*matdowndown + gup*matupup + gdu*matdownup + ksoc*matsoc + Omega*matW,nev=specnum,which=:SR)

    return lambda

end

function diagonaliseHtotW(Msize0::Int64, Np::Int64, gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega0::Float64, Omega1::Float64, NOmega::Int64, specnum::Int64)

    matho, matdowndown, matupup, matdownup, matsoc, matW = createHtotal(Msize0,Np)
    # save("data_Htot.jld", "matho", matho, "matsoc", matsoc, "matW", matW, "mat1", mat1, "mat2", mat2, "mat3", mat3)

    arrayOmega = LinRange(Omega0,Omega1,NOmega)
    arraylambda = zeros(ComplexF64,NOmega,specnum)
    arrayspect = zeros(ComplexF64,NOmega,specnum-1)

    for jj = 1:NOmega
        arraylambda[jj,:], ~ = eigs(matho + gdown*matdowndown + gup*matupup + gdu*matdownup + ksoc*matsoc + arrayOmega[jj]*matW,nev=specnum,which=:SR)
        arrayspect[jj,:] .= arraylambda[jj,2:end] .- arraylambda[jj,1]
    end

    return arrayOmega, arraylambda, arrayspect
    # save("data_arraylambda.jld", "arrayOmega", arrayOmega, "arraylambda", arraylambda)

end










function createHtotal0(Msize0::Int64, Np::Int64)
# gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega::Float64,
# specnum::Int64

    println("time to create Fock states")
    @time begin
        # define cutoff of energy in Fock states
        Msize = Msize0*2
        Enecutoff = Msize0 - 1 + Int64(Np/2)
        matp = zeros(Int64,Msize+1,Np+1)
        pascaltriangle!(Msize,Np,matp) # the size is Msize+1 times Np+1
        # note the indices are m+1 and n+1 for N^m_n
        indvec = cutMsizeEne(Msize0,Np,matp,Enecutoff) #Ene0minumhalf
        maxmatpcut = length(indvec)
    end

    println("time to construct single particle Hamiltonian")
    @time begin
        # Hamiltonian
        # construct the single-particle Hamiltonian
        # mat0 = spzeros(ComplexF64,maxmatpcut,maxmatpcut)
        matho = spzeros(Float64,maxmatpcut,maxmatpcut)
        matsoc = spzeros(ComplexF64,maxmatpcut,maxmatpcut)
        matW = spzeros(Float64,maxmatpcut,maxmatpcut)
        # Hsocfunccutoff!(indvec,Msize0,Np,matp,ksoc,Omega,mat0)
        Hsocfunccutoffk1W1!(indvec,Msize0,Np,matp,matho,matsoc,matW)
    end

    println("time to construct interaction Hamiltonian")
    @time begin
        # construct the interaction Hamiltonian
        mat1 = zeros(Float64,maxmatpcut,maxmatpcut)
        mat2 = zeros(Float64,maxmatpcut,maxmatpcut)
        mat3 = zeros(Float64,maxmatpcut,maxmatpcut)
        Hintfunccutoff2!(indvec,Msize0,Np,matp,mat1,mat2,mat3)
    end

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

function diagonaliseHtotsinglewithdata(gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega::Float64, specnum::Int64)

    # matho,matsoc,matW,mat1,mat2,mat3 = createHtotal(Msize0,Np)
    matho = load("data_Htot_M40.jld")["matho"]
    matsoc = load("data_Htot_M40.jld")["matsoc"]
    matW = load("data_Htot_M40.jld")["matW"]
    mat1 =load("data_Htot_M40.jld")["mat1"]
    mat2 =load("data_Htot_M40.jld")["mat2"]
    mat3 =load("data_Htot_M40.jld")["mat3"]

    lambda, ~ = eigs(matho+ksoc*matsoc+Omega*matW+gdown*mat1+gup*mat2+gdu*mat3,nev=specnum,which=:SR)

    return lambda

end

function diagonaliseHtotWwithdata(gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega0::Float64, Omega1::Float64, NOmega::Int64, specnum::Int64)

    # matho, matsoc, matW, mat1, mat2, mat3 = createHtotal(Msize0,Np)
    # save("data_Htot.jld", "matho", matho, "matsoc", matsoc, "matW", matW, "mat1", mat1, "mat2", mat2, "mat3", mat3)

    matho = load("data_Htot_M50.jld")["matho"]
    matsoc = load("data_Htot_M50.jld")["matsoc"]
    matW = load("data_Htot_M50.jld")["matW"]
    mat1 =load("data_Htot_M50.jld")["mat1"]
    mat2 =load("data_Htot_M50.jld")["mat2"]
    mat3 =load("data_Htot_M50.jld")["mat3"]

    arrayOmega = LinRange(Omega0,Omega1,NOmega)
    arraylambda = zeros(ComplexF64,NOmega,specnum)

    for jj = 1:NOmega
        arraylambda[jj,:], ~ = eigs(matho+ksoc*matsoc+arrayOmega[jj]*matW+gdown*mat1+gup*mat2+gdu*mat3,nev=specnum,which=:SR)
    end

    # return arraylambda
    # save("data_arraylambda.jld", "arrayOmega", arrayOmega, "arraylambda", arraylambda)

    return arraylambda, arrayOmega

end
