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

    spect = real(lambda .- lambda[1])

    return lambda, spect

end

function diagonaliseHtotW(Msize0::Int64, Np::Int64, gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega0::Float64, Omega1::Float64, NOmega::Int64, specnum::Int64)

    println("constructoing the Hamiltonian ...")
    @time begin
        matho, matdowndown, matupup, matdownup, matsoc, matW = createHtotal(Msize0,Np)
        # save("data_Htot.jld", "matho", matho, "matsoc", matsoc, "matW", matW, "mat1", mat1, "mat2", mat2, "mat3", mat3)
    end

    println("diagonalising the Hamiltonian for different Omega ...")
    arrayOmega = LinRange(Omega0,Omega1,NOmega)
    arraylambda = zeros(ComplexF64,NOmega,specnum)
    arrayspect = zeros(ComplexF64,NOmega,specnum-1)
    mat0 = matho + gdown*matdowndown + gup*matupup + gdu*matdownup + ksoc*matsoc

    for jj = 1:NOmega
        @time begin
            arraylambda[jj,:], _ = eigs(mat0 + arrayOmega[jj]*matW,nev=specnum,which=:SR)
            arrayspect[jj,:] .= arraylambda[jj,2:end] .- arraylambda[jj,1]
            println(jj)
        end
    end

    save("data_spectrum.jld", "arrayOmega", arrayOmega, "arraylambda", arraylambda, "arrayspect", arrayspect)

    return arrayOmega, arraylambda, arrayspect

end
