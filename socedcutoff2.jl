using Arpack, SparseArrays, LinearAlgebra
using JLD
using Combinatorics
# using Plots
# using ArnoldiMethod
# using FLoops

# define functions used here
include("pascaltriangle.jl")
include("cutMsizeEnespinless.jl")
include("cutMsizeEnespinmixed.jl")
include("cutMsizeEnespinmixed2.jl")
include("Hsocfunccutoffk1W1.jl")
include("Hintfunccutoff2.jl")
include("coefficientonebodysummary.jl")
include("paircorrelation2.jl")
# include("paircorrelation2_test.jl")
include("changebasis.jl")

function test_size(Msize0::Int64, Np::Int64)

    if Np != 3
       error("This code is specific for Np=3.")
    end

    # create Fock basis
    # for down down and up up
    Enecutoff = Msize0 - 1 + Np/2
    # Enecutoff = 13.5
    matp = zeros(Int64,Msize0+1,Np+1)
    pascaltriangle!(Msize0,Np,matp) # note the indices are m+1 and n+1 for N^m_n
    indvec = cutMsizeEnespinless(Msize0,Np,matp,Enecutoff)
    maxmatpcut = length(indvec)

    # for down down up
    # Enecutoff2 = Msize0 - 1 + (Np-1)/2
    # matp2 = zeros(Int64,Msize0+1,Np-1)
    matp20 = zeros(Int64,Msize0+1,Np-1+1) # Np-1=2
    matp21 = zeros(Int64,Msize0+1,1+1)
    pascaltriangle!(Msize0,Np-1,matp20)
    pascaltriangle!(Msize0,1,matp21)
    indvec2 = cutMsizeEnespinmixed(Msize0,Np,matp20,matp21,Enecutoff,1)
    # indvec2 = cutMsizeEnespinmixed2(Msize0,Np,matp20,matp21,Enecutoff2,1)
    # indvec3 = cutMsizeEnespinmixed(Msize0,Np,matp20,matp21,Enecutoff,1)
    maxmatpcut2 = length(indvec2)
    # maxmatpcut3 = length(indvec3)

    return (maxmatpcut+maxmatpcut2)*2

end

function createHtotal(Msize0::Int64, Np::Int64)

    if Np != 3
       error("This code is specific for Np=3.")
    end

    # create Fock basis
    # for down down and up up
    Enecutoff = Msize0 - 1 + Np/2
    # Enecutoff = 13.5
    matp = zeros(Int64,Msize0+1,Np+1)
    pascaltriangle!(Msize0,Np,matp) # note the indices are m+1 and n+1 for N^m_n
    indvec = cutMsizeEnespinless(Msize0,Np,matp,Enecutoff)
    maxmatpcut = length(indvec)

    # for down down up
    # Enecutoff2 = Msize0 - 1 + (Np-1)/2
    # matp2 = zeros(Int64,Msize0+1,Np-1)
    matp20 = zeros(Int64,Msize0+1,Np-1+1) # Np-1=2
    matp21 = zeros(Int64,Msize0+1,1+1)
    pascaltriangle!(Msize0,Np-1,matp20)
    pascaltriangle!(Msize0,1,matp21)
    indvec2 = cutMsizeEnespinmixed(Msize0,Np,matp20,matp21,Enecutoff,1)
    # indvec2 = cutMsizeEnespinmixed2(Msize0,Np,matp20,matp21,Enecutoff2,1)
    # indvec3 = cutMsizeEnespinmixed(Msize0,Np,matp20,matp21,Enecutoff,1)
    maxmatpcut2 = length(indvec2)
    # maxmatpcut3 = length(indvec3)
    
    # for up up down
    # matp30 = zeros(Int64,Msize0+1,Np-2+1) # Np-2=1
    # matp31 = zeros(Int64,Msize0+1,2+1)
    # pascaltriangle!(Msize0,Np-2,matp30)
    # pascaltriangle!(Msize0,2,matp31)
    # indvec3 = cutMsizeEnespinmixed(Msize0,Np,matp30,matp31,Enecutoff,2)
    # maxmatpcut3 = length(indvec3) # maxmatpcut3 == maxmatpcut2

    # construct transformation from 2nd quantisation basis to spin basis (only for downdownup and upupdown)
    mat_from2ndto1st_downup, mat_from1sttospin_downdownup, mat_from1sttospin_upupdown = changefrom2ndtospin_downup(indvec2,Msize0,Np,matp20,matp21)
    return mat_from2ndto1st_downup, mat_from1sttospin_downdownup, mat_from1sttospin_upupdown

    # construct single particle Hamiltonian
    matho = spzeros(Float64,maxmatpcut+maxmatpcut2*2+maxmatpcut,maxmatpcut+maxmatpcut2*2+maxmatpcut)
    matsoc = spzeros(Float64,maxmatpcut+maxmatpcut2*2+maxmatpcut,maxmatpcut+maxmatpcut2*2+maxmatpcut)
    matW = spzeros(Float64,maxmatpcut+maxmatpcut2*2+maxmatpcut,maxmatpcut+maxmatpcut2*2+maxmatpcut)
    Hsocfunccutoffk1W1!(indvec,indvec2,Msize0,Np,matp,matp20,matp21,matho,matsoc,matW)

    # construct interaction Hamiltonian
    matdowndown = spzeros(Float64,maxmatpcut+maxmatpcut2*2+maxmatpcut,maxmatpcut+maxmatpcut2*2+maxmatpcut)
    matupup = spzeros(Float64,maxmatpcut+maxmatpcut2*2+maxmatpcut,maxmatpcut+maxmatpcut2*2+maxmatpcut)
    matdownup = spzeros(Float64,maxmatpcut+maxmatpcut2*2+maxmatpcut,maxmatpcut+maxmatpcut2*2+maxmatpcut)
    Hintfunccutoff2!(indvec,indvec2,Msize0,Np,matp,matp20,matp21,matdowndown,matupup,matdownup)

    return matho, matdowndown, matupup, matdownup, matsoc, matW

end

function mat_from2ndtospin_downup(Msize0::Int64, Np::Int64)

    # Since we care about population for spin basis, we do not touch down down down and up up up

    if Np != 3
       error("This code is specific for Np=3.")
    end

    # create Fock basis
    # for down down and up up
    Enecutoff = Msize0 - 1 + Np/2
    # Enecutoff = 13.5
    matp = zeros(Int64,Msize0+1,Np+1)
    pascaltriangle!(Msize0,Np,matp) # note the indices are m+1 and n+1 for N^m_n
    indvec = cutMsizeEnespinless(Msize0,Np,matp,Enecutoff)
    maxmatpcut = length(indvec)

    # for down down up
    # Enecutoff2 = Msize0 - 1 + (Np-1)/2
    # matp2 = zeros(Int64,Msize0+1,Np-1)
    matp20 = zeros(Int64,Msize0+1,Np-1+1) # Np-1=2
    matp21 = zeros(Int64,Msize0+1,1+1)
    pascaltriangle!(Msize0,Np-1,matp20)
    pascaltriangle!(Msize0,1,matp21)
    indvec2 = cutMsizeEnespinmixed(Msize0,Np,matp20,matp21,Enecutoff,1)
    # indvec2 = cutMsizeEnespinmixed2(Msize0,Np,matp20,matp21,Enecutoff2,1)
    # indvec3 = cutMsizeEnespinmixed(Msize0,Np,matp20,matp21,Enecutoff,1)
    maxmatpcut2 = length(indvec2)
    # maxmatpcut3 = length(indvec3)
    
    # for up up down
    # matp30 = zeros(Int64,Msize0+1,Np-2+1) # Np-2=1
    # matp31 = zeros(Int64,Msize0+1,2+1)
    # pascaltriangle!(Msize0,Np-2,matp30)
    # pascaltriangle!(Msize0,2,matp31)
    # indvec3 = cutMsizeEnespinmixed(Msize0,Np,matp30,matp31,Enecutoff,2)
    # maxmatpcut3 = length(indvec3) # maxmatpcut3 == maxmatpcut2

    # construct transformation from 2nd quantisation basis to spin basis (only for downdownup and upupdown)
    mat_from2ndto1st_downup, mat_from1sttospin_downdownup, mat_from1sttospin_upupdown, mat_S2, mat_M2, mat_M4, mat_S3, mat_M1, mat_M3 = changefrom2ndtospin_downup(indvec2,Msize0,Np,matp20,matp21)
    # index_spin, basis_1st = changefrom2ndtospin_downup(indvec2,Msize0,Np,matp20,matp21)
    # return index_spin, basis_1st

    # save
    # save("data_mat_from2ndtospin_downup_$Msize0.Np3.jld", "Msize0", Msize0, "Np", Np, "mat_from2ndto1st_downup", mat_from2ndto1st_downup, "mat_from1sttospin_downdownup", mat_from1sttospin_downdownup, "mat_from1sttospin_upupdown", mat_from1sttospin_upupdown, "mat_S2", mat_S2, "mat_M2", mat_M2, "mat_M4", mat_M4, "mat_S3", mat_S3, "mat_M1", mat_M1, "mat_M3", mat_M3, "maxmatpcut", maxmatpcut, "maxmatpcut2", maxmatpcut2)
    
    return mat_from2ndto1st_downup, mat_from1sttospin_downdownup, mat_from1sttospin_upupdown, mat_S2, mat_M2, mat_M4, mat_S3, mat_M1, mat_M3, maxmatpcut, maxmatpcut2

end

function population_SM(phi::Union{Vector{Float64},Vector{ComplexF64}}, mat_from2ndto1st_downup::Union{Matrix{Float64},SparseMatrixCSC{Float64}}, mat_from1sttospin_downdownup::Union{Matrix{Float64},SparseMatrixCSC{Float64}}, mat_from1sttospin_upupdown::Union{Matrix{Float64},SparseMatrixCSC{Float64}}, mat_S2::SparseMatrixCSC{Int64}, mat_M2::SparseMatrixCSC{Int64}, mat_M4::SparseMatrixCSC{Int64}, mat_S3::SparseMatrixCSC{Int64}, mat_M1::SparseMatrixCSC{Int64}, mat_M3::SparseMatrixCSC{Int64}, maxmatpcut::Int64, maxmatpcut2::Int64)

    # construct transformation from 2nd quantisation basis to spin basis (only for downdownup and upupdown)
    # mat_from2ndto1st_downup, mat_from1sttospin_downdownup, mat_from1sttospin_upupdown, mat_S2, mat_M2, mat_M4, mat_S3, mat_M1, mat_M3, maxmatpcut, maxmatpcut2 = mat_from2ndtospin_downup(Msize0,Np)
    
    # population of S1 and S4
    pop_S1 = sum(abs.(phi[1:maxmatpcut]).^2)
    pop_S4 = sum(abs.(phi[end-maxmatpcut+1:end]).^2)

    # population of S2, S4, M1, M2, M3, M4
    phi_S2M2M4 = mat_from1sttospin_downdownup*(mat_from2ndto1st_downup*phi[maxmatpcut+1:maxmatpcut+maxmatpcut2])
    pop_S2 = sum(abs.(mat_S2*phi_S2M2M4).^2)
    pop_M2 = sum(abs.(mat_M2*phi_S2M2M4).^2)
    pop_M4 = sum(abs.(mat_M4*phi_S2M2M4).^2)
    phi_S3M1M3 = mat_from1sttospin_upupdown*(mat_from2ndto1st_downup*phi[maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2*2])
    pop_S3 = sum(abs.(mat_S3*phi_S3M1M3).^2)
    pop_M1 = sum(abs.(mat_M1*phi_S3M1M3).^2)
    pop_M3 = sum(abs.(mat_M3*phi_S3M1M3).^2)
    
    return pop_S1, pop_S2, pop_S3, pop_S4, pop_M1, pop_M2, pop_M3, pop_M4

end

function onebodydensitymatrix!(Msize0::Int64, Np::Int64, psi::Vector{ComplexF64}, rhoij::Matrix{ComplexF64}, rhoijdown::Matrix{ComplexF64}, rhoijup::Matrix{ComplexF64})

    # create Fock basis
    # for down down and up up
    Enecutoff = Msize0 - 1 + Np/2
    matp = zeros(Int64,Msize0+1,Np+1)
    pascaltriangle!(Msize0,Np,matp) # note the indices are m+1 and n+1 for N^m_n
    indvec = cutMsizeEnespinless(Msize0,Np,matp,Enecutoff)
    maxmatpcut = length(indvec)

    # for down down up
    matp20 = zeros(Int64,Msize0+1,Np-1+1) # Np-1=2
    matp21 = zeros(Int64,Msize0+1,1+1)
    pascaltriangle!(Msize0,Np-1,matp20)
    pascaltriangle!(Msize0,1,matp21)
    indvec2 = cutMsizeEnespinmixed(Msize0,Np,matp20,matp21,Enecutoff,1)
    maxmatpcut2 = length(indvec2)

    Msize = Msize0*2
    psijj = zeros(Float64,maxmatpcut+maxmatpcut2+maxmatpcut2+maxmatpcut)
    arraypsijj = zeros(Float64,maxmatpcut+maxmatpcut2+maxmatpcut2+maxmatpcut,Msize)
    # psiii = zeros(Float64,maxmatpcut+maxmatpcut2+maxmatpcut2+maxmatpcut)
    # rhoij = zeros(ComplexF64,maxmatpcut+maxmatpcut2+maxmatpcut2+maxmatpcut,maxmatpcut+maxmatpcut2+maxmatpcut2+maxmatpcut)
    # psi1 = copy(psi)

    rhoij .= 0.0

    for jj = 1:Msize
        coefficientonebody!(jj,psijj,Msize0,Np,maxmatpcut,maxmatpcut2,indvec,indvec2,matp,matp20,matp21)
        arraypsijj[:,jj] = psijj
        for ii = 1:jj
            rhoij[ii,jj] = (arraypsijj[:,ii].*psi)'*(arraypsijj[:,jj].*psi)
        end
    end
    rhoij .= rhoij/Np
    rhoij .= rhoij + rhoij' - diagm(diag(rhoij))

    rhoijdown .= rhoij[1:Msize0,1:Msize0]
    rhoijup .= rhoij[Msize0+1:Msize,Msize0+1:Msize]

    # rhoij .= Hermitian(rhoij)

    # return rhoij

end

function diagonaliseHtotsingle(Msize0::Int64, Np::Int64, gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega::Float64, specnum::Int64)

    matho, matdowndown, matupup, matdownup, matsoc, matW = createHtotal(Msize0,Np)
    # size_matho = size(matho)
    # lambda = zeros(ComplexF64,specnum)
    # phi = zeros(ComplexF64,specnum,size_matho[1])

    lambda, phi = eigs(matho + gdown*matdowndown + gup*matupup + gdu*matdownup + 1im*ksoc*matsoc + Omega*matW,nev=specnum,which=:SR)
    # lambda, phi = eigs(matho + gdown*matdowndown + gup*matupup + gdu*matdownup,nev=specnum,which=:SR)

    # decomp, history = partialschur(matho + gdown*matdowndown + gup*matupup + gdu*matdownup,nev=specnum,which=SR())

    # spect = real(lambda .- lambda[1])

    return lambda
    # return decomp

    # ~140 sec for Msize0=50 and specnum=6

end

function diagonaliseHtotspinpop_eigs(Msize0::Int64, Np::Int64, gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega::Float64, specnum::Int64)

    matho, matdowndown, matupup, matdownup, matsoc, matW = createHtotal(Msize0,Np)
    lambda, phi = eigs(matho + gdown*matdowndown + gup*matupup + gdu*matdownup + 1im*ksoc*matsoc + Omega*matW,nev=specnum,which=:SR)
    # lambda, phi = eigs(matho + gdown*matdowndown + gup*matupup + gdu*matdownup,nev=specnum,which=:SR)

    # for down3
    Enecutoff = Msize0 - 1 + Np/2
    matp = zeros(Int64,Msize0+1,Np+1)
    pascaltriangle!(Msize0,Np,matp) # note the indices are m+1 and n+1 for N^m_n
    indvec = cutMsizeEnespinless(Msize0,Np,matp,Enecutoff)
    maxmatpcut = length(indvec)

    # for down2up1
    matp20 = zeros(Int64,Msize0+1,Np-1+1) # Np-1=2
    matp21 = zeros(Int64,Msize0+1,1+1)
    pascaltriangle!(Msize0,Np-1,matp20)
    pascaltriangle!(Msize0,1,matp21)
    indvec2 = cutMsizeEnespinmixed(Msize0,Np,matp20,matp21,Enecutoff,1)
    maxmatpcut2 = length(indvec2)

    popdown3 = sum(real(abs.(phi[1:maxmatpcut,:]).^2),dims=1)'
    popdown2up1 = sum(real(abs.(phi[maxmatpcut+1:maxmatpcut+maxmatpcut2,:]).^2),dims=1)'
    popdown1up2 = sum(real(abs.(phi[maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2,:]).^2),dims=1)'
    popup3 = sum(real(abs.(phi[maxmatpcut+maxmatpcut2+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2+maxmatpcut,:]).^2),dims=1)'
    norm = popdown3 + popdown2up1 + popdown1up2 + popup3

    results = [lambda popdown3 popdown2up1 popdown1up2 popup3]
    # return results

    # rhoij = zeros(ComplexF64,Msize0*2,Msize0*2)
    # onebodydensitymatrix!(Msize0,Np,phi[:,1],rhoij)
    # lambdaconden, phiconden = eigs(rhoij,nev=5,which=:LR)

    rhoij = zeros(ComplexF64,Msize0*2,Msize0*2)
    rhoijdown = zeros(ComplexF64,Msize0,Msize0)
    rhoijup = zeros(ComplexF64,Msize0,Msize0)
    onebodydensitymatrix!(Msize0,Np,phi[:,1],rhoij,rhoijdown,rhoijup)
    lambdacondendown, phicondendown = eigs(rhoijdown,nev=5,which=:LR)
    lambdacondenup, phicondenup = eigs(rhoijup,nev=5,which=:LR)

    return lambdacondendown, phicondendown, lambdacondenup, phicondenup
    # return lambdaconden, phiconden
    # return rhoij

end

function diagonaliseHtotspinpop_test(Msize0::Int64, Np::Int64, gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega::Float64, specnum::Int64)

    # for down3
    Enecutoff = Msize0 - 1 + Np/2
    matp = zeros(Int64,Msize0+1,Np+1)
    pascaltriangle!(Msize0,Np,matp) # note the indices are m+1 and n+1 for N^m_n
    indvec = cutMsizeEnespinless(Msize0,Np,matp,Enecutoff)
    maxmatpcut = length(indvec)

    # for down2up1
    matp20 = zeros(Int64,Msize0+1,Np-1+1) # Np-1=2
    matp21 = zeros(Int64,Msize0+1,1+1)
    pascaltriangle!(Msize0,Np-1,matp20)
    pascaltriangle!(Msize0,1,matp21)
    indvec2 = cutMsizeEnespinmixed(Msize0,Np,matp20,matp21,Enecutoff,1)
    maxmatpcut2 = length(indvec2)

    matho, matdowndown, matupup, matdownup, matsoc, matW = createHtotal(Msize0,Np)
    mattot = Array(matho + gdown*matdowndown + gup*matupup + gdu*matdownup + 1im*ksoc*matsoc + Omega*matW)
    # mattot = Array(matho + gdown*matdowndown + gup*matupup + gdu*matdownup)
    mattot1 = Hermitian(mattot[1:maxmatpcut,1:maxmatpcut])
    mattot2 = Hermitian(mattot[maxmatpcut+1:maxmatpcut+maxmatpcut2,maxmatpcut+1:maxmatpcut+maxmatpcut2])
    mattot3 = Hermitian(mattot[maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2,maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2])
    mattot4 = Hermitian(mattot[maxmatpcut+maxmatpcut2+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2+maxmatpcut,maxmatpcut+maxmatpcut2+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2+maxmatpcut])
    lambda1, phi1 = eigen(mattot1,1:specnum)
    lambda2, phi2 = eigen(mattot2,1:specnum)
    lambda3, phi3 = eigen(mattot3,1:specnum)
    lambda4, phi4 = eigen(mattot4,1:specnum)

    results = [lambda1 lambda2 lambda3 lambda4]
    return results

end

function diagonaliseHtotspinpop(Msize0::Int64, Np::Int64, gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega::Float64, specnum::Int64)

    matho, matdowndown, matupup, matdownup, matsoc, matW = createHtotal(Msize0,Np)
    mattot = Hermitian(Array(matho + gdown*matdowndown + gup*matupup + gdu*matdownup + 1im*ksoc*matsoc + Omega*matW))
    # mattot = Hermitian(Array(matho + gdown*matdowndown + gup*matupup + gdu*matdownup))
    lambda, phi = eigen(mattot,1:specnum)
    # lambda, phi = eigs(mattot,nev=specnum,which=:SR)
    # lambda, phi = eigs(matho + gdown*matdowndown + gup*matupup + gdu*matdownup + 1im*ksoc*matsoc + Omega*matW,nev=specnum,which=:SR)
    # lambda, phi = eigs(matho + gdown*matdowndown + gup*matupup + gdu*matdownup,nev=specnum,which=:SR)

    # for down3
    Enecutoff = Msize0 - 1 + Np/2
    matp = zeros(Int64,Msize0+1,Np+1)
    pascaltriangle!(Msize0,Np,matp) # note the indices are m+1 and n+1 for N^m_n
    indvec = cutMsizeEnespinless(Msize0,Np,matp,Enecutoff)
    maxmatpcut = length(indvec)

    # for down2up1
    matp20 = zeros(Int64,Msize0+1,Np-1+1) # Np-1=2
    matp21 = zeros(Int64,Msize0+1,1+1)
    pascaltriangle!(Msize0,Np-1,matp20)
    pascaltriangle!(Msize0,1,matp21)
    indvec2 = cutMsizeEnespinmixed(Msize0,Np,matp20,matp21,Enecutoff,1)
    maxmatpcut2 = length(indvec2)

    popdown3 = sum(real(abs.(phi[1:maxmatpcut,:]).^2),dims=1)'
    popdown2up1 = sum(real(abs.(phi[maxmatpcut+1:maxmatpcut+maxmatpcut2,:]).^2),dims=1)'
    popdown1up2 = sum(real(abs.(phi[maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2,:]).^2),dims=1)'
    popup3 = sum(real(abs.(phi[maxmatpcut+maxmatpcut2+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2+maxmatpcut,:]).^2),dims=1)'
    norm = popdown3 + popdown2up1 + popdown1up2 + popup3

    results = [lambda popdown3 popdown2up1 popdown1up2 popup3]
    # return results

    rhoij = zeros(ComplexF64,Msize0*2,Msize0*2)
    rhoijdown = zeros(ComplexF64,Msize0,Msize0)
    rhoijup = zeros(ComplexF64,Msize0,Msize0)
    onebodydensitymatrix!(Msize0,Np,phi[:,1],rhoij,rhoijdown,rhoijup)
    lambdacondendown, phicondendown = eigs(rhoijdown,nev=5,which=:LR)
    lambdacondenup, phicondenup = eigs(rhoijup,nev=5,which=:LR)

    # lambdaconden, phiconden = eigs(rhoij,nev=5,which=:LR)
    # rhoij1 = Hermitian(rhoij)
    # lambdaconden, phiconden = eigen(rhoij1)

    # lambdaconden .= abs.(lambdaconden)

    return lambdacondendown, phicondendown, lambdacondenup, phicondenup

    # return lambdaconden
    # return lambdaconden, phiconden
    # return rhoij

end

function saveHtot(Msize0::Int64, Np::Int64)

    matho, matdowndown, matupup, matdownup, matsoc, matW = createHtotal(Msize0,Np)

    # save
    save("data_Htot$Msize0.Np3.jld", "Msize0", Msize0, "matho", matho, "matdowndown", matdowndown, "matupup", matupup, "matdownup", matdownup, "matsoc", matsoc, "matW", matW)

end

function diagonalisesavedHtot(Msize0::Int64, Np::Int64, gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega::Float64, specnum::Int64)

    # matho=load("data_Htot90_Np3.jld")["matho"]
    # matdowndown=load("data_Htot90_Np3.jld")["matdowndown"]
    # matdownup=load("data_Htot90_Np3.jld")["matdownup"]
    # matupup=load("data_Htot90_Np3.jld")["matupup"]
    # matsoc=load("data_Htot90_Np3.jld")["matsoc"]
    # matW=load("data_Htot90_Np3.jld")["matW"]

    matho=load("data_Htot$Msize0.Np3.jld")["matho"]
    matdowndown=load("data_Htot$Msize0.Np3.jld")["matdowndown"]
    matdownup=load("data_Htot$Msize0.Np3.jld")["matdownup"]
    matupup=load("data_Htot$Msize0.Np3.jld")["matupup"]
    matsoc=load("data_Htot$Msize0.Np3.jld")["matsoc"]
    matW=load("data_Htot$Msize0.Np3.jld")["matW"]

    # lambda, phi = eigs(matho + gdown*matdowndown + gup*matupup + gdu*matdownup,nev=specnum,which=:SR)
    lambda, phi = eigs(matho + gdown*matdowndown + gup*matupup + gdu*matdownup + 1im*ksoc*matsoc + Omega*matW,nev=specnum,which=:SR)

    # spect = real(lambda .- lambda[1])

    # for down3
    Enecutoff = Msize0 - 1 + Np/2
    matp = zeros(Int64,Msize0+1,Np+1)
    pascaltriangle!(Msize0,Np,matp) # note the indices are m+1 and n+1 for N^m_n
    indvec = cutMsizeEnespinless(Msize0,Np,matp,Enecutoff)
    maxmatpcut = length(indvec)

    # for down2up1
    matp20 = zeros(Int64,Msize0+1,Np-1+1) # Np-1=2
    matp21 = zeros(Int64,Msize0+1,1+1)
    pascaltriangle!(Msize0,Np-1,matp20)
    pascaltriangle!(Msize0,1,matp21)
    indvec2 = cutMsizeEnespinmixed(Msize0,Np,matp20,matp21,Enecutoff,1)
    maxmatpcut2 = length(indvec2)

    popdown3 = sum(real(abs.(phi[1:maxmatpcut,:]).^2),dims=1)'
    popdown2up1 = sum(real(abs.(phi[maxmatpcut+1:maxmatpcut+maxmatpcut2,:]).^2),dims=1)'
    popdown1up2 = sum(real(abs.(phi[maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2,:]).^2),dims=1)'
    popup3 = sum(real(abs.(phi[maxmatpcut+maxmatpcut2+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2+maxmatpcut,:]).^2),dims=1)'
    norm = popdown3 + popdown2up1 + popdown1up2 + popup3

    results = [lambda popdown3 popdown2up1 popdown1up2 popup3 norm]

    # return lambda, spect
    return results

end

# a good parameter set is Msize=50 and ksoc<=4

function diagonalisesavedHtotdiffW(Msize0::Int64, Np::Int64, gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega0::Float64, Omega1::Float64, NOmega::Int64, specnum::Int64)

    # matho=load("data_Htot90_Np3.jld")["matho"]
    # matdowndown=load("data_Htot90_Np3.jld")["matdowndown"]
    # matdownup=load("data_Htot90_Np3.jld")["matdownup"]
    # matupup=load("data_Htot90_Np3.jld")["matupup"]
    # matsoc=load("data_Htot90_Np3.jld")["matsoc"]
    # matW=load("data_Htot90_Np3.jld")["matW"]

    matho=load("data_Htot$Msize0.Np3.jld")["matho"]
    matdowndown=load("data_Htot$Msize0.Np3.jld")["matdowndown"]
    matdownup=load("data_Htot$Msize0.Np3.jld")["matdownup"]
    matupup=load("data_Htot$Msize0.Np3.jld")["matupup"]
    matsoc=load("data_Htot$Msize0.Np3.jld")["matsoc"]
    matW=load("data_Htot$Msize0.Np3.jld")["matW"]

    # for down3
    Enecutoff = Msize0 - 1 + Np/2
    matp = zeros(Int64,Msize0+1,Np+1)
    pascaltriangle!(Msize0,Np,matp) # note the indices are m+1 and n+1 for N^m_n
    indvec = cutMsizeEnespinless(Msize0,Np,matp,Enecutoff)
    maxmatpcut = length(indvec)

    # for down2up1
    matp20 = zeros(Int64,Msize0+1,Np-1+1) # Np-1=2
    matp21 = zeros(Int64,Msize0+1,1+1)
    pascaltriangle!(Msize0,Np-1,matp20)
    pascaltriangle!(Msize0,1,matp21)
    indvec2 = cutMsizeEnespinmixed(Msize0,Np,matp20,matp21,Enecutoff,1)
    maxmatpcut2 = length(indvec2)

    arrayOmega = LinRange(Omega0,Omega1,NOmega)
    arraylambda = zeros(ComplexF64,NOmega,specnum)
    arrayspect = zeros(ComplexF64,NOmega,specnum-1)
    arraypopdown3 = zeros(Float64,NOmega,specnum)
    arraypopdown2up1 = zeros(Float64,NOmega,specnum)
    arraypopdown1up2 = zeros(Float64,NOmega,specnum)
    arraypopup3 = zeros(Float64,NOmega,specnum)
    mat0 = matho + gdown*matdowndown + gup*matupup + gdu*matdownup + 1im*ksoc*matsoc
    mat1 = copy(mat0)
    phi = zeros(ComplexF64,maxmatpcut*2+maxmatpcut2*2,specnum)

    println("diagonalising the Hamiltonian for different Omega ...")
    for jj = 1:NOmega
        @time begin

            mat1 .= mat0 + arrayOmega[jj]*matW
            arraylambda[jj,:], phi = eigs(mat1,nev=specnum,which=:SR)
            arrayspect[jj,:] .= arraylambda[jj,2:end] .- arraylambda[jj,1]

            arraypopdown3[jj,:] = sum(abs.(phi[1:maxmatpcut,:]).^2,dims=1)'
            arraypopdown2up1[jj,:] = sum(abs.(phi[maxmatpcut+1:maxmatpcut+maxmatpcut2,:]).^2,dims=1)'
            arraypopdown1up2[jj,:] = sum(abs.(phi[maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2,:]).^2,dims=1)'
            arraypopup3[jj,:] = sum(abs.(phi[maxmatpcut+maxmatpcut2+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2+maxmatpcut,:]).^2,dims=1)'

            println(jj)

        end
    end

    save("data_spectrum.jld", "arrayOmega", arrayOmega, "arraylambda", arraylambda, "arrayspect", arrayspect, "arraypopdown3", arraypopdown3, "arraypopdown2up1", arraypopdown2up1, "arraypopdown1up2", arraypopdown1up2, "arraypopup3", arraypopup3)

    return arrayOmega, arraylambda, arrayspect, arraypopdown3, arraypopdown2up1, arraypopdown1up2, arraypopup3

end

function diagonalisesavedHtotdiffW_gdownup(Msize0::Int64, Np::Int64, gdown0::Float64, gdown1::Float64, Ng::Int64, gdu::Float64, ksoc::Float64, Omega0::Float64, Omega1::Float64, NOmega::Int64, specnum::Int64)

    # matho=load("data_Htot90_Np3.jld")["matho"]
    # matdowndown=load("data_Htot90_Np3.jld")["matdowndown"]
    # matdownup=load("data_Htot90_Np3.jld")["matdownup"]
    # matupup=load("data_Htot90_Np3.jld")["matupup"]
    # matsoc=load("data_Htot90_Np3.jld")["matsoc"]
    # matW=load("data_Htot90_Np3.jld")["matW"]

    matho=load("data_Htot$Msize0.Np3.jld")["matho"]
    matdowndown=load("data_Htot$Msize0.Np3.jld")["matdowndown"]
    matdownup=load("data_Htot$Msize0.Np3.jld")["matdownup"]
    matupup=load("data_Htot$Msize0.Np3.jld")["matupup"]
    matsoc=load("data_Htot$Msize0.Np3.jld")["matsoc"]
    matW=load("data_Htot$Msize0.Np3.jld")["matW"]

    # for down3
    Enecutoff = Msize0 - 1 + Np/2
    matp = zeros(Int64,Msize0+1,Np+1)
    pascaltriangle!(Msize0,Np,matp) # note the indices are m+1 and n+1 for N^m_n
    indvec = cutMsizeEnespinless(Msize0,Np,matp,Enecutoff)
    maxmatpcut = length(indvec)

    # for down2up1
    matp20 = zeros(Int64,Msize0+1,Np-1+1) # Np-1=2
    matp21 = zeros(Int64,Msize0+1,1+1)
    pascaltriangle!(Msize0,Np-1,matp20)
    pascaltriangle!(Msize0,1,matp21)
    indvec2 = cutMsizeEnespinmixed(Msize0,Np,matp20,matp21,Enecutoff,1)
    maxmatpcut2 = length(indvec2)

    arrayOmega = LinRange(Omega0,Omega1,NOmega)
    arraygdown = LinRange(gdown0,gdown1,Ng)

    arraylambda = zeros(ComplexF64,specnum,NOmega,Ng)
    arrayspect = zeros(ComplexF64,specnum-1,NOmega,Ng)
    arraypopdown3 = zeros(Float64,specnum,NOmega,Ng)
    arraypopdown2up1 = zeros(Float64,specnum,NOmega,Ng)
    arraypopdown1up2 = zeros(Float64,specnum,NOmega,Ng)
    arraypopup3 = zeros(Float64,specnum,NOmega,Ng)

    arraypopS1 = zeros(Float64,NOmega,Ng)
    arraypopS2 = zeros(Float64,NOmega,Ng)
    arraypopS3 = zeros(Float64,NOmega,Ng)
    arraypopS4 = zeros(Float64,NOmega,Ng)
    arraypopM1 = zeros(Float64,NOmega,Ng)
    arraypopM2 = zeros(Float64,NOmega,Ng)
    arraypopM3 = zeros(Float64,NOmega,Ng)
    arraypopM4 = zeros(Float64,NOmega,Ng)

    arrayenergyGStot = zeros(ComplexF64,NOmega,Ng)
    arrayenergyGSint = zeros(ComplexF64,NOmega,Ng)

    mat0 = matho + 1im*ksoc*matsoc
    mat1 = copy(mat0)
    matint = copy(mat0)
    phi = zeros(ComplexF64,maxmatpcut*2+maxmatpcut2*2,specnum)

    # construct transformation from 2nd quantisation basis to spin basis (only for downdownup and upupdown)
    mat_from2ndto1st_downup, mat_from1sttospin_downdownup, mat_from1sttospin_upupdown, mat_S2, mat_M2, mat_M4, mat_S3, mat_M1, mat_M3, maxmatpcut, maxmatpcut2 = mat_from2ndtospin_downup(Msize0,Np)

    println("diagonalising the Hamiltonian for different Omega ...")

    for jjg = 1:Ng

        gdownjjg = arraygdown[jjg]
        gupjjg = gdownjjg

        mat0 .= matho + gdownjjg*matdowndown + gupjjg*matupup + gdu*matdownup + 1im*ksoc*matsoc
        matint .= gdownjjg*matdowndown + gupjjg*matupup + gdu*matdownup

        println("jjg=",jjg)

        for jj = 1:NOmega
            @time begin

                mat1 .= mat0 + arrayOmega[jj]*matW
                arraylambda[:,jj,jjg], phi = eigs(mat1,nev=specnum,which=:SR)
                arrayspect[:,jj,jjg] .= arraylambda[2:end,jj,jjg] .- arraylambda[1,jj,jjg]

                arraypopdown3[:,jj,jjg] = sum(abs.(phi[1:maxmatpcut,:]).^2,dims=1)'
                arraypopdown2up1[:,jj,jjg] = sum(abs.(phi[maxmatpcut+1:maxmatpcut+maxmatpcut2,:]).^2,dims=1)'
                arraypopdown1up2[:,jj,jjg] = sum(abs.(phi[maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2,:]).^2,dims=1)'
                arraypopup3[:,jj,jjg] = sum(abs.(phi[maxmatpcut+maxmatpcut2+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2+maxmatpcut,:]).^2,dims=1)'

                # energy
                arrayenergyGStot[jj,jjg] = phi[:,1]'*mat1*phi[:,1] # same as arraylambda[1]
                arrayenergyGSint[jj,jjg] = phi[:,1]'*matint*phi[:,1]

                # population of S, M
                arraypopS1[jj,jjg], arraypopS2[jj,jjg], arraypopS3[jj,jjg], arraypopS4[jj,jjg], arraypopM1[jj,jjg], arraypopM2[jj,jjg], arraypopM3[jj,jjg], arraypopM4[jj,jjg] = population_SM(phi[:,1],mat_from2ndto1st_downup,mat_from1sttospin_downdownup,mat_from1sttospin_upupdown,mat_S2,mat_M2,mat_M4,mat_S3,mat_M1,mat_M3,maxmatpcut,maxmatpcut2)

                println("jj=",jj)

            end
        end

    end

    # indksoc = Int64(ksoc)
    # save("data_spectrum_ene_gdownup_jjg_ksoc$indksoc.jld", "arrayOmega", arrayOmega, "arraygdown", arraygdown, "ksoc", ksoc, "arraylambda", arraylambda, "arrayspect", arrayspect, "arraypopdown3", arraypopdown3, "arraypopdown2up1", arraypopdown2up1, "arraypopdown1up2", arraypopdown1up2, "arraypopup3", arraypopup3, "arrayenergyGStot", arrayenergyGStot, "arrayenergyGSint", arrayenergyGSint)

    return arrayOmega, arraygdown, ksoc, arraylambda, arrayspect, arraypopdown3, arraypopdown2up1, arraypopdown1up2, arraypopup3, arrayenergyGStot, arrayenergyGSint, arraypopS1, arraypopS2, arraypopS3, arraypopS4, arraypopM1, arraypopM2, arraypopM3, arraypopM4

end

function diagonalisesavedHtotdiffW_gdu(Msize0::Int64, Np::Int64, gdu0::Float64, gdu1::Float64, Ng::Int64, gdownup::Float64, ksoc::Float64, Omega0::Float64, Omega1::Float64, NOmega::Int64, specnum::Int64)

    matho=load("data_Htot$Msize0.Np3.jld")["matho"]
    matdowndown=load("data_Htot$Msize0.Np3.jld")["matdowndown"]
    matdownup=load("data_Htot$Msize0.Np3.jld")["matdownup"]
    matupup=load("data_Htot$Msize0.Np3.jld")["matupup"]
    matsoc=load("data_Htot$Msize0.Np3.jld")["matsoc"]
    matW=load("data_Htot$Msize0.Np3.jld")["matW"]

    # matho=load("data_Htot90_Np3.jld")["matho"]
    # matdowndown=load("data_Htot90_Np3.jld")["matdowndown"]
    # matdownup=load("data_Htot90_Np3.jld")["matdownup"]
    # matupup=load("data_Htot90_Np3.jld")["matupup"]
    # matsoc=load("data_Htot90_Np3.jld")["matsoc"]
    # matW=load("data_Htot90_Np3.jld")["matW"]
    # matho, matdowndown, matupup, matdownup, matsoc, matW = createHtotal(Msize0,Np)

    # for down3
    Enecutoff = Msize0 - 1 + Np/2
    matp = zeros(Int64,Msize0+1,Np+1)
    pascaltriangle!(Msize0,Np,matp) # note the indices are m+1 and n+1 for N^m_n
    indvec = cutMsizeEnespinless(Msize0,Np,matp,Enecutoff)
    maxmatpcut = length(indvec)

    # for down2up1
    matp20 = zeros(Int64,Msize0+1,Np-1+1) # Np-1=2
    matp21 = zeros(Int64,Msize0+1,1+1)
    pascaltriangle!(Msize0,Np-1,matp20)
    pascaltriangle!(Msize0,1,matp21)
    indvec2 = cutMsizeEnespinmixed(Msize0,Np,matp20,matp21,Enecutoff,1)
    maxmatpcut2 = length(indvec2)

    arrayOmega = LinRange(Omega0,Omega1,NOmega)
    arraygdu = LinRange(gdu0,gdu1,Ng)

    arraylambda = zeros(ComplexF64,specnum,NOmega,Ng)
    arrayspect = zeros(ComplexF64,specnum-1,NOmega,Ng)
    arraypopdown3 = zeros(Float64,specnum,NOmega,Ng)
    arraypopdown2up1 = zeros(Float64,specnum,NOmega,Ng)
    arraypopdown1up2 = zeros(Float64,specnum,NOmega,Ng)
    arraypopup3 = zeros(Float64,specnum,NOmega,Ng)

    arraypopS1 = zeros(Float64,NOmega,Ng)
    arraypopS2 = zeros(Float64,NOmega,Ng)
    arraypopS3 = zeros(Float64,NOmega,Ng)
    arraypopS4 = zeros(Float64,NOmega,Ng)
    arraypopM1 = zeros(Float64,NOmega,Ng)
    arraypopM2 = zeros(Float64,NOmega,Ng)
    arraypopM3 = zeros(Float64,NOmega,Ng)
    arraypopM4 = zeros(Float64,NOmega,Ng)
   
    arrayenergyGStot = zeros(ComplexF64,NOmega,Ng)
    arrayenergyGSint = zeros(ComplexF64,NOmega,Ng)

    # mat0 = matho + gdown*matdowndown + gup*matupup + gdu*matdownup + 1im*ksoc*matsoc
    mat0 = matho + 1im*ksoc*matsoc
    mat1 = copy(mat0)
    matint = copy(mat0)
    phi = zeros(ComplexF64,maxmatpcut*2+maxmatpcut2*2,specnum)

    # construct transformation from 2nd quantisation basis to spin basis (only for downdownup and upupdown)
    mat_from2ndto1st_downup, mat_from1sttospin_downdownup, mat_from1sttospin_upupdown, mat_S2, mat_M2, mat_M4, mat_S3, mat_M1, mat_M3, maxmatpcut, maxmatpcut2 = mat_from2ndtospin_downup(Msize0,Np)

    println("diagonalising the Hamiltonian for different Omega ...")

    for jjg = 1:Ng

        gdujjg = arraygdu[jjg]

        mat0 .= matho + gdownup*matdowndown + gdownup*matupup + gdujjg*matdownup + 1im*ksoc*matsoc
        matint .= gdownup*matdowndown + gdownup*matupup + gdujjg*matdownup

        println("jjg=",jjg)

        for jj = 1:NOmega
            @time begin

                mat1 .= mat0 + arrayOmega[jj]*matW
                arraylambda[:,jj,jjg], phi = eigs(mat1,nev=specnum,which=:SR)
                arrayspect[:,jj,jjg] .= arraylambda[2:end,jj,jjg] .- arraylambda[1,jj,jjg]

                arraypopdown3[:,jj,jjg] = sum(abs.(phi[1:maxmatpcut,:]).^2,dims=1)'
                arraypopdown2up1[:,jj,jjg] = sum(abs.(phi[maxmatpcut+1:maxmatpcut+maxmatpcut2,:]).^2,dims=1)'
                arraypopdown1up2[:,jj,jjg] = sum(abs.(phi[maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2,:]).^2,dims=1)'
                arraypopup3[:,jj,jjg] = sum(abs.(phi[maxmatpcut+maxmatpcut2+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2+maxmatpcut,:]).^2,dims=1)'

                # energy
                arrayenergyGStot[jj,jjg] = phi[:,1]'*mat1*phi[:,1] # same as arraylambda[1]
                arrayenergyGSint[jj,jjg] = phi[:,1]'*matint*phi[:,1]

                # population of S, M
                arraypopS1[jj,jjg], arraypopS2[jj,jjg], arraypopS3[jj,jjg], arraypopS4[jj,jjg], arraypopM1[jj,jjg], arraypopM2[jj,jjg], arraypopM3[jj,jjg], arraypopM4[jj,jjg] = population_SM(phi[:,1],mat_from2ndto1st_downup,mat_from1sttospin_downdownup,mat_from1sttospin_upupdown,mat_S2,mat_M2,mat_M4,mat_S3,mat_M1,mat_M3,maxmatpcut,maxmatpcut2)

                println("jj=",jj)

            end
        end

    end

    # indksoc = Int64(ksoc)
    # save("data_spectrum_ene_gdu_jjg_ksoc$indksoc.jld", "arrayOmega", arrayOmega, "arraygdu", arraygdu, "ksoc", ksoc, "arraylambda", arraylambda, "arrayspect", arrayspect, "arraypopdown3", arraypopdown3, "arraypopdown2up1", arraypopdown2up1, "arraypopdown1up2", arraypopdown1up2, "arraypopup3", arraypopup3, "arrayenergyGStot", arrayenergyGStot, "arrayenergyGSint", arrayenergyGSint)

    return arrayOmega, arraygdu, ksoc, arraylambda, arrayspect, arraypopdown3, arraypopdown2up1, arraypopdown1up2, arraypopup3, arrayenergyGStot, arrayenergyGSint, arraypopS1, arraypopS2, arraypopS3, arraypopS4, arraypopM1, arraypopM2, arraypopM3, arraypopM4

end

function plot_spinpop(arrayOmega,arraypopdown3, arraypopdown2up1, arraypopdown1up2, arraypopup3,arraypopS1,arraypopS2,arraypopS3,arraypopS4,arraypopM1,arraypopM2,arraypopM3,arraypopM4,jj1,jj2)

    NOmega = length(arrayOmega)

    plot(arrayOmega,ones(NOmega)*3/4,color=:black,ls=:dot,lw=2)
    plot!(arrayOmega,ones(NOmega)*1/4,color=:black,ls=:dot,lw=2)
    plot!([arrayOmega[jj2],arrayOmega[jj2]+1e-5],[-10,10],color=:black,lw=2,ls=:dot)

    plot!(arrayOmega,arraypopdown3[1,:,jj1]+arraypopup3[1,:,jj1],color=:grey,lw=4)
    plot!(arrayOmega,arraypopdown2up1[1,:,jj1]+arraypopdown1up2[1,:,jj1],color=:black,lw=4)
    plot!(arrayOmega,arraypopM1[:,jj1]+arraypopM2[:,jj1]+arraypopM3[:,jj1]+arraypopM4[:,jj1],lw=3,color=:red,ls=:dot)
    plot!(arrayOmega,arraypopS2[:,jj1]+arraypopS3[:,jj1],lw=3,ls=:dash,color=:blue)
    # plot!(arrayOmega,arraypopS1[:,jj1]+arraypopS4[:,jj1],lw=2)

    plot!(legend=:none)
    
    xlims!((9,41))
    ylims!((-0.05,1.05))
    plot!(aspect_ratio=32/1.1)

end

function plot_eneint(arrayOmega,arrayenergyGSint_g1,arrayenergyGSint_g2,arrayenergyGSint_g10)

    plot(arrayOmega,real(arrayenergyGSint_g1[:,1]),lw=3)
    plot!(arrayOmega,real(arrayenergyGSint_g2[:,1]),lw=3,ls=:dash)
    plot!(arrayOmega,real(arrayenergyGSint_g10[:,1]),lw=3,ls=:dot)
    
    plot!(legend=:none)
    
    xlims!((9,41))
    ylims!((0.0,0.51))
    plot!(aspect_ratio=32/0.51)

end

function plot_arrayspect(arrayOmega,arrayspect,jj1,jj2)

    plot([2*4^2,2*4^2+1e-5],[-10,10],color=:black,lw=2,ls=:dash)
    plot!([arrayOmega[jj2],arrayOmega[jj2]+1e-5],[-10,10],color=:black,lw=2,ls=:dot)

    plot!(arrayOmega,real(arrayspect[1,:,jj1]),color=:blue,lw=4)
    plot!(arrayOmega,real(arrayspect[2,:,jj1]),color=:blue,lw=4)
    plot!(arrayOmega,real(arrayspect[3,:,jj1]),color=:blue,lw=4)
    plot!(arrayOmega,real(arrayspect[4,:,jj1]),color=:blue,lw=4)
    plot!(arrayOmega,real(arrayspect[5,:,jj1]),color=:blue,lw=4)

    plot!(legend=:none)
    
    xlims!((9,41))
    ylims!((-0.05,1.2))
    plot!(aspect_ratio=32/1.25)

end

function plottest0(arraypopdown2up1)

    test = arraypopdown2up1[1,:,:].-3/8
    size0 = size(arraypopdown2up1[1,:,:])

    for jj = 1:size0[1]
        for ii = 1:size0[2]
            if test[jj,ii] < 0.0
                test[jj,ii] = NaN
            end
        end
    end

    return test

end

function plottest1(arraypopdown3)

    test = arraypopdown3[1,:,:].-1/8
    size0 = size(arraypopdown3[1,:,:])

    for jj = 1:size0[1]
        for ii = 1:size0[2]
            if test[jj,ii] < 0.0
                test[jj,ii] = NaN
            end
        end
    end

    return test

end

function plotenergyerror()

    arrayk = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 5.5, 6.0]

    energyg0 = [1.5000000000000175, 1.5, 1.4999999999998943, 1.5000000000001439, 1.500031718214867, 1.5999068416323183, 2.0359027534663454, 3.024184652434016]

    energyg1 = [2.4197868676065117, 2.420144673581986, 2.421841420906276, 2.4254587195358415, 2.433722389339515, 2.518617064258919, 2.906251219386661, 3.861473180998914]

    energyg10 = [4.223469871621669, 4.226212679791851, 4.239303553268347, 4.267740113334156, 4.336313301135704, 4.751150726067593, 5.422254241434274, 6.622550679838369]

    # plot(arrayk,energyg0)

    return arrayk, energyg0, energyg1, energyg10

end

function plotpairTG(Lx,Nx)

    xrange = LinRange(-Lx,Lx,Nx)
    yrange = LinRange(-Lx,Lx,Nx)
    fun = zeros(Float64,Nx,Nx)
    
    for jjx = 1:Nx
        xx = xrange[jjx]
        for jjy = 1:Nx
            yy = yrange[jjy]
            fun[jjx,jjy] = exp(-xx^2-yy^2)*(xx-yy)^2*(2*xx^2+2*yy^2+8*xx*yy+4*xx^2*yy^2+3)/(6*pi)
            #fun[jjx,jjy] = exp(-2*(xx^2+yy^2))*(xx-yy)^2*(3+4xx^2+4*yy^2+16*xx^2*yy^2+16*xx*yy)
        end
    end
    # fun[:,:] .= fun[:,:]*sqrt(pi/2)/16/(3*sqrt(pi)^3/32/sqrt(2))

    return xrange, yrange, fun

end

function diagonaliseH_onebody_test(Msize0::Int64, Np::Int64, gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega::Float64, specnum::Int64, Lx::Float64, Nx::Int64)

    println("constructoing the Hamiltonian ...")
    @time begin
        matho, matdowndown, matupup, matdownup, matsoc, matW = createHtotal(Msize0,Np)
        # save("data_Htot.jld", "matho", matho, "matsoc", matsoc, "matW", matW, "mat1", mat1, "mat2", mat2, "mat3", mat3)
    end

    # for down3
    Enecutoff = Msize0 - 1 + Np/2
    matp = zeros(Int64,Msize0+1,Np+1)
    pascaltriangle!(Msize0,Np,matp) # note the indices are m+1 and n+1 for N^m_n
    indvec = cutMsizeEnespinless(Msize0,Np,matp,Enecutoff)
    maxmatpcut = length(indvec)

    # for down2up1
    matp20 = zeros(Int64,Msize0+1,Np-1+1) # Np-1=2
    matp21 = zeros(Int64,Msize0+1,1+1)
    pascaltriangle!(Msize0,Np-1,matp20)
    pascaltriangle!(Msize0,1,matp21)
    indvec2 = cutMsizeEnespinmixed(Msize0,Np,matp20,matp21,Enecutoff,1)
    maxmatpcut2 = length(indvec2)

    arraylambda = zeros(ComplexF64,specnum)
    arrayspect = zeros(ComplexF64,specnum-1)

    # arraypopdown3 = zeros(Float64,specnum,NOmega,Ng)
    # arraypopdown2up1 = zeros(Float64,specnum,NOmega,Ng)
    # arraypopdown1up2 = zeros(Float64,specnum,NOmega,Ng)
    # arraypopup3 = zeros(Float64,specnum,NOmega,Ng)

    # arraypopdown3GSspatial = zeros(Float64,maxmatpcut,NOmega,Ng)
    # arraypopdown2up1GSspatial = zeros(Float64,maxmatpcut2,NOmega,Ng)
    # arraypopdown1up2GSspatial = zeros(Float64,maxmatpcut2,NOmega,Ng)
    # arraypopup3GSspatial = zeros(Float64,maxmatpcut,NOmega,Ng)

    rhoij = zeros(ComplexF64,Msize0*2,Msize0*2)
    rhoijdown = zeros(ComplexF64,Msize0,Msize0)
    rhoijup = zeros(ComplexF64,Msize0,Msize0)
    lambdaconden = zeros(ComplexF64,specnum)
    phiconden = zeros(ComplexF64,Msize0*2,specnum)

    # mat0 = matho + gdown*matdowndown + gup*matupup + gdu*matdownup + 1im*ksoc*matsoc
    mat0 = matho + 1im*ksoc*matsoc
    # mat1 = spzeros(ComplexF64,maxmatpcut+maxmatpcut2*2+maxmatpcut,maxmatpcut+maxmatpcut2*2+maxmatpcut)
    mat1 = copy(mat0)
    psi = zeros(ComplexF64,maxmatpcut*2+maxmatpcut2*2,specnum)
    xrange = LinRange(-Lx,Lx,Nx)
    yrange = LinRange(-Lx,Lx,Nx)

    # Mpairdown3int = zeros(Int64,maxmatpcut+maxmatpcut2,maxmatpcut+maxmatpcut2,12)
    # Mpairdown3int1 = zeros(Int64,maxmatpcut,maxmatpcut,12)
    # Mpairdown3int2 = zeros(Int64,maxmatpcut2,maxmatpcut2,12)
    # Mpairdown2up1int = zeros(Int64,maxmatpcut2,maxmatpcut2,12)
    # Mpairdown1up2int = zeros(Int64,maxmatpcut2,maxmatpcut2,12)
    # Mpairup3int = zeros(Int64,maxmatpcut+maxmatpcut2,maxmatpcut+maxmatpcut2,12)
    # Mpairup3int1 = zeros(Int64,maxmatpcut2,maxmatpcut2,12)
    # Mpairup3int2 = zeros(Int64,maxmatpcut,maxmatpcut,12)
    # Mpairdown3int = spzeros(Int64,(maxmatpcut+maxmatpcut2)^2,12)
    # Mpairdown2up1int = spzeros(Int64,maxmatpcut2^2,12)

    # Mpairdown3float = zeros(Float64,maxmatpcut+maxmatpcut2,maxmatpcut+maxmatpcut2,3)
    # Mpairdown3float1 = zeros(Float64,maxmatpcut,maxmatpcut,3)
    # Mpairdown3float2 = zeros(Float64,maxmatpcut2,maxmatpcut2,3)
    # Mpairdown2up1float = zeros(Float64,maxmatpcut2,maxmatpcut2,3)
    # Mpairdown1up2float = zeros(Float64,maxmatpcut2,maxmatpcut2,3)
    # Mpairup3float = zeros(Float64,maxmatpcut+maxmatpcut2,maxmatpcut+maxmatpcut2,3)
    # Mpairup3float1 = zeros(Float64,maxmatpcut2,maxmatpcut2,3)
    # Mpairup3float2 = zeros(Float64,maxmatpcut,maxmatpcut,3)
    # Mpairdown3float = spzeros(Float64,(maxmatpcut+maxmatpcut2)^2,3)
    # Mpairdown2up1float = spzeros(Float64,maxmatpcut2^2,3)

    println("diagonalising the Hamiltonian ...")
    @time begin
        mat0 .= matho + gdown*matdowndown + gup*matupup + gdu*matdownup + 1im*ksoc*matsoc
        mat1 .= mat0 + Omega*matW
        # mattot = Hermitian(Array(matho + gdown*matdowndown + gup*matupup + gdu*matdownup + 1im*ksoc*matsoc + Omega*matW))
        arraylambda, psi = eigs(mat1,nev=specnum,which=:SR)
        # arraylambda, psi = eigen(mattot,1:specnum)
        arrayspect .= arraylambda[2:end] .- arraylambda[1]
    end

    # println("calculating one body density matrix ...")
    # @time begin
    #     onebodydensitymatrix!(Msize0,Np,phi[:,1],rhoij,rhoijdown,rhoijup)
    #     lambdacondendown, phicondendown = eigs(rhoijdown,nev=5,which=:LR)
    #     lambdacondenup, phicondenup = eigs(rhoijup,nev=5,which=:LR)
    #     # onebodydensitymatrix!(Msize0,Np,phi[:,1],rhoij)
    #     # lambdaconden, phiconden = eigs(rhoij,nev=specnum,which=:LR)
    # end

    # psi .= 0.0
    # psi[:,1] = ones(maxmatpcut*2+maxmatpcut2*2)
    # psi[1,1] = 1.0
    # psi[2,1] = 1.0
    # psi[3,1] = 1.0
    # psi[:,1] = psi[:,1]/sqrt(sum(abs.(psi[:,1]).^2))

    println("calculating pair correlation ...")
    @time begin
        fun_nudown, fun_nudu, fun_nuup = paircorrelation_fun(indvec,indvec2,Msize0,Np,matp,matp20,matp21,psi[:,1],xrange,yrange)
        # fun_nu2 = paircorrelation_fun(indvec,indvec2,Msize0,Np,matp,matp20,matp21,psi[:,2],xrange,yrange)
        # fun_nu3 = paircorrelation_fun(indvec,indvec2,Msize0,Np,matp,matp20,matp21,psi[:,3],xrange,yrange)
        # fun_nu4 = paircorrelation_fun(indvec,indvec2,Msize0,Np,matp,matp20,matp21,psi[:,4],xrange,yrange)
    end

    # return arraylambda, arrayspect, lambdacondendown, phicondendown, lambdacondenup, phicondenup
    return arraylambda, psi, xrange, yrange, fun_nudown, fun_nudu, fun_nuup

end

function diagonalisesavedHtotdiffW_gdownup_onebody(Msize0::Int64, Np::Int64, gdown0::Float64, gdown1::Float64, Ng::Int64, gdu::Float64, ksoc::Float64, Omega0::Float64, Omega1::Float64, NOmega::Int64, specnum::Int64)

    # println("constructoing the Hamiltonian ...")
    # @time begin
    #     matho, matdowndown, matupup, matdownup, matsoc, matW = createHtotal(Msize0,Np)
    #     # save("data_Htot.jld", "matho", matho, "matsoc", matsoc, "matW", matW, "mat1", mat1, "mat2", mat2, "mat3", mat3)
    # end

    matho=load("data_Htot90_Np3.jld")["matho"]
    matdowndown=load("data_Htot90_Np3.jld")["matdowndown"]
    matdownup=load("data_Htot90_Np3.jld")["matdownup"]
    matupup=load("data_Htot90_Np3.jld")["matupup"]
    matsoc=load("data_Htot90_Np3.jld")["matsoc"]
    matW=load("data_Htot90_Np3.jld")["matW"]
    # matho, matdowndown, matupup, matdownup, matsoc, matW = createHtotal(Msize0,Np)

    # for down3
    Enecutoff = Msize0 - 1 + Np/2
    matp = zeros(Int64,Msize0+1,Np+1)
    pascaltriangle!(Msize0,Np,matp) # note the indices are m+1 and n+1 for N^m_n
    indvec = cutMsizeEnespinless(Msize0,Np,matp,Enecutoff)
    maxmatpcut = length(indvec)

    # for down2up1
    matp20 = zeros(Int64,Msize0+1,Np-1+1) # Np-1=2
    matp21 = zeros(Int64,Msize0+1,1+1)
    pascaltriangle!(Msize0,Np-1,matp20)
    pascaltriangle!(Msize0,1,matp21)
    indvec2 = cutMsizeEnespinmixed(Msize0,Np,matp20,matp21,Enecutoff,1)
    maxmatpcut2 = length(indvec2)

    arrayOmega = LinRange(Omega0,Omega1,NOmega)
    arraygdown = LinRange(gdown0,gdown1,Ng)

    arraylambda = zeros(ComplexF64,specnum,NOmega,Ng)
    arrayspect = zeros(ComplexF64,specnum-1,NOmega,Ng)

    # arraypopdown3 = zeros(Float64,specnum,NOmega,Ng)
    # arraypopdown2up1 = zeros(Float64,specnum,NOmega,Ng)
    # arraypopdown1up2 = zeros(Float64,specnum,NOmega,Ng)
    # arraypopup3 = zeros(Float64,specnum,NOmega,Ng)

    # arraypopdown3GSspatial = zeros(Float64,maxmatpcut,NOmega,Ng)
    # arraypopdown2up1GSspatial = zeros(Float64,maxmatpcut2,NOmega,Ng)
    # arraypopdown1up2GSspatial = zeros(Float64,maxmatpcut2,NOmega,Ng)
    # arraypopup3GSspatial = zeros(Float64,maxmatpcut,NOmega,Ng)

    rhoij = zeros(ComplexF64,Msize0*2,Msize0*2)
    rhoijdown = zeros(ComplexF64,Msize0,Msize0)
    rhoijup = zeros(ComplexF64,Msize0,Msize0)

    arraylambdacondendown = zeros(ComplexF64,specnum,NOmega,Ng)
    arraylambdacondenup = zeros(ComplexF64,specnum,NOmega,Ng)
    arrayphicondendown = zeros(ComplexF64,Msize0,specnum,NOmega,Ng)
    arrayphicondenup = zeros(ComplexF64,Msize0,specnum,NOmega,Ng)

    # mat0 = matho + gdown*matdowndown + gup*matupup + gdu*matdownup + 1im*ksoc*matsoc
    mat0 = matho + 1im*ksoc*matsoc
    # mat1 = spzeros(ComplexF64,maxmatpcut+maxmatpcut2*2+maxmatpcut,maxmatpcut+maxmatpcut2*2+maxmatpcut)
    mat1 = copy(mat0)
    phi = zeros(ComplexF64,maxmatpcut*2+maxmatpcut2*2,specnum)

    # println("diagonalising the Hamiltonian for different Omega ...")

    for jjg = 1:Ng

        gdownjjg = arraygdown[jjg]
        gupjjg = gdownjjg

        mat0 .= matho + gdownjjg*matdowndown + gupjjg*matupup + gdu*matdownup + 1im*ksoc*matsoc

        println("jjg=",jjg)

        for jj = 1:NOmega # parfor
            @time begin

                mat1 .= mat0 + arrayOmega[jj]*matW
                arraylambda[:,jj,jjg], phi = eigs(mat1,nev=specnum,which=:SR)
                arrayspect[:,jj,jjg] .= arraylambda[2:end,jj,jjg] .- arraylambda[1,jj,jjg]

                # onebodydensitymatrix!(Msize0,Np,phi[:,1],rhoij)
                # lambdaconden, phiconden = eigs(rhoij,nev=specnum,which=:LR)

                rhoij .= 0.0
                rhoijdown .= 0.0
                rhoijup .= 0.0
                onebodydensitymatrix!(Msize0,Np,phi[:,1],rhoij,rhoijdown,rhoijup)
                lambdacondendown, phicondendown = eigs(rhoijdown,nev=specnum,which=:LR)
                lambdacondenup, phicondenup = eigs(rhoijup,nev=specnum,which=:LR)

                arraylambdacondendown[:,jj,jjg] = lambdacondendown
                arraylambdacondenup[:,jj,jjg] = lambdacondenup
                arrayphicondendown[:,:,jj,jjg] = phicondendown
                arrayphicondenup[:,:,jj,jjg] = phicondenup

                # arraypopdown3[:,jj,jjg] = sum(abs.(phi[1:maxmatpcut,:]).^2,dims=1)'
                # arraypopdown2up1[:,jj,jjg] = sum(abs.(phi[maxmatpcut+1:maxmatpcut+maxmatpcut2,:]).^2,dims=1)'
                # arraypopdown1up2[:,jj,jjg] = sum(abs.(phi[maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2,:]).^2,dims=1)'
                # arraypopup3[:,jj,jjg] = sum(abs.(phi[maxmatpcut+maxmatpcut2+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2+maxmatpcut,:]).^2,dims=1)'

                # arraypopdown3GSspatial[:,jj,jjg] = phi[1:maxmatpcut,1]
                # arraypopdown2up1GSspatial[:,jj,jjg] = phi[maxmatpcut+1:maxmatpcut+maxmatpcut2,1]
                # arraypopdown1up2GSspatial[:,jj,jjg] = phi[maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2,1]
                # arraypopup3GSspatial[:,jj,jjg] = phi[maxmatpcut+maxmatpcut2+maxmatpcut2+1:end,1]

                println("jj=",jj)

            end
        end

    end

    indksoc = Int64(ksoc)
    save("data_spectrum_onebody_gdownup_jjg_ksoc$indksoc.jld", "arrayOmega", arrayOmega, "arraygdown", arraygdown, "ksoc", ksoc, "arraylambda", arraylambda, "arrayspect", arrayspect, "arraylambdacondendown", arraylambdacondendown, "arraylambdacondenup", arraylambdacondenup, "arrayphicondendown", arrayphicondendown, "arrayphicondenup", arrayphicondenup)

    # return arrayOmega, arraygdown, ksoc, arraylambda, arrayspect, arraypopdown3, arraypopdown2up1, arraypopdown1up2, arraypopup3

end

function diagonalisesavedHtotdiffW_gdu_onebody(Msize0::Int64, Np::Int64, gdu0::Float64, gdu1::Float64, Ng::Int64, gdownup::Float64, ksoc::Float64, Omega0::Float64, Omega1::Float64, NOmega::Int64, specnum::Int64)

    # println("constructoing the Hamiltonian ...")
    # @time begin
    #     matho, matdowndown, matupup, matdownup, matsoc, matW = createHtotal(Msize0,Np)
    #     # save("data_Htot.jld", "matho", matho, "matsoc", matsoc, "matW", matW, "mat1", mat1, "mat2", mat2, "mat3", mat3)
    # end

    matho=load("data_Htot90_Np3.jld")["matho"]
    matdowndown=load("data_Htot90_Np3.jld")["matdowndown"]
    matdownup=load("data_Htot90_Np3.jld")["matdownup"]
    matupup=load("data_Htot90_Np3.jld")["matupup"]
    matsoc=load("data_Htot90_Np3.jld")["matsoc"]
    matW=load("data_Htot90_Np3.jld")["matW"]
    # matho, matdowndown, matupup, matdownup, matsoc, matW = createHtotal(Msize0,Np)

    # for down3
    Enecutoff = Msize0 - 1 + Np/2
    matp = zeros(Int64,Msize0+1,Np+1)
    pascaltriangle!(Msize0,Np,matp) # note the indices are m+1 and n+1 for N^m_n
    indvec = cutMsizeEnespinless(Msize0,Np,matp,Enecutoff)
    maxmatpcut = length(indvec)

    # for down2up1
    matp20 = zeros(Int64,Msize0+1,Np-1+1) # Np-1=2
    matp21 = zeros(Int64,Msize0+1,1+1)
    pascaltriangle!(Msize0,Np-1,matp20)
    pascaltriangle!(Msize0,1,matp21)
    indvec2 = cutMsizeEnespinmixed(Msize0,Np,matp20,matp21,Enecutoff,1)
    maxmatpcut2 = length(indvec2)

    arrayOmega = LinRange(Omega0,Omega1,NOmega)
    arraygdu = LinRange(gdu0,gdu1,Ng)

    arraylambda = zeros(ComplexF64,specnum,NOmega,Ng)
    arrayspect = zeros(ComplexF64,specnum-1,NOmega,Ng)

    # arraypopdown3 = zeros(Float64,specnum,NOmega,Ng)
    # arraypopdown2up1 = zeros(Float64,specnum,NOmega,Ng)
    # arraypopdown1up2 = zeros(Float64,specnum,NOmega,Ng)
    # arraypopup3 = zeros(Float64,specnum,NOmega,Ng)

    # arraypopdown3GSspatial = zeros(Float64,maxmatpcut,NOmega,Ng)
    # arraypopdown2up1GSspatial = zeros(Float64,maxmatpcut2,NOmega,Ng)
    # arraypopdown1up2GSspatial = zeros(Float64,maxmatpcut2,NOmega,Ng)
    # arraypopup3GSspatial = zeros(Float64,maxmatpcut,NOmega,Ng)

    # arrayenergyGStot = zeros(Float64,NOmega,Ng)
    # arrayenergyGSint = zeros(Float64,NOmega,Ng)

    rhoij = zeros(ComplexF64,Msize0*2,Msize0*2)
    rhoijdown = zeros(ComplexF64,Msize0,Msize0)
    rhoijup = zeros(ComplexF64,Msize0,Msize0)

    arraylambdacondendown = zeros(ComplexF64,specnum,NOmega,Ng)
    arraylambdacondenup = zeros(ComplexF64,specnum,NOmega,Ng)
    arrayphicondendown = zeros(ComplexF64,Msize0,specnum,NOmega,Ng)
    arrayphicondenup = zeros(ComplexF64,Msize0,specnum,NOmega,Ng)

    # mat0 = matho + gdown*matdowndown + gup*matupup + gdu*matdownup + 1im*ksoc*matsoc
    mat0 = matho + 1im*ksoc*matsoc
    # mat1 = spzeros(ComplexF64,maxmatpcut+maxmatpcut2*2+maxmatpcut,maxmatpcut+maxmatpcut2*2+maxmatpcut)
    mat1 = copy(mat0)
    # matint = copy(mat0)
    phi = zeros(ComplexF64,maxmatpcut*2+maxmatpcut2*2,specnum)

    # println("diagonalising the Hamiltonian for different Omega ...")

    for jjg = 1:Ng

        gdujjg = arraygdu[jjg]
        # gupjjg = gdownjjg

        mat0 .= matho + gdownup*matdowndown + gdownup*matupup + gdujjg*matdownup + 1im*ksoc*matsoc
        # matint .= gdownup*matdowndown + gdownup*matupup + gdujjg*matdownup

        println("jjg=",jjg)

        for jj = 1:NOmega # parfor
            @time begin

                mat1 .= mat0 + arrayOmega[jj]*matW
                arraylambda[:,jj,jjg], phi = eigs(mat1,nev=specnum,which=:SR)
                arrayspect[:,jj,jjg] .= arraylambda[2:end,jj,jjg] .- arraylambda[1,jj,jjg]

                # onebodydensitymatrix!(Msize0,Np,phi[:,1],rhoij)
                # lambdaconden, phiconden = eigs(rhoij,nev=specnum,which=:LR)

                rhoij .= 0.0
                rhoijdown .= 0.0
                rhoijup .= 0.0
                onebodydensitymatrix!(Msize0,Np,phi[:,1],rhoij,rhoijdown,rhoijup)
                lambdacondendown, phicondendown = eigs(rhoijdown,nev=specnum,which=:LR)
                lambdacondenup, phicondenup = eigs(rhoijup,nev=specnum,which=:LR)

                arraylambdacondendown[:,jj,jjg] = lambdacondendown
                arraylambdacondenup[:,jj,jjg] = lambdacondenup
                arrayphicondendown[:,:,jj,jjg] = phicondendown
                arrayphicondenup[:,:,jj,jjg] = phicondenup

                # arraypopdown3[:,jj,jjg] = sum(abs.(phi[1:maxmatpcut,:]).^2,dims=1)'
                # arraypopdown2up1[:,jj,jjg] = sum(abs.(phi[maxmatpcut+1:maxmatpcut+maxmatpcut2,:]).^2,dims=1)'
                # arraypopdown1up2[:,jj,jjg] = sum(abs.(phi[maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2,:]).^2,dims=1)'
                # arraypopup3[:,jj,jjg] = sum(abs.(phi[maxmatpcut+maxmatpcut2+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2+maxmatpcut,:]).^2,dims=1)'

                # arraypopdown3GSspatial[:,jj,jjg] = phi[1:maxmatpcut,1]
                # arraypopdown2up1GSspatial[:,jj,jjg] = phi[maxmatpcut+1:maxmatpcut+maxmatpcut2,1]
                # arraypopdown1up2GSspatial[:,jj,jjg] = phi[maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2,1]
                # arraypopup3GSspatial[:,jj,jjg] = phi[maxmatpcut+maxmatpcut2+maxmatpcut2+1:end,1]

                # energy
                # arrayenergyGStot[jj,jjg] = phi[:,1]'*mat1*phi[:,1] # same as arraylambda[1]
                # arrayenergyGSint[jj,jjg] = phi[:,1]'*matint*phi[:,1]

                println("jj=",jj)

            end
        end

    end

    indksoc = Int64(ksoc)
    save("data_spectrum_onebody_gdu_jjg_ksoc$indksoc.jld", "arrayOmega", arrayOmega, "arraygdu", arraygdu, "ksoc", ksoc, "arraylambda", arraylambda, "arrayspect", arrayspect, "arraylambdacondendown", arraylambdacondendown, "arraylambdacondenup", arraylambdacondenup, "arrayphicondendown", arrayphicondendown, "arrayphicondenup", arrayphicondenup)

    # return arrayOmega, arraygdu, ksoc, arraylambda, arrayspect, arraypopdown3, arraypopdown2up1, arraypopdown1up2, arraypopup3

end

function diagonaliseH_paircorrelation_test(Msize0::Int64, Np::Int64, gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega::Float64, specnum::Int64, Lx::Float64, Nx::Int64)

    # construct Hamiltonian
    if Msize0 == 50
       matho=load("data_Htot$Msize0.Np3.jld")["matho"]
       matdowndown=load("data_Htot$Msize0.Np3.jld")["matdowndown"]
       matdownup=load("data_Htot$Msize0.Np3.jld")["matdownup"]
       matupup=load("data_Htot$Msize0.Np3.jld")["matupup"]
       matsoc=load("data_Htot$Msize0.Np3.jld")["matsoc"]
       matW=load("data_Htot$Msize0.Np3.jld")["matW"]
    else
       matho, matdowndown, matupup, matdownup, matsoc, matW = createHtotal(Msize0,Np)
    end

    # for down3
    Enecutoff = Msize0 - 1 + Np/2
    matp = zeros(Int64,Msize0+1,Np+1)
    pascaltriangle!(Msize0,Np,matp) # note the indices are m+1 and n+1 for N^m_n
    indvec = cutMsizeEnespinless(Msize0,Np,matp,Enecutoff)
    maxmatpcut = length(indvec)

    # for down2up1
    matp20 = zeros(Int64,Msize0+1,Np-1+1) # Np-1=2
    matp21 = zeros(Int64,Msize0+1,1+1)
    pascaltriangle!(Msize0,Np-1,matp20)
    pascaltriangle!(Msize0,1,matp21)
    indvec2 = cutMsizeEnespinmixed(Msize0,Np,matp20,matp21,Enecutoff,1)
    maxmatpcut2 = length(indvec2)

    arraylambda = zeros(ComplexF64,specnum)
    arrayspect = zeros(ComplexF64,specnum-1)

    psi = zeros(ComplexF64,maxmatpcut*2+maxmatpcut2*2,specnum)
    xrange = LinRange(-Lx,Lx,Nx)
    yrange = LinRange(-Lx,Lx,Nx)

    indgdown = gdown #Int64(gdown)
    indgup = gup #Int64(gup)
    indgdu = gdu #Int64(gdu)
    indksoc = ksoc #Int64(ksoc)
    indOmega = Omega #Int64(Omega)
    indNx = Nx

    println("diagonalising the Hamiltonian ...")
    @time begin
        # mattot = Hermitian(Array(matho + gdown*matdowndown + gup*matupup + gdu*matdownup + 1im*ksoc*matsoc + Omega*matW))
        mattot = matho + gdown*matdowndown + gup*matupup + gdu*matdownup + 1im*ksoc*matsoc + Omega*matW
        arraylambda, psi = eigs(mattot,nev=specnum,which=:SR)
        arrayspect .= arraylambda[2:end] .- arraylambda[1]
    end

    # energy
    mat0 = matho + 1im*ksoc*matsoc
    energyEk = psi[:,1]'*mat0*psi[:,1] + ksoc^2/2*3
    println("kinetic energy is", energyEk)


    println("calculating pair correlation ...")
    @time begin
        fun_nudown, fun_nudu, fun_nuup = paircorrelation_fun(indvec,indvec2,Msize0,Np,matp,matp20,matp21,psi[:,1],xrange,yrange)
    end

    fun = fun_nudown + fun_nudu + fun_nuup

    return arraylambda, arrayspect, xrange, yrange, fun_nudown, fun_nudu, fun_nuup, fun, psi[:,1]
    # save("data_paircorre_gdown$indgdown.gup$indgup.gdu$indgdu.ksoc$indksoc.Omega$indOmega.Nx$indNx.jld", "arraylambda", arraylambda, "arrayspect", arrayspect, "xrange", xrange, "yrange", yrange, "fun_nudown", fun_nudown, "fun_nudu", fun_nudu, "fun_nuup", fun_nuup)

end

function diagonaliseH_paircorrelation_arrayOmega(Msize0::Int64, Np::Int64, gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega::Float64, specnum::Int64, Lx::Float64, Nx::Int64, Omega0::Float64, Omega1::Float64, NOmega::Int64)

    arrayOmega = LinRange(Omega0,Omega1,NOmega)
    # arrayOmega = [10.0, ]
    # NOmega = length(arrayOmega)
    
    xrange = LinRange(-Lx,Lx,Nx)
    yrange = LinRange(-Lx,Lx,Nx)
    fun_nudown_Omega = zeros(Float64,Nx,Nx,NOmega)
    fun_nudu_Omega = zeros(Float64,Nx,Nx,NOmega)
    fun_nuup_Omega = zeros(Float64,Nx,Nx,NOmega)
    fun_Omega = zeros(Float64,Nx,Nx,NOmega)

    for jj = 1:NOmega
        arraylambda, arrayspect, xrange, yrange, fun_nudown_Omega[:,:,jj], fun_nudu_Omega[:,:,jj], fun_nuup_Omega[:,:,jj], fun_Omega[:,:,jj], psi = diagonaliseH_paircorrelation_test(Msize0,Np,gdown,gup,gdu,ksoc,arrayOmega[jj],specnum,Lx,Nx)
        println("jj=",jj)
    end

    return xrange, yrange, fun_nudown_Omega, fun_nudu_Omega, fun_nuup_Omega, fun_Omega

end

function diagonaliseH_paircorrelation(Msize0::Int64, Np::Int64, gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega::Float64, specnum::Int64, Lx::Float64, Nx::Int64)

    # construct Hamiltonian
    if Msize0 == 50
       matho=load("data_Htot$Msize0.Np3.jld")["matho"]
       matdowndown=load("data_Htot$Msize0.Np3.jld")["matdowndown"]
       matdownup=load("data_Htot$Msize0.Np3.jld")["matdownup"]
       matupup=load("data_Htot$Msize0.Np3.jld")["matupup"]
       matsoc=load("data_Htot$Msize0.Np3.jld")["matsoc"]
       matW=load("data_Htot$Msize0.Np3.jld")["matW"]
    else
       matho, matdowndown, matupup, matdownup, matsoc, matW = createHtotal(Msize0,Np)
    end
    # matho=load("data_Htot90_Np3.jld")["matho"]
    # matdowndown=load("data_Htot90_Np3.jld")["matdowndown"]
    # matdownup=load("data_Htot90_Np3.jld")["matdownup"]
    # matupup=load("data_Htot90_Np3.jld")["matupup"]
    # matsoc=load("data_Htot90_Np3.jld")["matsoc"]
    # matW=load("data_Htot90_Np3.jld")["matW"]

    # for down3
    Enecutoff = Msize0 - 1 + Np/2
    matp = zeros(Int64,Msize0+1,Np+1)
    pascaltriangle!(Msize0,Np,matp) # note the indices are m+1 and n+1 for N^m_n
    indvec = cutMsizeEnespinless(Msize0,Np,matp,Enecutoff)
    maxmatpcut = length(indvec)

    # for down2up1
    matp20 = zeros(Int64,Msize0+1,Np-1+1) # Np-1=2
    matp21 = zeros(Int64,Msize0+1,1+1)
    pascaltriangle!(Msize0,Np-1,matp20)
    pascaltriangle!(Msize0,1,matp21)
    indvec2 = cutMsizeEnespinmixed(Msize0,Np,matp20,matp21,Enecutoff,1)
    maxmatpcut2 = length(indvec2)

    arraylambda = zeros(ComplexF64,specnum)
    arrayspect = zeros(ComplexF64,specnum-1)

    psi = zeros(ComplexF64,maxmatpcut*2+maxmatpcut2*2,specnum)
    xrange = LinRange(-Lx,Lx,Nx)
    yrange = LinRange(-Lx,Lx,Nx)

    indgdown = gdown #Int64(gdown)
    indgup = gup #Int64(gup)
    indgdu = gdu #Int64(gdu)
    indksoc = ksoc #Int64(ksoc)
    indOmega = Omega #Int64(Omega)
    indNx = Nx

    println("diagonalising the Hamiltonian ...")
    @time begin
        # mattot = Hermitian(Array(matho + gdown*matdowndown + gup*matupup + gdu*matdownup + 1im*ksoc*matsoc + Omega*matW))
        mattot = matho + gdown*matdowndown + gup*matupup + gdu*matdownup + 1im*ksoc*matsoc + Omega*matW
        arraylambda, psi = eigs(mattot,nev=specnum,which=:SR)
        arrayspect .= arraylambda[2:end] .- arraylambda[1]
    end

    println("calculating pair correlation ...")
    @time begin
        fun_nudown, fun_nudu, fun_nuup = paircorrelation_fun(indvec,indvec2,Msize0,Np,matp,matp20,matp21,psi[:,1],xrange,yrange)
    end

    fun_nu = fun_nudown + fun_nudu + fun_nuup

    # return arraylambda, arrayspect, xrange, yrange, fun_nudown, fun_nudu, fun_nuup
    save("data_paircorre_gdown$indgdown.gup$indgup.gdu$indgdu.ksoc$indksoc.Omega$indOmega.Nx$indNx.jld", "arraylambda", arraylambda, "arrayspect", arrayspect, "xrange", xrange, "yrange", yrange, "fun_nudown", fun_nudown, "fun_nudu", fun_nudu, "fun_nuup", fun_nuup, "fun_nu", fun_nu)

end

function diagonaliseH_paircorrelation_arrayOmegag12(Msize0::Int64, Np::Int64, gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega0::Float64, Omega1::Float64, NOmega::Int64, specnum::Int64, Lx::Float64, Nx::Int64)

    # construct Hamiltonian
    matho=load("data_Htot90_Np3.jld")["matho"]
    matdowndown=load("data_Htot90_Np3.jld")["matdowndown"]
    matdownup=load("data_Htot90_Np3.jld")["matdownup"]
    matupup=load("data_Htot90_Np3.jld")["matupup"]
    matsoc=load("data_Htot90_Np3.jld")["matsoc"]
    matW=load("data_Htot90_Np3.jld")["matW"]

    # for down3
    Enecutoff = Msize0 - 1 + Np/2
    matp = zeros(Int64,Msize0+1,Np+1)
    pascaltriangle!(Msize0,Np,matp) # note the indices are m+1 and n+1 for N^m_n
    indvec = cutMsizeEnespinless(Msize0,Np,matp,Enecutoff)
    maxmatpcut = length(indvec)

    # for down2up1
    matp20 = zeros(Int64,Msize0+1,Np-1+1) # Np-1=2
    matp21 = zeros(Int64,Msize0+1,1+1)
    pascaltriangle!(Msize0,Np-1,matp20)
    pascaltriangle!(Msize0,1,matp21)
    indvec2 = cutMsizeEnespinmixed(Msize0,Np,matp20,matp21,Enecutoff,1)
    maxmatpcut2 = length(indvec2)

    arraylambda = zeros(ComplexF64,specnum)
    # arrayspect = zeros(ComplexF64,specnum-1)

    rhoij = zeros(ComplexF64,Msize0*2,Msize0*2)
    rhoijdown = zeros(ComplexF64,Msize0,Msize0)
    rhoijup = zeros(ComplexF64,Msize0,Msize0)
    lambdaconden = zeros(ComplexF64,specnum)
    phiconden = zeros(ComplexF64,Msize0*2,specnum)

    psi = zeros(ComplexF64,maxmatpcut*2+maxmatpcut2*2,specnum)
    xrange = LinRange(-Lx,Lx,Nx)
    yrange = LinRange(-Lx,Lx,Nx)

    indgdown = Int64(gdown)
    indgup = Int64(gup)
    indgdu = Int64(gdu)
    indksoc = Int64(ksoc)
    # indOmega = Int64(Omega)
    indNx = Nx

    arrayOmega = LinRange(Omega0,Omega1,NOmega)
    mattot = copy(matho + gdown*matdowndown + gup*matupup + gdu*matdownup + 1im*ksoc*matsoc)
    fun_nudown_Omega = zeros(Float64,Nx,Nx,NOmega)
    fun_nudu_Omega = zeros(Float64,Nx,Nx,NOmega)
    fun_nuup_Omega = zeros(Float64,Nx,Nx,NOmega)

    for jjOmega = 1:NOmega

        @time begin

        Omega = arrayOmega[jjOmega]

        mattot .= matho + gdown*matdowndown + gup*matupup + gdu*matdownup + 1im*ksoc*matsoc + Omega*matW
        arraylambda, psi = eigs(mattot,nev=specnum,which=:SR)
        # arrayspect .= arraylambda[2:end] .- arraylambda[1]

        fun_nudown, fun_nudu, fun_nuup = paircorrelation_fun(indvec,indvec2,Msize0,Np,matp,matp20,matp21,psi[:,1],xrange,yrange)

        fun_nudown_Omega[:,:,jjOmega] .= fun_nudown
        fun_nudu_Omega[:,:,jjOmega] .= fun_nudu
        fun_nuup_Omega[:,:,jjOmega] .= fun_nuup

        end

    end

    # return arraylambda, arrayspect, xrange, yrange, fun_nudown, fun_nudu, fun_nuup
    save("data_paircorre_Omega_gdown$indgdown.gup$indgup.gdu$indgdu.ksoc$indksoc.Nx$indNx.jld", "xrange", xrange, "yrange", yrange, "arrayOmega", arrayOmega, "fun_nudown_Omega", fun_nudown_Omega, "fun_nudu_Omega", fun_nudu_Omega, "fun_nuup_Omega", fun_nuup_Omega)

end

function diagonaliseH_paircorrelation_arrayOmegag12_parfor(Msize0::Int64, Np::Int64, gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega0::Float64, Omega1::Float64, NOmega::Int64, specnum::Int64, Lx::Float64, Nx::Int64)

    # construct Hamiltonian
    matho=load("data_Htot90_Np3.jld")["matho"]
    matdowndown=load("data_Htot90_Np3.jld")["matdowndown"]
    matdownup=load("data_Htot90_Np3.jld")["matdownup"]
    matupup=load("data_Htot90_Np3.jld")["matupup"]
    matsoc=load("data_Htot90_Np3.jld")["matsoc"]
    matW=load("data_Htot90_Np3.jld")["matW"]

    # for down3
    Enecutoff = Msize0 - 1 + Np/2
    matp = zeros(Int64,Msize0+1,Np+1)
    pascaltriangle!(Msize0,Np,matp) # note the indices are m+1 and n+1 for N^m_n
    indvec = cutMsizeEnespinless(Msize0,Np,matp,Enecutoff)
    maxmatpcut = length(indvec)

    # for down2up1
    matp20 = zeros(Int64,Msize0+1,Np-1+1) # Np-1=2
    matp21 = zeros(Int64,Msize0+1,1+1)
    pascaltriangle!(Msize0,Np-1,matp20)
    pascaltriangle!(Msize0,1,matp21)
    indvec2 = cutMsizeEnespinmixed(Msize0,Np,matp20,matp21,Enecutoff,1)
    maxmatpcut2 = length(indvec2)

    arraylambda = zeros(ComplexF64,specnum)
    # arrayspect = zeros(ComplexF64,specnum-1)

    rhoij = zeros(ComplexF64,Msize0*2,Msize0*2)
    rhoijdown = zeros(ComplexF64,Msize0,Msize0)
    rhoijup = zeros(ComplexF64,Msize0,Msize0)
    lambdaconden = zeros(ComplexF64,specnum)
    phiconden = zeros(ComplexF64,Msize0*2,specnum)

    psi = zeros(ComplexF64,maxmatpcut*2+maxmatpcut2*2,specnum)
    xrange = LinRange(-Lx,Lx,Nx)
    yrange = LinRange(-Lx,Lx,Nx)

    indgdown = Int64(gdown)
    indgup = Int64(gup)
    indgdu = Int64(gdu)
    indksoc = Int64(ksoc)
    # indOmega = Int64(Omega)
    indNx = Nx

    arrayOmega = LinRange(Omega0,Omega1,NOmega)
    # mattot = spzeros(ComplexF64,maxmatpcut*2+maxmatpcut2*2,maxmatpcut*2+maxmatpcut2*2) #copy(matho + gdown*matdowndown + gup*matupup + gdu*matdownup + 1im*ksoc*matsoc)
    fun_nudown_Omega = zeros(Float64,Nx,Nx,NOmega)
    fun_nudu_Omega = zeros(Float64,Nx,Nx,NOmega)
    fun_nuup_Omega = zeros(Float64,Nx,Nx,NOmega)

    Threads.@threads for jjOmega = 1:NOmega

        # @time begin

        # Omega = arrayOmega[jjOmega]

        mattot = spzeros(ComplexF64,maxmatpcut*2+maxmatpcut2*2,maxmatpcut*2+maxmatpcut2*2)

        mattot .= matho + gdown*matdowndown + gup*matupup + gdu*matdownup + 1im*ksoc*matsoc + arrayOmega[jjOmega]*matW
        arraylambda, psi = eigs(mattot,nev=specnum,which=:SR)
        # arrayspect .= arraylambda[2:end] .- arraylambda[1]

        fun_nudown, fun_nudu, fun_nuup = paircorrelation_fun(indvec,indvec2,Msize0,Np,matp,matp20,matp21,psi[:,1],xrange,yrange)

        fun_nudown_Omega[:,:,jjOmega] .= fun_nudown
        fun_nudu_Omega[:,:,jjOmega] .= fun_nudu
        fun_nuup_Omega[:,:,jjOmega] .= fun_nuup

        # end

    end

    # return arraylambda, arrayspect, xrange, yrange, fun_nudown, fun_nudu, fun_nuup
    save("data_paircorre_Omega_parfor_gdown$indgdown.gup$indgup.gdu$indgdu.ksoc$indksoc.Nx$indNx.jld", "xrange", xrange, "yrange", yrange, "arrayOmega", arrayOmega, "fun_nudown_Omega", fun_nudown_Omega, "fun_nudu_Omega", fun_nudu_Omega, "fun_nuup_Omega", fun_nuup_Omega)

end

function TG_paircorrelation_Np3(Nx::Int64,Lx::Float64,xi0::Float64)

    # funTG = zeros(Float64,Nx,Nx)

    # xrange = LinRange(-Lx,Lx,Nx)
    # yrange = LinRange(-Lx,Lx,Nx)

    # xi0 = 2*pi/10;
    xrange = LinRange(-Lx,Lx,Nx)
    dx = xrange[2]-xrange[1]

    psi0 = 1.0/sqrt(sqrt(pi))*exp.(-(xrange./xi0).^2/2)
    psi1 = 1.0/sqrt(sqrt(pi))*exp.(-(xrange./xi0).^2/2)/sqrt(2).*(xrange./xi0)*2
    psi2 = 1.0/sqrt(sqrt(pi))*exp.(-(xrange./xi0).^2/2)/sqrt(8).*(4*(xrange./xi0).^2 .- 2)

    funTG = psi0.^2*(psi1.^2)' + psi1.^2*(psi2.^2)' + psi2.^2*(psi0.^2)' +
            psi1.^2*(psi0.^2)' + psi2.^2*(psi1.^2)' + psi0.^2*(psi2.^2)' -
            2*(psi0.*psi1)*(psi0.*psi1)' -
            2*(psi1.*psi2)*(psi1.*psi2)' -
            2*(psi0.*psi2)*(psi0.*psi2)'

    # funTG .= funTG/(sum(funTG)*dx^2)
    funTG .= funTG/3/2/xi0^2

    return funTG

end

function sweepingxi0_paircorrelation(Nx::Int64,Lx::Float64,fun_nu::Matrix{Float64},xi0_0::Float64,xi0_1::Float64,Nxi0::Int64)

    xrange = LinRange(-Lx,Lx,Nx)
    dx = xrange[2]-xrange[1]

    arrayxi0 = LinRange(xi0_0,xi0_1,Nxi0)
    diff = zeros(Float64,Nxi0)

    fun_nu .= fun_nu/(sum(fun_nu)*dx^2)

    for jj = 1:Nxi0
        xi0 = arrayxi0[jj]
        funTG = TG_paircorrelation_Np3(Nx,Lx,xi0)
        funTG .= funTG/(sum(funTG)*dx^2)
        diff[jj] = sum(abs.(fun_nu - funTG))*dx^2
    end

    return arrayxi0, diff

end

function maxk(a::Vector{Float64}, k::Int64)
    b = partialsortperm(a, 1:k, rev=true)
    return collect(zip(b, a[b]))
end

function max2k(a::Matrix, k::Int64)

    if length(size(a)) != 2
       error("input should be 2d matrix")
    end

    length_v = length(a)
    v = zeros(Float64,length_v)
    size_v = size(a)

    for jj = 1:length_v
        jj1 = mod1(jj,size_v[1])
        jj2 = (jj-1)÷size_v[1]+1
        v[jj] = a[jj1,jj2]
    end

    b = partialsortperm(v, 1:k, rev=true)
    ind_b = zeros(Int64,k,2)

    for jj = 1:k
        jj1 = mod1(b[jj],size_v[1])
        jj2 = (b[jj]-1)÷size_v[1]+1
        ind_b[jj,:] = [jj1,jj2]
    end

    return v[b], ind_b

end

function charapara_widthdepth(fun_nu,Nx,xrange)

    peak2 = maxk(fun_nu[:,Int64((Nx-1)/2)+1],2)
    width = 0.0
    depth = 0.0
    charapara = 0.0

    if peak2[1][2] - peak2[2][2] < 10^(-5)
       width = abs(xrange[peak2[1][1]] - xrange[peak2[2][1]])
       depth = peak2[1][2]-minimum(fun_nu[min(peak2[1][1],peak2[2][1]):max(peak2[1][1],peak2[2][1]),Int64((Nx-1)/2)+1])
       charapara = width*depth
    end

    return charapara

end

function charapara_width(fun_nu,Nx,xrange)

    peak2 = maxk(fun_nu[:,Int64((Nx-1)/2)+1],2)
    width = 0.0

    if peak2[1][2] - peak2[2][2] < 10^(-5)
       width = abs(xrange[peak2[1][1]] - xrange[peak2[2][1]])
    else
       println("there is only one peak")
       width = 0.0
    end

    return width

end
