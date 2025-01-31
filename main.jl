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
include("changebasis.jl")

#--- main code
# createHtotal: create the total Hamiltonian
# diagonaliseHtot... and diagonalisesavedHtot...: compute the spectrum and the spin population

#--- quick comment
# run "diagonaliseHtotsingle" to get the spectrum for a parameter set

######### create the total Hamiltonian

function createHtotal(Msize0::Int64, Np::Int64)

    # construct the total Hamiltonian

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

    # for down down up and up up down
    matp20 = zeros(Int64,Msize0+1,Np-1+1) # Np-1=2
    matp21 = zeros(Int64,Msize0+1,1+1)
    pascaltriangle!(Msize0,Np-1,matp20)
    pascaltriangle!(Msize0,1,matp21)
    indvec2 = cutMsizeEnespinmixed(Msize0,Np,matp20,matp21,Enecutoff,1)
    maxmatpcut2 = length(indvec2)
    
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

function saveHtot(Msize0::Int64, Np::Int64)

    # save data of the total Hamiltonian

    matho, matdowndown, matupup, matdownup, matsoc, matW = createHtotal(Msize0,Np)

    # save
    save("data_Htot$Msize0.Np3.jld", "Msize0", Msize0, "matho", matho, "matdowndown", matdowndown, "matupup", matupup, "matdownup", matdownup, "matsoc", matsoc, "matW", matW)

end

######### create the transformation that converts the basis for 2nd quantisation to the spin basis

function mat_from2ndtospin_downup(Msize0::Int64, Np::Int64)

    # convert down down up and up up down into spin basis
    # we do not touch down down down and up up up

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

    # for down down up and up up down
    # Enecutoff2 = Msize0 - 1 + (Np-1)/2
    # matp2 = zeros(Int64,Msize0+1,Np-1)
    matp20 = zeros(Int64,Msize0+1,Np-1+1) # Np-1=2
    matp21 = zeros(Int64,Msize0+1,1+1)
    pascaltriangle!(Msize0,Np-1,matp20)
    pascaltriangle!(Msize0,1,matp21)
    indvec2 = cutMsizeEnespinmixed(Msize0,Np,matp20,matp21,Enecutoff,1)
    maxmatpcut2 = length(indvec2)
    
    # construct transformation from 2nd quantisation basis to spin basis (only for downdownup and upupdown)
    mat_from2ndto1st_downup, mat_from1sttospin_downdownup, mat_from1sttospin_upupdown, mat_S2, mat_M2, mat_M4, mat_S3, mat_M1, mat_M3 = changefrom2ndtospin_downup(indvec2,Msize0,Np,matp20,matp21)
    
    # save
    # save("data_mat_from2ndtospin_downup_$Msize0.Np3.jld", "Msize0", Msize0, "Np", Np, "mat_from2ndto1st_downup", mat_from2ndto1st_downup, "mat_from1sttospin_downdownup", mat_from1sttospin_downdownup, "mat_from1sttospin_upupdown", mat_from1sttospin_upupdown, "mat_S2", mat_S2, "mat_M2", mat_M2, "mat_M4", mat_M4, "mat_S3", mat_S3, "mat_M1", mat_M1, "mat_M3", mat_M3, "maxmatpcut", maxmatpcut, "maxmatpcut2", maxmatpcut2)
    
    return mat_from2ndto1st_downup, mat_from1sttospin_downdownup, mat_from1sttospin_upupdown, mat_S2, mat_M2, mat_M4, mat_S3, mat_M1, mat_M3, maxmatpcut, maxmatpcut2

end

######### compute the spectrum and the spin population

function diagonaliseHtotsingle(Msize0::Int64, Np::Int64, gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega::Float64, specnum::Int64)

    # diagonalise the total Hamiltonian and extract the spectrum

    matho, matdowndown, matupup, matdownup, matsoc, matW = createHtotal(Msize0,Np)
    lambda, phi = eigs(matho + gdown*matdowndown + gup*matupup + gdu*matdownup + 1im*ksoc*matsoc + Omega*matW,nev=specnum,which=:SR)
    return lambda

end

function diagonaliseHtotspinpop_eigs(Msize0::Int64, Np::Int64, gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega::Float64, specnum::Int64)

    # diagonalise the total Hamiltonian and compute the spin population

    matho, matdowndown, matupup, matdownup, matsoc, matW = createHtotal(Msize0,Np)
    lambda, phi = eigs(matho + gdown*matdowndown + gup*matupup + gdu*matdownup + 1im*ksoc*matsoc + Omega*matW,nev=specnum,which=:SR)
    # lambda, phi = eigs(matho + gdown*matdowndown + gup*matupup + gdu*matdownup,nev=specnum,which=:SR)
    # lambda, phi = eigs(mattot,nev=specnum,which=:SR)
    # lambda, phi = eigen(mattot,1:specnum)

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
    return results

end

function diagonalisesavedHtot(Msize0::Int64, Np::Int64, gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega::Float64, specnum::Int64)

    # import data of the total Hamiltonian and compute the spin population

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

function diagonalisesavedHtotdiffW(Msize0::Int64, Np::Int64, gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega0::Float64, Omega1::Float64, NOmega::Int64, specnum::Int64)

    # import data of the total Hamiltonian and compute the spin population for different Omega

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

    # import data of the total Hamiltonian and compute the spin population for different Omega and different gdown

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

    # import data of the total Hamiltonian and compute the spin population for different Omega and different gdu

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

function population_SM(phi::Union{Vector{Float64},Vector{ComplexF64}}, mat_from2ndto1st_downup::Union{Matrix{Float64},SparseMatrixCSC{Float64}}, mat_from1sttospin_downdownup::Union{Matrix{Float64},SparseMatrixCSC{Float64}}, mat_from1sttospin_upupdown::Union{Matrix{Float64},SparseMatrixCSC{Float64}}, mat_S2::SparseMatrixCSC{Int64}, mat_M2::SparseMatrixCSC{Int64}, mat_M4::SparseMatrixCSC{Int64}, mat_S3::SparseMatrixCSC{Int64}, mat_M1::SparseMatrixCSC{Int64}, mat_M3::SparseMatrixCSC{Int64}, maxmatpcut::Int64, maxmatpcut2::Int64)

    # compute the population of the spin basis

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

######### compute pair correlation

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
        arraylambda, arrayspect, xrange, yrange, fun_nudown_Omega[:,:,jj], fun_nudu_Omega[:,:,jj], fun_nuup_Omega[:,:,jj], fun_Omega[:,:,jj], psi = diagonaliseH_paircorrelation(Msize0,Np,gdown,gup,gdu,ksoc,arrayOmega[jj],specnum,Lx,Nx)
        println("jj=",jj)
    end

    return xrange, yrange, fun_nudown_Omega, fun_nudu_Omega, fun_nuup_Omega, fun_Omega

end

function TG_paircorrelation_Np3(Nx::Int64,Lx::Float64,xi0::Float64)

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

######### plot

function plot_spinpop(arrayOmega,arraypopdown3, arraypopdown2up1, arraypopdown1up2, arraypopup3,arraypopS1,arraypopS2,arraypopS3,arraypopS4,arraypopM1,arraypopM2,arraypopM3,arraypopM4,jj1,jj2)

    # plot the spin population

    NOmega = length(arrayOmega)

    plot(arrayOmega,ones(NOmega)*3/4,color=:black,ls=:dot,lw=2)
    plot!(arrayOmega,ones(NOmega)*1/4,color=:black,ls=:dot,lw=2)
    # plot!([arrayOmega[jj2],arrayOmega[jj2]+1e-5],[-10,10],color=:black,lw=2,ls=:dot)

    plot!(arrayOmega,arraypopdown3[1,:,jj1]+arraypopup3[1,:,jj1],color=:grey,lw=4)
    plot!(arrayOmega,arraypopdown2up1[1,:,jj1]+arraypopdown1up2[1,:,jj1],color=:black,lw=4)
    # plot!(arrayOmega,arraypopM1[:,jj1]+arraypopM2[:,jj1]+arraypopM3[:,jj1]+arraypopM4[:,jj1],lw=3,color=:red,ls=:dot)
    # plot!(arrayOmega,arraypopS2[:,jj1]+arraypopS3[:,jj1],lw=3,ls=:dash,color=:blue)
    # plot!(arrayOmega,arraypopS1[:,jj1]+arraypopS4[:,jj1],lw=2)

    plot!(legend=:none)
    
    xlims!((9,41))
    ylims!((-0.05,1.05))
    plot!(aspect_ratio=32/1.1)

end

function plot_eneint(arrayOmega,arrayenergyGSint_g1,arrayenergyGSint_g2,arrayenergyGSint_g10)

    # plot the interaction energy

    plot(arrayOmega,real(arrayenergyGSint_g1[:,1]),lw=3)
    plot!(arrayOmega,real(arrayenergyGSint_g2[:,1]),lw=3,ls=:dash)
    plot!(arrayOmega,real(arrayenergyGSint_g10[:,1]),lw=3,ls=:dot)
    
    plot!(legend=:none)
    
    xlims!((9,41))
    ylims!((0.0,0.51))
    plot!(aspect_ratio=32/0.51)

end

function plot_arrayspect(arrayOmega,arrayspect,jj1,jj2)

    # plot the spectrum

    arrayspect1 = copy(arrayspect)
    jj3 = 167 #184 #129 #151 #179
    arrayspect1[3,jj3:end,jj1] = arrayspect[4,jj3:end,jj1]
    arrayspect1[4,jj3:end,jj1] = arrayspect[3,jj3:end,jj1]

    plot([2*4^2,2*4^2+1e-5],[-10,10],color=:black,lw=2,ls=:dash)
    plot!([arrayOmega[jj2],arrayOmega[jj2]+1e-5],[-10,10],color=:black,lw=2,ls=:dot)

    plot!(arrayOmega,real(arrayspect1[1,:,jj1]),lw=11,color=:black)
    plot!(arrayOmega,real(arrayspect1[2,:,jj1]),lw=9,color=RGB(8/255,104/255,172/255))
    plot!(arrayOmega,real(arrayspect1[3,:,jj1]),lw=7,color=RGB(67/255,162/255,202/255))
    plot!(arrayOmega,real(arrayspect1[4,:,jj1]),lw=6,color=RGB(123/255,204/255,196/255))
    plot!(arrayOmega,real(arrayspect1[5,:,jj1]),lw=4,color=RGB(186/255,228/255,188/255)) 

    plot!(legend=:none)
    
    xlims!((9,41))
    ylims!((-0.05,1.2))
    plot!(aspect_ratio=32/1.25)

end