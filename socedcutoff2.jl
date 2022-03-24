using Arpack, SparseArrays, LinearAlgebra
using JLD
# using ArnoldiMethod
# using FLoops

# define functions used here
include("pascaltriangle.jl")
include("cutMsizeEnespinless.jl")
include("cutMsizeEnespinmixed.jl")
include("cutMsizeEnespinmixed2.jl")
include("Hsocfunccutoffk1W1.jl")
include("Hintfunccutoff2.jl")

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

    return results

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

    return results

end

function saveHtot(Msize0::Int64, Np::Int64)

    matho, matdowndown, matupup, matdownup, matsoc, matW = createHtotal(Msize0,Np)

    # save
    save("data_Htot120_Np3.jld", "Msize0", Msize0, "matho", matho, "matdowndown", matdowndown, "matupup", matupup, "matdownup", matdownup, "matsoc", matsoc, "matW", matW)

end

function diagonalisesavedHtot(Msize0::Int64, Np::Int64, gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega::Float64, specnum::Int64)

    matho=load("data_Htot90_Np3.jld")["matho"]
    matdowndown=load("data_Htot90_Np3.jld")["matdowndown"]
    matdownup=load("data_Htot90_Np3.jld")["matdownup"]
    matupup=load("data_Htot90_Np3.jld")["matupup"]
    matsoc=load("data_Htot90_Np3.jld")["matsoc"]
    matW=load("data_Htot90_Np3.jld")["matW"]

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
    # mat1 = spzeros(ComplexF64,maxmatpcut+maxmatpcut2*2+maxmatpcut,maxmatpcut+maxmatpcut2*2+maxmatpcut)
    mat1 = copy(mat0)
    phi = zeros(ComplexF64,maxmatpcut*2+maxmatpcut2*2,specnum)

    println("diagonalising the Hamiltonian for different Omega ...")
    for jj = 1:NOmega # parfor
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

function diagonalisesavedHtotdiffW_cluster(Msize0::Int64, Np::Int64, gdown0::Float64, gdown1::Float64, Ng::Int64, gdu::Float64, ksoc::Float64, Omega0::Float64, Omega1::Float64, NOmega::Int64, specnum::Int64)

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

                arraypopdown3[:,jj,jjg] = sum(abs.(phi[1:maxmatpcut,:]).^2,dims=1)'
                arraypopdown2up1[:,jj,jjg] = sum(abs.(phi[maxmatpcut+1:maxmatpcut+maxmatpcut2,:]).^2,dims=1)'
                arraypopdown1up2[:,jj,jjg] = sum(abs.(phi[maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2,:]).^2,dims=1)'
                arraypopup3[:,jj,jjg] = sum(abs.(phi[maxmatpcut+maxmatpcut2+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2+maxmatpcut,:]).^2,dims=1)'

                println("jj=",jj)

            end
        end

    end

    save("data_spectrum_jjg.jld", "arrayOmega", arrayOmega, "arraylambda", arraylambda, "arrayspect", arrayspect, "arraypopdown3", arraypopdown3, "arraypopdown2up1", arraypopdown2up1, "arraypopdown1up2", arraypopdown1up2, "arraypopup3", arraypopup3)

    # return arrayOmega, arraylambda, arrayspect, arraypopdown3, arraypopdown2up1, arraypopdown1up2, arraypopup3

end
