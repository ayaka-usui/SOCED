using Arpack, SparseArrays, LinearAlgebra
using Plots

# define functions used here
include("ades.jl")
include("acre.jl")
include("pascaltriangle.jl")
include("in2b.jl")
# include("b2in.jl")
include("epsilon.jl")
include("cutMsizeEne.jl")
include("Hintfunccutoff.jl")
include("delta.jl")

function Hsocfunccutoff!(indvec::Vector{Int64}, Msize0::Int64, Np::Int64, matp::Matrix{Int64}, energyijHOmat::Vector{Float64}, energyijSOCmat::Matrix{ComplexF64}, energyijWmat::Matrix{Float64}, ksoc::Float64, Omega::Float64, matHO::SparseMatrixCSC{Float64}, matSOC::SparseMatrixCSC{ComplexF64}, matW::SparseMatrixCSC{Float64})

    Msize = Msize0*2
    maxmatpcut = length(indvec)

    # matp = zeros(Int,Msize+1,Np+1);
    # pascaltriangle!(Msize,Np,matp) # the size is Msize+1 times Np+1
    # maxmatp = matp[Msize+1,Np+1] # the indices are m+1 and n+1 for N^m_ns

    # defines vectors and matrices
    # Hsoc .= 0; #Hsoc = spzeros(ComplexF64,maxmatpcut,maxmatpcut); #spzeros(ComplexF64,maxmatp,maxmatp);
    matHO .= 0.0
    matSOC .= 0.0
    matW .= 0.0

    # define a matrix for the Hamiltonian
    for nn = 1:maxmatpcut #maxmatp

        vecmbnn = spzeros(Int64,Msize+1)
        in2b!(indvec[nn],Msize,Np,matp,vecmbnn) #vecmbnn = in2b(indvec[nn],Msize,Np,matp) #in2b(nn,Msize,Np)
        energyij = 0.

        # free terms
        for jj = 1:Msize

            vecmbnnj = spzeros(Int64,Msize+1)
            ades!(jj,vecmbnn,vecmbnnj) #vecmbnnj = ades(jj,vecmbnn)
            if vecmbnnj[Msize+1] == 0
               continue
            end

            for ii = 1:Msize

                vecmbnnij = spzeros(Int64,Msize+1)
                acre!(ii,vecmbnnj,vecmbnnij) #vecmbnnij = acre(ii,vecmbnnj)

                # energyij = epsilon(ii,jj,ksoc,Omega)
                # energyij = energyijmat[ii,jj]
                energyij = energyijHOmat[ii]*delta(ii-jj) + energyijSOCmat[ii,jj] + energyijWmat[ii,jj]
                if isapprox(energyij,0) #energyij == 0 # energy is Int and zero if ii==jj
                   continue #vecmbnnij[1+Msize] = 0
                end

                for mm = 1:nn

                    vecmbmm = spzeros(Int64,Msize+1)
                    in2b!(indvec[mm],Msize,Np,matp,vecmbmm) #vecmbmm = in2b(indvec[mm],Msize,Np)

                    if vecmbnnij[1:Msize] == vecmbmm[1:Msize]

                       # Hsoc[mm,nn] = Hsoc[mm,nn] + energyij*sqrt(vecmbnnij[Msize+1])
                       matHO[mm,nn] = matHO[mm,nn] + energyijHOmat[ii]*delta(ii-jj)*sqrt(vecmbnnij[Msize+1])
                       matSOC[mm,nn] = matSOC[mm,nn] + energyijSOCmat[ii,jj]*sqrt(vecmbnnij[Msize+1])
                       matW[mm,nn] = matW[mm,nn] + energyijWmat[ii,jj]*sqrt(vecmbnnij[Msize+1])

                    end

                end

            end

        end

        # # Interactions
        # for ll = 1:Msize
        #
        # end

    end

    # since the Hamiltonian is helmitian
    # Hsoc .= Hsoc + Hsoc' - spdiagm(diag(Hsoc))
    # matHO .=
    matSOC .= matSOC + matSOC' - spdiagm(diag(matSOC))
    matW .= matW + matW' - spdiagm(diag(matW))

    # return Hsoc

end

function main(gdown::Float64, gup::Float64, gdu::Float64, ksoc::Float64, Omega0::Float64, Omega1::Float64, NOmega::Int64, Msize0::Int64, Np::Int64, Ene0minumhalf::Int64, specnum::Int64)

    # define cutoff of energy in Fock states
    Msize = Msize0*2
    matp = zeros(Int64,Msize+1,Np+1)
    pascaltriangle!(Msize,Np,matp) # the size is Msize+1 times Np+1
    # note the indices are m+1 and n+1 for N^m_n
    indvec = cutMsizeEne(Msize0,Np,matp,Ene0minumhalf)
    maxmatpcut = length(indvec)

    # interaction Hamiltonian
    # mat1, mat2, mat3 = Hintfunccutoff(indvec,Msize0,Np)
    mat1 = spzeros(Float64,maxmatpcut,maxmatpcut)
    mat2 = spzeros(Float64,maxmatpcut,maxmatpcut)
    mat3 = spzeros(Float64,maxmatpcut,maxmatpcut)
    # Hintfunccutoff!(indvec,Msize0,Np,matp,mat1,mat2,mat3)
    # lambda, phi = eigs(mat0+gdown*mat1+gup*mat2+gdu*mat3,nev=specnum,which=:SR)

    # single-particle Hamiltonian
    # calculate and save energyij
    energyijmat = zeros(Float64,Msize,Msize);
    energyijHOmat = zeros(Float64,Msize);
    energyijSOCmat = zeros(ComplexF64,Msize,Msize);
    energyijWmat = zeros(Float64,Msize,Msize);
    for ii = 1:Msize
        for jj = 1:ii

            # energyijmat[ii,jj] = epsilon(ii,jj,ksoc,Omega)

            if ii == jj
               energyijHOmat[ii] = epsilon(ii,jj,0.0,0.0)
            elseif iseven(ii+jj)
               energyijSOCmat[ii,jj] = epsilon(ii,jj,1.0,0.0)
            elseif abs(ii-jj) == 1 && ceil(Int64,ii/2) == ceil(Int64,jj/2) #nii == njj # && isodd(ii+jj)
               energyijWmat[ii,jj] = epsilon(ii,jj,0.0,1.0)
            end

        end
    end
    # energyijmat .= energyijmat + energyijmat' + spdiagm(diag(energyijmat))
    energyijSOCmat .= energyijSOCmat + energyijSOCmat' + spdiagm(diag(energyijSOCmat))
    energyijWmat .= energyijWmat + energyijWmat' + spdiagm(diag(energyijWmat))

    matHO = spzeros(Float64,maxmatpcut,maxmatpcut)
    matSOC = spzeros(ComplexF64,maxmatpcut,maxmatpcut)
    matW = spzeros(Float64,maxmatpcut,maxmatpcut)
    Hsocfunccutoff!(indvec,Msize0,Np,matp,energyijHOmat,energyijSOCmat,energyijWmat,1.0,1.0,matHO,matSOC,matW)

    # make an array for different Omega
    arrayOmega = LinRange(Omega0,Omega1,NOmega)
    lambdajj = zeros(Float64,specnum,NOmega)

    # different Omega
    for jj = 1:NOmega

        Omegajj = arrayOmega[jj]

        # diagonalise
        lambda, phi = eigs(matHO+matSOC*ksoc+matW*Omegajj+gdown*mat1+gup*mat2+gdu*mat3,nev=specnum,which=:SR)
        lambdajj[:,jj] .= real(lambda[:])

    end

    return lambdajj, arrayOmega


end
