using Arpack, SparseArrays, LinearAlgebra

# define functions used here
include("ades.jl")
include("acre.jl")
include("pascaltriangle.jl")
include("in2b.jl")
include("b2in.jl")
include("epsilon.jl")

function Hsocfunc(Msize0::Int64, Np::Int64, ksoc::Float64, Omega::Float64)

    Msize = Msize0*2

    matp = zeros(Int,Msize+1,Np+1);
    pascaltriangle!(Msize,Np,matp) # the size is Msize+1 times Np+1
    maxmatp = matp[Msize+1,Np+1] # the indices are m+1 and n+1 for N^m_ns

    # defines vectors and matrices
    vecmbnn = spzeros(Int64,Msize+1);
    vecmbnnj = spzeros(Int64,Msize+1);
    vecmbnnij = spzeros(Int64,Msize+1);
    Hsoc = spzeros(ComplexF64,maxmatp,maxmatp);

    # define a matrix for the Hamiltonian
    for nn = 1:maxmatp

        vecmbnn = in2b(nn,Msize,Np)
        energyij = 0.

        for jj = 1:Msize

            vecmbnnj = ades(jj,vecmbnn)
            if vecmbnnj[Msize+1] == 0
               continue
            end

            for ii = 1:Msize

                vecmbnnij = acre(ii,vecmbnnj)

                energyij = epsilon(ii,jj,ksoc,Omega)
                if isapprox(energyij,0) #energyij == 0 # energy is Int and zero if ii==jj
                   continue #vecmbnnij[1+Msize] = 0
                end

                for mm = 1:nn

                    vecmbmm = in2b(mm,Msize,Np)
                    if vecmbnnij[1:Msize] == vecmbmm[1:Msize]
                       Hsoc[mm,nn] = Hsoc[mm,nn] + energyij*sqrt(vecmbnnij[Msize+1])
                    end

                end

            end

        end

    end

    # since the Hamiltonian is helmitian
    Hsoc = Hsoc + Hsoc' - spdiagm(diag(Hsoc))

    return Hsoc

end

function diagonaliseHsoc(Hsoc::SparseMatrixCSC{ComplexF64})

    lambda, phi = eigs(Hsoc,nev=6,which=:SR)

    return lambda, phi

end
