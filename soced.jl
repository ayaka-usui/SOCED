using Arpack, SparseArrays, LinearAlgebra

# define functions used here
include("ades.jl")
include("acre.jl")
include("pascaltriangle.jl")
include("in2b.jl")
include("b2in.jl")
include("epsilon.jl")

function Hsoc(Msize0::Int64,Np::Int64,ksoc::Float64,Omega::Float64)

    Msize = Msize0*2
    matp = pascaltriangle(Msize,Np) # the size is Msize+1 times Np+1
    maxmatp = matp[Msize+1,Np+1] # the indices are m+1 and n+1 for N^m_ns

    # defines vectors and matrices
    vecmb = sparse(zeros(Int64,Msize+1));
    vecmbnn = sparse(zeros(Int64,Msize+1));
    Hsoc = sparse(zeros(Float64,maxmatp,maxmatp));

    # define a matrix for the Hamiltonian
    for nn = 1:maxmatp

        vecmbnn = in2b(nn,Msize,Np)
        energyijsum = 0.
        energyij = 0.

        for jj = 1:Msize

            vecmbnnj = ades(jj,vecmbnn)
            if vecmbnnj[Msize+1] == 0
               continue
            end

        for ii = 1:Msize

            vecmbnnij = acre(ii,vecmbnnj)

            energyij = epsilon(ii,jj,Msize0,ksoc,Omega)
            if isapprox(energyij,0) #energyij == 0 # energy is Int and zero if ii==jj
               continue #vecmbnnij[1+Msize] = 0
            end

            for mm = nn:maxmatp
                vecmbmm = in2b(mm,Msize,Np)
                Hsoc[mm,nn] = energyijsum
            end

        end


    end

    # since the Hamiltonian is elmitian
    Hsoc = Hsoc + Hsoc' - diagm(diag(Hsoc))

    return Hsoc

end
