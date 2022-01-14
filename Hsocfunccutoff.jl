include("ades.jl")
include("acre.jl")
# include("pascaltriangle.jl")
include("in2b.jl")
# include("b2in.jl")
include("epsilon.jl")

function Hsocfunccutoff!(indvec::Vector{Int64}, Msize0::Int64, Np::Int64, matp::Matrix{Int64}, ksoc::Float64, Omega::Float64, Hsoc::SparseMatrixCSC{ComplexF64})

    Msize = Msize0*2
    maxmatpcut = length(indvec)
    # matp = zeros(Int,Msize+1,Np+1);
    # pascaltriangle!(Msize,Np,matp) # the size is Msize+1 times Np+1
    # maxmatp = matp[Msize+1,Np+1] # the indices are m+1 and n+1 for N^m_ns

    # calculate and save energyij
    energyijmat = zeros(ComplexF64,Msize,Msize); #zeros(Float64,Msize,Msize);
    for ii = 1:Msize
        for jj = 1:Msize
            energyijmat[ii,jj] = epsilon(ii,jj,ksoc,Omega)
        end
    end

    # defines vectors and matrices
    Hsoc .= 0; #Hsoc = spzeros(ComplexF64,maxmatpcut,maxmatpcut); #spzeros(ComplexF64,maxmatp,maxmatp);
    vecmbnn = spzeros(Int64,Msize+1)
    vecmbnnj = spzeros(Int64,Msize+1)
    vecmbnnij = spzeros(Int64,Msize+1)
    vecmbmm = spzeros(Int64,Msize+1)

    # define a matrix for the Hamiltonian
    for nn = 1:maxmatpcut #maxmatp # parfor

        vecmbnn .= spzeros(Int64,Msize+1)
        in2b!(indvec[nn],Msize,Np,matp,vecmbnn) #vecmbnn = in2b(indvec[nn],Msize,Np,matp) #in2b(nn,Msize,Np)
        # energyij = 0.

        # free terms
        for jj = 1:Msize

            vecmbnnj .= spzeros(Int64,Msize+1)
            ades!(jj,vecmbnn,vecmbnnj) #vecmbnnj = ades(jj,vecmbnn)
            if vecmbnnj[Msize+1] == 0 # go back if a|0> = 0
               continue
            end

            for ii = 1:Msize

                vecmbnnij .= spzeros(Int64,Msize+1)
                acre!(ii,vecmbnnj,vecmbnnij) #vecmbnnij = acre(ii,vecmbnnj)

                # energyij = epsilon(ii,jj,ksoc,Omega)
                energyij = energyijmat[ii,jj]
                if isapprox(energyij,0) #energyij == 0 # energy is Int and zero if ii==jj
                   continue #vecmbnnij[1+Msize] = 0
                end

                for mm = 1:nn

                    vecmbmm .= spzeros(Int64,Msize+1)
                    in2b!(indvec[mm],Msize,Np,matp,vecmbmm) #vecmbmm = in2b(indvec[mm],Msize,Np)

                    if vecmbnnij[1:Msize] == vecmbmm[1:Msize]
                       Hsoc[mm,nn] = Hsoc[mm,nn] + energyij*sqrt(vecmbnnij[Msize+1])
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
    Hsoc .= Hsoc + Hsoc' - spdiagm(diag(Hsoc))

end
