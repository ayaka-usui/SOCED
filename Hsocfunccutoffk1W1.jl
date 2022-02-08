include("epsilon.jl")
include("in2b.jl")

include("ades.jl")
include("acre.jl")

function Hsocfunccutoffk1W1!(indvec::Vector{Int64}, Msize0::Int64, Np::Int64, matp::Matrix{Int64}, Hho::SparseMatrixCSC{Float64}, Hsoc::SparseMatrixCSC{ComplexF64}, HW::SparseMatrixCSC{Float64})

    maxmatpcut = length(indvec)

    # calculate and save energyij
    energyijmat = zeros(ComplexF64,Msize0,Msize0);
    for ii = 1:Msize0
        for jj = 1:Msize0
            energyijmat[ii,jj] = epsilon(ii,jj,1.0,1.0)
        end
    end

    # reset matrices
    Hho .= 0
    Hsoc .= 0
    HW .= 0

    # defines vectors
    # vecmbnn = spzeros(Int64,Msize0+1)
    # vecmbnnj = spzeros(Int64,Msize0+1)
    # vecmbnnij = spzeros(Int64,Msize0+1)
    # vecmbmm = spzeros(Int64,Msize0+1)

    # define a matrix for the Hamiltonian
    for nn = 1:maxmatpcut # parfor
        for mm = 1:nn

            in2b!(indvec[nn],Msize0,Np,matp,vecmbnn)

        end
    end









    for nn = 1:maxmatpcut #maxmatp # parfor

        vecmbnn .= spzeros(Int64,Msize+1)
        in2b!(indvec[nn],Msize,Np,matp,vecmbnn,1) #vecmbnn = in2b(indvec[nn],Msize,Np,matp) #in2b(nn,Msize,Np)
        # energyij = 0.

        # free terms
        for jj = 1:Msize

            vecmbnnj .= spzeros(Int64,Msize+1)
            ades!(jj,vecmbnn,vecmbnnj,1) #vecmbnnj = ades(jj,vecmbnn)
            if vecmbnnj[Msize+1] == 0 # go back if a|0> = 0
               continue
            end

            for ii = 1:Msize

                vecmbnnij .= spzeros(Int64,Msize+1)
                acre!(ii,vecmbnnj,vecmbnnij,1) #vecmbnnij = acre(ii,vecmbnnj)

                # energyij = epsilon(ii,jj,ksoc,Omega)
                energyij = energyijmat[ii,jj]
                if isapprox(energyij,0) #energyij == 0 # energy is Int and zero if ii==jj
                   continue #vecmbnnij[1+Msize] = 0
                end

                for mm = 1:nn

                    vecmbmm .= spzeros(Int64,Msize+1)
                    in2b!(indvec[mm],Msize,Np,matp,vecmbmm,1) #vecmbmm = in2b(indvec[mm],Msize,Np)

                    if vecmbnnij[1:Msize] == vecmbmm[1:Msize]

                       if ii == jj
                          Hho[mm,nn] = Hho[mm,nn] + energyij*sqrt(vecmbnnij[Msize+1])
                       elseif iseven(ii+jj) && abs(ii-jj) == 2   #abs(nii-njj) == 1
                          Hsoc[mm,nn] = Hsoc[mm,nn] + energyij*sqrt(vecmbnnij[Msize+1])
                       elseif abs(ii-jj) == 1 && ceil(Int64,ii/2) == ceil(Int64,jj/2) #nii == njj
                          HW[mm,nn] = HW[mm,nn] + energyij*sqrt(vecmbnnij[Msize+1])
                       end

                    end

                end

            end

        end

    end

    # since the Hamiltonian is helmitian
    # Hho .= Hho + Hho' - spdiagm(diag(Hho))
    Hsoc .= Hsoc + Hsoc' - spdiagm(diag(Hsoc))
    HW .= HW + HW' - spdiagm(diag(HW))

end
