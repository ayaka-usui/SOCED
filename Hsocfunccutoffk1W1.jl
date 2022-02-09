include("in2bind.jl")

# include("epsilon.jl")
# include("ades.jl")
# include("acre.jl")

function Hsocfunccutoffk1W1!(indvec::Vector{Int64}, Msize0::Int64, Np::Int64, matp::Matrix{Int64}, Hho::SparseMatrixCSC{Float64}, Hsoc::SparseMatrixCSC{ComplexF64}, HW::SparseMatrixCSC{Float64})

    maxmatpcut = length(indvec)

    # # calculate and save energyij
    # energyijmat = zeros(ComplexF64,Msize0,Msize0);
    # for ii = 1:Msize0
    #     for jj = 1:Msize0
    #         energyijmat[ii,jj] = epsilon(ii,jj,1.0,1.0)
    #     end
    # end

    # reset matrices
    Hho .= 0
    Hsoc .= 0
    HW .= 0

    # defines vectors
    vecmbindnn = zeros(Int64,Np)
    vecmbindmm = zeros(Int64,Np)

    # define a matrix for the Hamiltonian
    for nn = 1:maxmatpcut # parfor

        # ket
        in2bind!(indvec[nn],Msize0,Np,matp,vecmbindnn)

        # bra
        mm = nn
        in2bind!(indvec[mm],Msize0,Np,matp,vecmbindmm)

        if vecmbindnn[1] == vecmbindnn[2]
           Hho[nn,nn] = (vecmbindnn[1]-1+1/2)*Np
        else
           Hho[nn,nn] = (vecmbindnn[1]-1+1/2) + (vecmbindnn[2]-1+1/2)
        end

    end

    # since the Hamiltonian is helmitian
    Hsoc .= Hsoc + Hsoc' - spdiagm(diag(Hsoc))
    HW .= HW + HW' - spdiagm(diag(HW))

end
