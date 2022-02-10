include("in2bind.jl")
include("delta.jl")

function epsilonsoc(ii::Int64,jj::Int64,common::Int64,Np::Int64,ksoc::Float64)

    if abs(jj-ii) != 1
       return 0.0
    end

    if jj == common # a|2> = sqrt(2)|1>
       epsilpn = sqrt(Np)
    elseif ii == common # a^+|1> = sqrt(2)|2>
       epsilpn = sqrt(Np)
    else
       epsilpn = 1.0
    end

    njj = jj - 1
    nii = ii - 1
    epsilpn = epsilpn*(-1im)*ksoc/sqrt(2)*(sqrt(njj+1)*delta(nii-njj-1) - sqrt(njj)*delta(nii-njj+1)) # for down down

    return epsilpn

end

function epsilonsoctest(ii::Int64,jj::Int64,ksoc::Float64) # for down down
    njj = jj - 1
    nii = ii - 1
    return (-1im)*ksoc/sqrt(2)*(sqrt(njj+1)*delta(nii-njj-1) - sqrt(njj)*delta(nii-njj+1))
end

function Hsocfunccutoffk1W1!(indvec::Vector{Int64}, indvec2::Vector{Int64}, Msize0::Int64, Np::Int64, matp::Matrix{Int64}, matp2::Matrix{Int64}, Hho::SparseMatrixCSC{Float64}, Hsoc::SparseMatrixCSC{ComplexF64}, HW::SparseMatrixCSC{Float64})

    maxmatpcut = length(indvec)
    maxmatpcut2 = length(indvec2)

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

    # define vectors
    vecmbindnn = zeros(Int64,Np)
    vecmbindmm = zeros(Int64,Np)

    # define a matrix for the Hamiltonian for down down
    for nn = 1:maxmatpcut # parfor

        # ket
        in2bind!(indvec[nn],Msize0,Np,matp,vecmbindnn)

        # bra
        # in2bind!(indvec[mm],Msize0,Np,matp,vecmbindmm)
        vecmbindmm .= vecmbindnn

        # Hho
        if vecmbindnn[1] == vecmbindnn[2]
           Hho[nn,nn] = (vecmbindnn[1]-1+1/2)*Np
        else
           Hho[nn,nn] = (vecmbindnn[1]-1+1/2) + (vecmbindnn[2]-1+1/2)
        end

        for mm = 1:nn-1

            # bra
            in2bind!(indvec[mm],Msize0,Np,matp,vecmbindmm)

            # Hsoc
            if vecmbindnn[1] == vecmbindmm[1]

               common = vecmbindnn[1]
               jj = vecmbindnn[2]
               ii = vecmbindmm[2] # ii != jj
               Hsoc[mm,nn] = epsilonsoc(ii,jj,common,Np,1.0)

            elseif vecmbindnn[1] == vecmbindmm[2]

               common = vecmbindnn[1]
               jj = vecmbindnn[2]
               ii = vecmbindmm[1] # ii != jj
               Hsoc[mm,nn] = epsilonsoc(ii,jj,common,Np,1.0)

            elseif vecmbindnn[2] == vecmbindmm[1]

               common = vecmbindnn[2]
               jj = vecmbindnn[1]
               ii = vecmbindmm[2] # ii != jj
               Hsoc[mm,nn] = epsilonsoc(ii,jj,common,Np,1.0)

            elseif vecmbindnn[2] == vecmbindmm[2]

               common = vecmbindnn[2]
               jj = vecmbindnn[1]
               ii = vecmbindmm[1] # ii != jj
               Hsoc[mm,nn] = epsilonsoc(ii,jj,common,Np,1.0)

            end

        end

    end

    # define the Hamiltonian for up up
    Hho[end-(maxmatpcut-1):end,end-(maxmatpcut-1):end] = Hho[1:maxmatpcut,1:maxmatpcut]
    Hsoc[end-(maxmatpcut-1):end,end-(maxmatpcut-1):end] = Hsoc[1:maxmatpcut,1:maxmatpcut]*(-1) # due to up up

    # define the Hamiltonian for down up
    vecmbindnn2 = zeros(Int64,Np)

    for nn = 1:maxmatpcut2

        # ket
        indket2 = mod(indvec2[nn],Msize0)
        if indket2 != 0
           indket1 = div(indvec2[nn],Msize0)+1
        else
           indket1 = div(indvec2[nn],Msize0)
           indket2 = Msize0
        end
        in2bind!(indket1,Msize0,Np-1,matp2,vecmbindnn2)
        vecmbindnn[1] = vecmbindnn2[1]
        in2bind!(indket2,Msize0,Np-1,matp2,vecmbindnn2)
        vecmbindnn[2] = vecmbindnn2[1]

        # Hho
        Hho[maxmatpcut+nn,maxmatpcut+nn] = (vecmbindnn[1]-1+1/2) + (vecmbindnn[2]-1+1/2)

        for mm = 1:nn-1

            # bra
            indbra2 = mod(indvec2[mm],Msize0)
            if indbra2 != 0
               indbra1 = div(indvec2[mm],Msize0)+1
            else
               indbra1 = div(indvec2[mm],Msize0)
               indbra2 = Msize0
            end
            in2bind!(indbra1,Msize0,Np-1,matp2,vecmbindnn2)
            vecmbindmm[1] = vecmbindnn2[1]
            in2bind!(indbra2,Msize0,Np-1,matp2,vecmbindnn2)
            vecmbindmm[2] = vecmbindnn2[1]

            # Hsoc
            if vecmbindnn[1] == vecmbindmm[1]

               jj = vecmbindnn[2]
               ii = vecmbindmm[2] # ii != jj
               Hsoc[maxmatpcut+mm,maxmatpcut+nn] = epsilonsoc(ii,jj,0,Np,1.0)*(-1) # up up

            elseif vecmbindnn[2] == vecmbindmm[2]

               jj = vecmbindnn[1]
               ii = vecmbindmm[1] # ii != jj
               Hsoc[maxmatpcut+mm,maxmatpcut+nn] = epsilonsoc(ii,jj,0,Np,1.0) # down down

            end

        end

    end

    # since the Hamiltonian is helmitian
    Hsoc .= Hsoc + Hsoc' - spdiagm(diag(Hsoc))
    HW .= HW + HW' - spdiagm(diag(HW))

end
