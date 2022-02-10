include("in2bind.jl")
include("delta.jl")

function epsilonsoc(nii::Int64,njj::Int64,ksoc::Float64)
    return -1im*ksoc/sqrt(2)*(sqrt(njj+1)*delta(nii-njj-1) - sqrt(njj)*delta(nii-njj+1)) # for down down
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

    # defines vectors
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

            if vecmbindnn[1] == vecmbindmm[1]

               common = vecmbindnn[1]
               jj = vecmbindnn[2]
               ii = vecmbindmm[2] # ii != jj

               if jj == common # a|2> = sqrt(2)|1>
                  Hsoc[mm,nn] = epsilonsoc(ii,jj,1.0)*sqrt(Np) # ksoc=1.0
               elseif ii == common # a^+|1> = sqrt(2)|2>
                  Hsoc[mm,nn] = epsilonsoc(ii,jj,1.0)*sqrt(Np)
               else
                  Hsoc[mm,nn] = epsilonsoc(ii,jj,1.0)
               end

            elseif vecmbindnn[1] == vecmbindmm[2]

               common = vecmbindnn[1]
               jj = vecmbindnn[2]
               ii = vecmbindmm[1] # ii != jj

               if jj == common
                  Hsoc[mm,nn] = epsilonsoc(ii,jj,1.0)*sqrt(Np)
               elseif ii == common
                  Hsoc[mm,nn] = epsilonsoc(ii,jj,1.0)*sqrt(Np)
               else
                  Hsoc[mm,nn] = epsilonsoc(ii,jj,1.0)
               end

            elseif vecmbindnn[2] == vecmbindmm[1]

               common = vecmbindnn[2]
               jj = vecmbindnn[1]
               ii = vecmbindmm[2] # ii != jj

               if jj == common
                  Hsoc[mm,nn] = epsilonsoc(ii,jj,1.0)*sqrt(Np)
               elseif ii == common
                  Hsoc[mm,nn] = epsilonsoc(ii,jj,1.0)*sqrt(Np)
               else
                  Hsoc[mm,nn] = epsilonsoc(ii,jj,1.0)
               end

            elseif vecmbindnn[2] == vecmbindmm[2]

               common = vecmbindnn[2]
               jj = vecmbindnn[1]
               ii = vecmbindmm[1] # ii != jj

               if jj == common
                  Hsoc[mm,nn] = epsilonsoc(ii,jj,1.0)*sqrt(Np)
               elseif ii == common
                  Hsoc[mm,nn] = epsilonsoc(ii,jj,1.0)*sqrt(Np)
               else
                  Hsoc[mm,nn] = epsilonsoc(ii,jj,1.0)
               end

            end

        end

    end

    # define the Hamiltonian for up up
    Hho[end-(maxmatpcut-1):end,end-(maxmatpcut-1):end] = Hho[1:maxmatpcut,1:maxmatpcut]

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

        for mm = 1:nn

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

            #

        end

    end

    # since the Hamiltonian is helmitian
    Hsoc .= Hsoc + Hsoc' - spdiagm(diag(Hsoc))
    HW .= HW + HW' - spdiagm(diag(HW))

end
