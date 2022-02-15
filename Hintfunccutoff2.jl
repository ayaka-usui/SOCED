using Arpack, SparseArrays, LinearAlgebra

# define functions used here
include("vijkl.jl")

function coefficientInt(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},vecmbindnn2::Vector{Int64},vecmbindmm2::Vector{Int64}Np::Int64)

    # vecmbindnn2 = zeros(Int64,2)
    # vecmbindmm2 = zeros(Int64,2)

    for kk = 1:Np
        for ll = 1:Np
            if vecmbindnn[kk] == vecmbindmm[ll]

               vecmbindnn2 .= vecmbindnn[filter(x->x!=kk,1:Np)]
               vecmbindmm2 .= vecmbindmm[filter(x->x!=ll,1:Np)]

               


            end
        end
    end


    # indices for Vijkl
    ind1 = [vecmbindnn[1]-1, vecmbindnn[2]-1, vecmbindmm[1]-1, vecmbindmm[2]-1]
    sort!(ind1,rev=true)
    ind2 = binomial(ind1[1]+3,4) + binomial(ind1[2]+2,3) + binomial(ind1[3]+1,2) + binomial(ind1[4],1) + 1


    ind0 = 0
    for kk = 1:Np-1
        for ll = kk+1:Np
            if vecmbindnn[kk] == vecmbindnn[ll]
               ind0 += 1
            end
        end
    end






    if vecmbindnn[1] == vecmbindnn[2]
       if vecmbindmm[1] == vecmbindmm[2]
          Hintdown[mm,nn] = vecV[ind2]*factorial(Np)
       else
          Hintdown[mm,nn] = vecV[ind2]*sqrt(factorial(Np))*2
       end
    else
       if vecmbindmm[1] == vecmbindmm[2]
          Hintdown[mm,nn] = vecV[ind2]*sqrt(factorial(Np))*2
       else
          Hintdown[mm,nn] = vecV[ind2]*4
       end
    end

end

function Hintfunccutoff2!(indvec::Vector{Int64}, indvec2::Vector{Int64}, indvec3::Vector{Int64}, Msize0::Int64, Np::Int64, matp::Matrix{Int64}, matp20::Matrix{Int64}, matp21::Matrix{Int64}, matp30::Matrix{Int64}, matp31::Matrix{Int64}, Hintdown::ST, Hintup::ST, Hintdu::ST) where ST <: Union{SparseMatrixCSC{Float64},Array{Float64}}

    maxmatpcut = length(indvec)
    maxmatpcut2 = length(indvec2)
    # maxmatpcut3 = length(indvec3) # maxmatpcut3 == maxmatpcut2

    # register data of Vijkl in a vector
    ind0 = 0
    vecV = zeros(Float64,binomial(Msize0+3,4))
    for n1 = 0:Msize0-1
        for n2 = 0:n1
            for n3 = 0:n2
                for n4 = 0:n3
                    ind0 += 1
                    vecV[ind0] = vijkl(n1,n2,n3,n4)
                end
            end
        end
    end

    # defines vectors and matrices
    Hintdown .= 0. #spzeros(Float64,maxmatpcut,maxmatpcut);
    Hintup .= 0. #spzeros(Float64,maxmatpcut,maxmatpcut);
    Hintdu .= 0. #spzeros(Float64,maxmatpcut,maxmatpcut);
    vecmbindnn = zeros(Int64,Np)
    vecmbindmm = zeros(Int64,Np)

    # define a matrix for the Hamiltonian for down down down
    for nn = 1:maxmatpcut # parfor

        # ket
        in2bind!(indvec[nn],Msize0,Np,matp,vecmbindnn)

        for mm = 1:nn

            # bra
            in2bind!(indvec[mm],Msize0,Np,matp,vecmbindmm)

            ind,coeff = coefficientInt(vecmbindnn,vecmbindmm,Np)
            Hintdown[mm,nn] = vecV[ind]*coeff





            if vecmbindnn[1] == vecmbindnn[2]
               if vecmbindmm[1] == vecmbindmm[2]
                  Hintdown[mm,nn] = vecV[ind2]*factorial(Np)
               else
                  Hintdown[mm,nn] = vecV[ind2]*sqrt(factorial(Np))*2
               end
            else
               if vecmbindmm[1] == vecmbindmm[2]
                  Hintdown[mm,nn] = vecV[ind2]*sqrt(factorial(Np))*2
               else
                  Hintdown[mm,nn] = vecV[ind2]*4
               end
            end

        end

    end

    # use conjectures for the lower triangle elements of Hint since it is hermite
    Hintdown .= Hintdown + Hintdown' - spdiagm(diag(Hintdown))

    # define Hintup
    Hintup[end-(maxmatpcut-1):end,end-(maxmatpcut-1):end] = Hintdown[1:maxmatpcut,1:maxmatpcut]

    # define Hintdu
    vecmbindnn2 = zeros(Int64,Np)

    for nn = 1:maxmatpcut2 # parfor

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

            ind1 = [vecmbindnn[1]-1, vecmbindnn[2]-1, vecmbindmm[1]-1, vecmbindmm[2]-1]
            sort!(ind1,rev=true)
            ind2 = binomial(ind1[1]+3,4) + binomial(ind1[2]+2,3) + binomial(ind1[3]+1,2) + binomial(ind1[4],1) + 1

            Hintdu[maxmatpcut+mm,maxmatpcut+nn] = vecV[ind2]*2

        end

    end

    # use conjectures for the lower triangle elements of Hint since it is hermite
    Hintdu .= Hintdu + Hintdu' - spdiagm(diag(Hintdu))

end
