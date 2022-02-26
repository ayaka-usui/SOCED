using Arpack, SparseArrays, LinearAlgebra

# define functions used here
include("vijkl.jl")

function coefficientInt(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},vecmbindnn3::Vector{Int64},vecmbindmm3::Vector{Int64},common::Vector{Int64},vecindcoeff::Matrix{Float64},Np::Int64)

    # vecmbindnn3 = zeros(Int64,2)
    # vecmbindmm3 = zeros(Int64,2)

    # a|n>= sqrt(n)|n-1>
    # a^+|n>= sqrt(n+1)|n+1>

    vecindcoeff .= 0.0 #zeros(Float64,3,2)
    ind0 = 0
    vecmbind0 = 0

    for pp = 1:Np

        # consider to edit when parfor is inplemented
        if vecmbindnn[pp] == vecmbind0
           continue
        end

        for qq = 1:Np

            if vecmbindnn[pp] == vecmbindmm[qq]

               vecmbindnn3 .= vecmbindnn[findall(x->x!=pp,1:Np)] # ll kk
               vecmbindmm3 .= vecmbindmm[findall(x->x!=qq,1:Np)] # jj ii
               common .= vecmbindnn[pp]

               element = 1.0

               # coefficients of operators for ket
               if vecmbindnn3[1] == vecmbindnn3[2] # a_{kk} a_{kk}
                  if common[1] == vecmbindnn3[1]
                     element = sqrt(Np*(Np-1))*element
                  else
                     element = sqrt(2*1)*element
                  end
               else # vecmbindnn3[1] != vecmbindnn3[2]
                     element = element*2 # a_{kk} a_{ll} + a_{ll} a_{kk}
                     if common[1] == vecmbindnn3[1]
                        element = sqrt(2)*element
                     elseif common[1] == vecmbindnn3[2]
                        element = sqrt(2)*element
                     else
                        element = 1.0*element
                     end
               end

               # coefficients of operators for bra
               if vecmbindmm3[1] == vecmbindmm3[2] # a_{ii} a_{ii}
                  if common[1] == vecmbindmm3[1]
                     element = sqrt(2*3)*element
                  else
                     element = sqrt(1*2)*element
                  end
               else # vecmbindmm3[1] != vecmbindmm3[2]
                  element = element*2 # a_{kk} a_{ll} + a_{ll} a_{kk}
                  if common[1] == vecmbindmm3[1]
                     element = sqrt(2)*element
                  elseif common[1] == vecmbindnn3[2]
                     element = sqrt(2)*element
                  else
                     element = 1.0*element
                  end
               end

               # indices for Vijkl
               ind1 = [vecmbindnn3[1]-1, vecmbindnn3[2]-1, vecmbindmm3[1]-1, vecmbindmm3[2]-1] # ii jj kk ll
               sort!(ind1,rev=true)
               ind2 = binomial(ind1[1]+3,4) + binomial(ind1[2]+2,3) + binomial(ind1[3]+1,2) + binomial(ind1[4],1) + 1

               ind0 += 1
               vecindcoeff[ind0,1] = ind2
               vecindcoeff[ind0,2] = element

               vecmbind0 = vecmbindnn[pp]

               break

            end

        end
    end

    return vecindcoeff, ind0

end

function coefficientInt2(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},vecmbindnn3::Vector{Int64},vecmbindmm3::Vector{Int64},common::Vector{Int64},vecindcoeff::Matrix{Float64},Np::Int64)

    # down down up for g12

    # vecmbindnn3 = zeros(Int64,2)
    # vecmbindmm3 = zeros(Int64,2)

    # a|n>= sqrt(n)|n-1>
    # a^+|n>= sqrt(n+1)|n+1>

    vecindcoeff .= 0.0 #zeros(Float64,3,2)
    ind0 = 0
    vecmbind0 = 0

    for pp = 1:Np-1

        # consider to edit when parfor is inplemented
        if vecmbindnn[pp] == vecmbind0
           continue
        end

        for qq = 1:Np-1

            if vecmbindnn[pp] == vecmbindmm[qq]

               vecmbindnn3 .= vecmbindnn[findall(x->x!=pp,1:Np)] # ll kk
               vecmbindmm3 .= vecmbindmm[findall(x->x!=qq,1:Np)] # jj ii
               common .= vecmbindnn[pp]

               element = 1.0

               # coefficients of operators for ket
               element = element*2 # a_{kk} a_{ll} + a_{ll} a_{kk}
               if common[1] == vecmbindnn3[1]
                  element = sqrt(2)*element
               else
                  element = 1.0*element
               end

               # coefficients of operators for bra
               element = element*2 # a_{kk} a_{ll} + a_{ll} a_{kk}
               if common[1] == vecmbindmm3[1]
                  element = sqrt(2)*element
               else
                  element = 1.0*element
               end

               # indices for Vijkl
               ind1 = [vecmbindnn3[1]-1, vecmbindnn3[2]-1, vecmbindmm3[1]-1, vecmbindmm3[2]-1] # ii jj kk ll
               sort!(ind1,rev=true)
               ind2 = binomial(ind1[1]+3,4) + binomial(ind1[2]+2,3) + binomial(ind1[3]+1,2) + binomial(ind1[4],1) + 1

               ind0 += 1
               vecindcoeff[ind0,1] = ind2
               vecindcoeff[ind0,2] = element

               vecmbind0 = vecmbindnn[pp]

               break

            end

        end
    end

    return vecindcoeff, ind0

end

function coefficientInt3(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},vecmbindnn3::Vector{Int64},vecmbindmm3::Vector{Int64},common::Vector{Int64},vecindcoeff::Matrix{Float64},Np::Int64)

    # down down up for g

    # vecmbindnn3 = zeros(Int64,2)
    # vecmbindmm3 = zeros(Int64,2)

    # a|n>= sqrt(n)|n-1>
    # a^+|n>= sqrt(n+1)|n+1>

    vecindcoeff .= 0.0 #zeros(Float64,3,2)
    ind0 = 0
    vecmbind0 = 0

    if vecmbindnn[end] != vecmbindmm[end]
       return vecindcoeff, ind0
    end

    vecmbindnn3 .= vecmbindnn[1:2] # ll kk
    vecmbindmm3 .= vecmbindmm[1:2] # jj ii



    



    for pp = 1:Np-1

        # consider to edit when parfor is inplemented
        if vecmbindnn[pp] == vecmbind0
           continue
        end

        for qq = 1:Np-1

            if vecmbindnn[pp] == vecmbindmm[qq]

               vecmbindnn3 .= vecmbindnn[findall(x->x!=pp,1:Np)] # ll kk
               vecmbindmm3 .= vecmbindmm[findall(x->x!=qq,1:Np)] # jj ii
               common .= vecmbindnn[pp]

               element = 1.0

               # coefficients of operators for ket
               element = element*2 # a_{kk} a_{ll} + a_{ll} a_{kk}
               if common[1] == vecmbindnn3[1]
                  element = sqrt(2)*element
               else
                  element = 1.0*element
               end

               # coefficients of operators for bra
               element = element*2 # a_{kk} a_{ll} + a_{ll} a_{kk}
               if common[1] == vecmbindmm3[1]
                  element = sqrt(2)*element
               else
                  element = 1.0*element
               end

               # indices for Vijkl
               ind1 = [vecmbindnn3[1]-1, vecmbindnn3[2]-1, vecmbindmm3[1]-1, vecmbindmm3[2]-1] # ii jj kk ll
               sort!(ind1,rev=true)
               ind2 = binomial(ind1[1]+3,4) + binomial(ind1[2]+2,3) + binomial(ind1[3]+1,2) + binomial(ind1[4],1) + 1

               ind0 += 1
               vecindcoeff[ind0,1] = ind2
               vecindcoeff[ind0,2] = element

               vecmbind0 = vecmbindnn[pp]

               break

            end

        end
    end

    return vecindcoeff, ind0

end

function Hintfunccutoff2!(indvec::Vector{Int64}, indvec2::Vector{Int64}, Msize0::Int64, Np::Int64, matp::Matrix{Int64}, matp20::Matrix{Int64}, matp21::Matrix{Int64}, Hintdown::ST, Hintup::ST, Hintdu::ST) where ST <: Union{SparseMatrixCSC{Float64},Array{Float64}}

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
    vecmbindnn3 = zeros(Int64,2)
    vecmbindmm3 = zeros(Int64,2)
    vecindcoeff = zeros(Float64,3,2)
    common = zeros(Int64,Np-1)

    # define a matrix for the Hamiltonian for down down down
    for nn = 1:maxmatpcut # parfor

        # ket
        in2bind!(indvec[nn],Msize0,Np,matp,vecmbindnn)

        for mm = 1:nn

            # bra
            in2bind!(indvec[mm],Msize0,Np,matp,vecmbindmm)

            vecindcoeff, ind0 = coefficientInt(vecmbindnn,vecmbindmm,vecmbindnn3,vecmbindmm3,common[1:Np-2],vecindcoeff,Np)
            ind2 = Int64.(vecindcoeff[1:ind0,1])
            Hintdown[mm,nn] = sum(vecV[ind2].*vecindcoeff[1:ind0,2])

        end

    end

    # define Hint for down down up
    maxmatp21 = matp21[Msize0+1,1+1]
    vecmbindnn2 = zeros(Int64,Np)
    vecmbindmm2 = zeros(Int64,Np)

    for nn = 1:maxmatpcut2 # parfor # down down up for ket

        # ket
        # indvec2[nn] = (nndown-1)*maxmatp21 + nnup
        indketup = mod(indvec2[nn],maxmatp21)
        if indketup != 0
           indketdown = div(indvec2[nn],maxmatp21) + 1
        else # indketup == 0
           indketdown = div(indvec2[nn],maxmatp21)
           indketup = maxmatp21
        end
        in2bind!(indketdown,Msize0,Np-1,matp20,vecmbindnn2)
        vecmbindnn[1:Np-1] = vecmbindnn2[1:Np-1]
        in2bind!(indketup,Msize0,1,matp21,vecmbindnn2)
        vecmbindnn[Np:Np] = vecmbindnn2[1:1]

        for mm = 1:nn # down down up for bra

            # bra
            # indvec2[mm] = (mmdown-1)*maxmatp21 + mmup
            indbraup = mod(indvec2[mm],maxmatp21)
            if indbraup != 0
               indbradown = div(indvec2[mm],maxmatp21) + 1
            else # indketup == 0
               indbradown = div(indvec2[mm],maxmatp21)
               indbraup = maxmatp21
            end
            in2bind!(indbradown,Msize0,Np-1,matp20,vecmbindmm2)
            vecmbindmm[1:Np-1] = vecmbindmm2[1:Np-1]
            in2bind!(indbraup,Msize0,1,matp21,vecmbindmm2)
            vecmbindmm[Np:Np] = vecmbindmm2[1:1]

            vecindcoeff, ind0 = coefficientInt2(vecmbindnn,vecmbindmm,vecmbindnn3,vecmbindmm3,common[1:Np-2],vecindcoeff,Np)
            ind2 = Int64.(vecindcoeff[1:ind0,1])
            Hintdu[maxmatpcut+mm,maxmatpcut+nn] = sum(vecV[ind2].*vecindcoeff[1:ind0,2])

            vecindcoeff, ind0 = coefficientInt3(vecmbindnn,vecmbindmm,vecmbindnn3,vecmbindmm3,common[1:Np-2],vecindcoeff,Np)
            ind2 = Int64.(vecindcoeff[1:ind0,1])
            Hintdown[maxmatpcut+mm,maxmatpcut+nn] = sum(vecV[ind2].*vecindcoeff[1:ind0,2])

        end

    end

    # define Hint for down up up
    # indNp = 2
    # maxmatp31 = matp31[Msize0+1,2+1]
    #
    # for nn = 1:maxmatpcut2 # parfor
    #
    #     # ket
    #     # indvec3[nn] = (nndown-1)*maxmatp31 + nnup
    #     indketup = mod(indvec3[nn],maxmatp31)
    #     if indketup != 0
    #        indketdown = div(indvec3[nn],maxmatp31) + 1
    #     else # indketup == 0
    #        indketdown = div(indvec3[nn],maxmatp31)
    #        indketup = maxmatp31
    #     end
    #     in2bind!(indketdown,Msize0,Np-indNp,matp30,vecmbindnn2)
    #     vecmbindnn[1:Np-indNp] = vecmbindnn2[1:Np-indNp]
    #     in2bind!(indketup,Msize0,indNp,matp31,vecmbindnn2)
    #     vecmbindnn[Np-indNp+1:Np] = vecmbindnn2[1:indNp]
    #
    #     for mm = 1:nn
    #
    #         # bra
    #         # indvec2[mm] = (mmdown-1)*maxmatp31 + mmup
    #         indbraup = mod(indvec3[mm],maxmatp31)
    #         if indbraup != 0
    #            indbradown = div(indvec3[mm],maxmatp31) + 1
    #         else # indketup == 0
    #            indbradown = div(indvec3[mm],maxmatp31)
    #            indbraup = maxmatp31
    #         end
    #         in2bind!(indbradown,Msize0,Np-indNp,matp30,vecmbindmm2)
    #         vecmbindmm[1:Np-indNp] = vecmbindmm2[1:Np-indNp]
    #         in2bind!(indbraup,Msize0,indNp,matp31,vecmbindmm2)
    #         vecmbindmm[Np-indNp+1:Np] = vecmbindmm2[1:indNp]
    #
    #         vecindcoeff, ind0 = coefficientInt3(vecmbindnn,vecmbindmm,vecmbindnn3,vecmbindmm3,common[1:Np-2],vecindcoeff,Np)
    #         ind2 = Int64.(vecindcoeff[1:ind0,1])
    #         Hintdu[maxmatpcut+maxmatpcut2+mm,maxmatpcut+maxmatpcut2+nn] = sum(vecV[ind2].*vecindcoeff[1:ind0,2])
    #
    #     end
    #
    # end

    # use conjectures for the lower triangle elements of Hint since it is hermite
    Hintdown .= Hintdown + Hintdown' - spdiagm(diag(Hintdown))
    Hintdu .= Hintdu + Hintdu' - spdiagm(diag(Hintdu))

    # define Hintup
    Hintup[end-(maxmatpcut-1):end,end-(maxmatpcut-1):end] = Hintdown[1:maxmatpcut,1:maxmatpcut]

    # define Hint for up up down
    Hintdu[maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2,maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2] = Hintdu[maxmatpcut+1:maxmatpcut+maxmatpcut2,maxmatpcut+1:maxmatpcut+maxmatpcut2]

end
