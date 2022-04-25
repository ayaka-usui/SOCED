using Arpack, SparseArrays, LinearAlgebra

# define functions used here
include("vijkl.jl")

function coefficientpair(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},vecmbindnn3::Vector{Int64},vecmbindmm3::Vector{Int64},vecindcoeff::Matrix{Float64},Np::Int64)

    # vecmbindnn3 = zeros(Int64,2)
    # vecmbindmm3 = zeros(Int64,2)

    # a|n>= sqrt(n)|n-1>
    # a^+|n>= sqrt(n+1)|n+1>

    vecindcoeff .= 0.0 #zeros(Float64,3,5)
    ind0 = 0
    vecmbind0 = 0
    # common = 0

    vecmbindnn3 .= 0
    vecmbindmm3 .= 0

    for pp = 1:Np

        # consider to edit when parfor is inplemented
        if vecmbindnn[pp] == vecmbind0
           continue
        end

        for qq = 1:Np

            if vecmbindnn[pp] == vecmbindmm[qq]

               vecmbindnn3 .= vecmbindnn[findall(x->x!=pp,1:Np)] # ll kk
               vecmbindmm3 .= vecmbindmm[findall(x->x!=qq,1:Np)] # jj ii
               common = vecmbindnn[pp]

               element = 1.0

               # coefficients of operators for ket
               if vecmbindnn3[1] == vecmbindnn3[2] # a_{kk} a_{kk}
                  if common == vecmbindnn3[1]
                     element = sqrt(Np*(Np-1))*element
                  else
                     element = sqrt(2*1)*element
                  end
               else # vecmbindnn3[1] != vecmbindnn3[2] # a_{kk} a_{ll}
                  if common == vecmbindnn3[1] || common == vecmbindnn3[2]
                     element = sqrt(2)*element
                  # else
                     # element = 1.0*element
                  end
                  element = element*2 # a_{kk} a_{ll} + a_{ll} a_{kk}
               end

               # coefficients of operators for bra
               if vecmbindmm3[1] == vecmbindmm3[2] # a^+_{ii} a^+_{ii}
                  if common == vecmbindmm3[1]
                     element = sqrt(2*3)*element
                  else
                     element = sqrt(1*2)*element
                  end
               else # vecmbindmm3[1] != vecmbindmm3[2] # a^+_{ii} a^+_{jj}
                  if common == vecmbindmm3[1] || common == vecmbindmm3[2]
                     element = sqrt(2)*element
                  # else
                     # element = 1.0*element
                  end
                  element = element*2 # a^+_{ii} a^+_{jj} + a^+_{jj} a^+_{ii}
               end

               # indices for Vijkl
               # ind1 = [vecmbindnn3[1]-1, vecmbindnn3[2]-1, vecmbindmm3[1]-1, vecmbindmm3[2]-1] # ii jj kk ll
               # sort!(ind1,rev=true)
               # ind2 = binomial(ind1[1]+3,4) + binomial(ind1[2]+2,3) + binomial(ind1[3]+1,2) + binomial(ind1[4],1) + 1

               ind0 += 1
               vecindcoeff[ind0,1:4] = vecmbindnn3 #ind2
               vecindcoeff[ind0,end] = element

               vecmbind0 = vecmbindnn[pp]

               break

            end

        end
    end

    return vecindcoeff, ind0

end

function coefficientpair2(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},vecmbindnn3::Vector{Int64},vecmbindmm3::Vector{Int64},vecindcoeff::Matrix{Float64},Np::Int64)

    # down down up for g12

    # vecmbindnn3 = zeros(Int64,2)
    # vecmbindmm3 = zeros(Int64,2)

    # a|n>= sqrt(n)|n-1>
    # a^+|n>= sqrt(n+1)|n+1>

    vecindcoeff .= 0.0 #zeros(Float64,3,2)
    ind0 = 0
    vecmbind0 = 0
    # common = 0

    vecmbindnn3 .= 0
    vecmbindmm3 .= 0

    for pp = 1:Np-1

        # consider to edit when parfor is implemented
        if vecmbindnn[pp] == vecmbind0
           continue
        end

        for qq = 1:Np-1

            if vecmbindnn[pp] == vecmbindmm[qq]

               vecmbindnn3 .= vecmbindnn[findall(x->x!=pp,1:Np)] # ll kk
               vecmbindmm3 .= vecmbindmm[findall(x->x!=qq,1:Np)] # jj ii
               common = vecmbindnn[pp]

               element = 1.0

               # coefficients of operators for ket
               if common == vecmbindnn3[1]
                  element = sqrt(2)*element
               # else
               #    element = 1.0*element
               end

               # coefficients of operators for bra
               if common == vecmbindmm3[1]
                  element = sqrt(2)*element
               # else
               #    element = 1.0*element
               end

               # a^+_{down}a^+_{up}a_{down}a_{up} + a^+_{up}a^+_{down}a_{up}a_{down}
               element = 2*element

               # indices for Vijkl
               # ind1 = [vecmbindnn3[1]-1, vecmbindnn3[2]-1, vecmbindmm3[1]-1, vecmbindmm3[2]-1] # ii jj kk ll
               # sort!(ind1,rev=true)
               # ind2 = binomial(ind1[1]+3,4) + binomial(ind1[2]+2,3) + binomial(ind1[3]+1,2) + binomial(ind1[4],1) + 1

               ind0 += 1
               vecindcoeff[ind0,1:4] = vecmbindnn3
               vecindcoeff[ind0,end] = element

               # vecindcoeff[ind0,1] = ind2
               # vecindcoeff[ind0,2] = element

               vecmbind0 = vecmbindnn[pp]

               break

            end

        end
    end

    return vecindcoeff, ind0

end

function coefficientpair3(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},vecmbindnn3::Vector{Int64},vecmbindmm3::Vector{Int64},vecindcoeff::Matrix{Float64},Np::Int64)

    # down down up for g

    # vecmbindnn3 = zeros(Int64,2)
    # vecmbindmm3 = zeros(Int64,2)

    # a|n>= sqrt(n)|n-1>
    # a^+|n>= sqrt(n+1)|n+1>

    vecindcoeff .= 0.0 #zeros(Float64,3,2)
    ind0 = 0
    # vecmbind0 = 0
    # common = 0

    if vecmbindnn[end] != vecmbindmm[end]
       return vecindcoeff, ind0
    end

    vecmbindnn3 .= vecmbindnn[1:2] # ll kk
    vecmbindmm3 .= vecmbindmm[1:2] # jj ii

    element = 1.0

    if vecmbindnn3[1] == vecmbindnn3[2]
       element = sqrt(2)*element
    else # vecmbindnn3[1] != vecmbindnn3[2]
       element = 2*element
    end
    if vecmbindmm3[1] == vecmbindmm3[2]
       element = sqrt(2)*element
    else # vecmbindmm3[1] != vecmbindmm3[2]
       element = 2*element
    end

    # indices for Vijkl
    # ind1 = [vecmbindnn3[1]-1, vecmbindnn3[2]-1, vecmbindmm3[1]-1, vecmbindmm3[2]-1] # ii jj kk ll
    # sort!(ind1,rev=true)
    # ind2 = binomial(ind1[1]+3,4) + binomial(ind1[2]+2,3) + binomial(ind1[3]+1,2) + binomial(ind1[4],1) + 1

    ind0 += 1
    vecindcoeff[ind0,1:4] = vecmbindnn3
    vecindcoeff[ind0,end] = element

    # vecindcoeff[ind0,1] = ind2
    # vecindcoeff[ind0,2] = element

    return vecindcoeff, ind0

end

function paircorrelation_fun!(indvec::Vector{Int64}, indvec2::Vector{Int64}, Msize0::Int64, Np::Int64, matp::Matrix{Int64}, matp20::Matrix{Int64}, matp21::Matrix{Int64}, Mpairdown::ST, Mpairup::ST, Mpairdu::ST, psi::Vector{ComplexF64}) where ST <: Union{SparseMatrixCSC{Float64},Array{Float64}}

    maxmatpcut = length(indvec)
    maxmatpcut2 = length(indvec2)

    # defines vectors and matrices
    Mpairdown .= 0. #spzeros(Float64,maxmatpcut,maxmatpcut,15);
    Mpairup .= 0. #spzeros(Float64,maxmatpcut,maxmatpcut,15);
    Mpairdu .= 0. #spzeros(Float64,maxmatpcut2,maxmatpcut2,15);
    vecmbindnn = zeros(Int64,Np)
    vecmbindmm = zeros(Int64,Np)
    vecmbindnn3 = zeros(Int64,2)
    vecmbindmm3 = zeros(Int64,2)
    vecindcoeff = zeros(Float64,3,5)

    # define a pair correlation for down down down
    for nn = 1:maxmatpcut # parfor

        # ket
        in2bind!(indvec[nn],Msize0,Np,matp,vecmbindnn)

        for mm = 1:nn

            # bra
            in2bind!(indvec[mm],Msize0,Np,matp,vecmbindmm)

            vecindcoeff, ind0 = coefficientpair(vecmbindnn,vecmbindmm,vecmbindnn3,vecmbindmm3,vecindcoeff,Np)
            Mpairdown[mm,nn,1:5] = vecindcoeff[1,1:5]
            if ind0 >= 2
               Mpairdown[mm,nn,5+1:5+5] = vecindcoeff[2,1:5]
            end
            if ind0 == 3
               Mpairdown[mm,nn,10+1:10+5] = vecindcoeff[3,1:5]
            end

            # ind2 = Int64.(vecindcoeff[1:ind0,1])
            # Hintdown[mm,nn] = sum(vecV[ind2].*vecindcoeff[1:ind0,2])

        end

    end

    # define pair correlation for down down up
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

        # a = findfirst(x->x==indvec2[nn],indvec3)

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

            # b = findfirst(x->x==indvec2[mm],indvec3)

            vecindcoeff, ind0 = coefficientpair2(vecmbindnn,vecmbindmm,vecmbindnn3,vecmbindmm3,vecindcoeff,Np)
            Mpairdu[maxmatpcut+mm,maxmatpcut+nn,1:5] = vecindcoeff[1,1:5]
            if ind0 >= 2
               Mpairdu[maxmatpcut+mm,maxmatpcut+nn,5+1:5+5] = vecindcoeff[2,1:5]
            end
            if ind0 == 3
               Mpairdu[maxmatpcut+mm,maxmatpcut+nn,10+1:10+5] = vecindcoeff[3,1:5]
            end
            # ind2 = Int64.(vecindcoeff[1:ind0,1])
            # Hintdu[maxmatpcut+mm,maxmatpcut+nn] = sum(vecV[ind2].*vecindcoeff[1:ind0,2])

            vecindcoeff, ind0 = coefficientpair3(vecmbindnn,vecmbindmm,vecmbindnn3,vecmbindmm3,vecindcoeff,Np)
            Mpairdown[maxmatpcut+mm,maxmatpcut+nn,1:5] = vecindcoeff[1,1:5]
            if ind0 >= 2
               Mpairdown[maxmatpcut+mm,maxmatpcut+nn,5+1:5+5] = vecindcoeff[2,1:5]
            end
            if ind0 == 3
               Mpairdown[maxmatpcut+mm,maxmatpcut+nn,10+1:10+5] = vecindcoeff[3,1:5]
            end
            # ind2 = Int64.(vecindcoeff[1:ind0,1])
            # Hintdown[maxmatpcut+mm,maxmatpcut+nn] = sum(vecV[ind2].*vecindcoeff[1:ind0,2])

        end

    end

    # pair correlation
    fun_nu .= 0.0 #zeros(ComplexF64,Nx,Ny)

    for jjy = 1:Ny
        for jjx = 1:Nx

            for kk1 = 1:maxmatpcut
                for kk0 = 1:maxmatpcut
                    fun_nu[jjx,jjy] += psi[kk1]*(Mpairdown[kk1,kk0,5]*psi[kk0]*funHO(jjx,jjy,Int64.(Mpairdown[kk1,kk0,1:4])) +
                                        Mpairdown[kk1,kk0,10]*psi[kk0]*funHO(jjx,jjy,Int64.(Mpairdown[kk1,kk0,5+1:5+4])) +
                                        Mpairdown[kk1,kk0,15]*psi[kk0]*funHO(jjx,jjy,Int64.(Mpairdown[kk1,kk0,10+1:10+4])))
                end
            end

            for kk1 = maxmatpcut+1:maxmatpcut+maxmatpcut2
                for kk0 = maxmatpcut+1:maxmatpcut+maxmatpcut2
                    fun_nu[jjx,jjy] += psi[kk1]*(Mpairdown[kk1,kk0,5]*psi[kk0]*funHO(jjx,jjy,Int64.(Mpairdown[kk1,kk0,1:4])) +
                                        Mpairdown[kk1,kk0,10]*psi[kk0]*funHO(jjx,jjy,Int64.(Mpairdown[kk1,kk0,5+1:5+4])) +
                                        Mpairdown[kk1,kk0,15]*psi[kk0]*funHO(jjx,jjy,Int64.(Mpairdown[kk1,kk0,10+1:10+4])))
                end
            end

        end
    end






    # use conjectures for the lower triangle elements of Hint since it is hermite
    Hintdown .= Hintdown + Hintdown' - spdiagm(diag(Hintdown))
    Hintdu .= Hintdu + Hintdu' - spdiagm(diag(Hintdu))

    # define Hintup
    Hintup[end-(maxmatpcut-1):end,end-(maxmatpcut-1):end] = Hintdown[1:maxmatpcut,1:maxmatpcut]
    Hintup[end-maxmatpcut-maxmatpcut2+1:end-maxmatpcut,end-maxmatpcut-maxmatpcut2+1:end-maxmatpcut] = Hintdown[1+maxmatpcut:maxmatpcut2+maxmatpcut,1+maxmatpcut:maxmatpcut2+maxmatpcut]

    # define Hint for up up down
    Hintdu[maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2,maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2] = Hintdu[maxmatpcut+1:maxmatpcut+maxmatpcut2,maxmatpcut+1:maxmatpcut+maxmatpcut2]

end
