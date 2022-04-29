
using Arpack, SparseArrays, LinearAlgebra
using Polynomials, SpecialPolynomials

# define functions used here
# include("vijkl.jl")

function setfunHO(xrange::LinRange{Float64},Msize0::Int64)

    # y = variable(Polynomial{Rational{Int}})
    Nx = length(xrange)
    outcome = zeros(Float64,Nx,Msize0)
    dx = abs(xrange[2]-xrange[1])

    for nn = 1:Msize0
        for jjx = 1:Nx
            outcome[jjx,nn] = basis(Hermite,nn-1)(xrange[jjx])*exp(-xrange[jjx]^2/2)
        end
        outcome[:,nn] = outcome[:,nn]/sqrt(sum(abs.(outcome[:,nn]).^2)*dx)
    end

    return outcome

end

function coefficientpair(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},vecmbindnn3::Vector{Int64},vecmbindmm3::Vector{Int64},vecindcoeff0::Matrix{Int64},vecindcoeff1::Vector{Float64},Np::Int64)

    # vecmbindnn3 = zeros(Int64,2)
    # vecmbindmm3 = zeros(Int64,2)

    # a|n>= sqrt(n)|n-1>
    # a^+|n>= sqrt(n+1)|n+1>

    vecindcoeff0 .= 0 #zeros(Float64,3,4)
    vecindcoeff1 .= 0.0 #zeros(Float64,3)
    ind0 = 0
    vecmbind0 = 0
    # common = 0

    vecmbindnn3 .= 0
    vecmbindmm3 .= 0

    for pp = 1:Np

        # consider to edit when parfor is implemented
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
                  # element = element*2 # a_{kk} a_{ll} + a_{ll} a_{kk}
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
                  # element = element*2 # a^+_{ii} a^+_{jj} + a^+_{jj} a^+_{ii}
               end

               # a^+_{ii} a^+_{jj} a_{kk} a_{ll}
               ind1 = [vecmbindmm3[2], vecmbindmm3[1], vecmbindnn3[2], vecmbindnn3[1]] # ii jj kk ll
               ind0 += 1
               vecindcoeff0[ind0,1:4] = ind1
               vecindcoeff1[ind0] = element

               # a^+_{ii} a^+_{jj} a_{ll} a_{kk}
               if vecmbindnn3[1] != vecmbindnn3[2]
                  ind1 .= [vecmbindmm3[2], vecmbindmm3[1], vecmbindnn3[1], vecmbindnn3[2]]
                  ind0 += 1
                  vecindcoeff0[ind0,1:4] = ind1
                  vecindcoeff1[ind0] = element
               end

               # a^+_{jj} a^+_{ii} a_{kk} a_{ll}
               if vecmbindmm3[1] != vecmbindmm3[2]
                  ind1 .= [vecmbindmm3[1], vecmbindmm3[2], vecmbindnn3[2], vecmbindnn3[1]]
                  ind0 += 1
                  vecindcoeff0[ind0,1:4] = ind1
                  vecindcoeff1[ind0] = element
               end

               # a^+_{jj} a^+_{ii} a_{ll} a_{kk}
               if vecmbindnn3[1] != vecmbindnn3[2] && vecmbindmm3[1] != vecmbindmm3[2]
                  ind1 .= [vecmbindmm3[1], vecmbindmm3[2], vecmbindnn3[1], vecmbindnn3[2]]
                  ind0 += 1
                  vecindcoeff0[ind0,1:4] = ind1
                  vecindcoeff1[ind0] = element
               end

               vecmbind0 = vecmbindnn[pp]

               break

            end

        end
    end

    return vecindcoeff0, vecindcoeff1, ind0

end

function coefficientpair2(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},vecmbindnn3::Vector{Int64},vecmbindmm3::Vector{Int64},vecindcoeff0::Matrix{Int64},vecindcoeff1::Vector{Float64},Np::Int64)

    # down down up for g12

    # vecmbindnn3 = zeros(Int64,2)
    # vecmbindmm3 = zeros(Int64,2)

    # a|n>= sqrt(n)|n-1>
    # a^+|n>= sqrt(n+1)|n+1>

    vecindcoeff0 .= 0 #zeros(Float64,3,2)
    vecindcoeff1 .= 0.0 #zeros(Float64,3)
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
               end

               # coefficients of operators for bra
               if common == vecmbindmm3[1]
                  element = sqrt(2)*element
               end

               # a^+_{down}a^+_{up}a_{down}a_{up} + a^+_{up}a^+_{down}a_{up}a_{down}
               # element = 2*element

               # a^+_{down}a^+_{up}a_{down}a_{up}
               ind1 = [vecmbindmm3[1], vecmbindmm3[2], vecmbindnn3[1], vecmbindnn3[2]] # ii jj kk ll
               ind0 += 1
               vecindcoeff0[ind0,1:4] = ind1
               vecindcoeff1[ind0,end] = element

               # a^+_{up}a^+_{down}a_{up}a_{down}
               ind1 .= [vecmbindmm3[2], vecmbindmm3[1], vecmbindnn3[2], vecmbindnn3[1]] # ii jj kk ll
               ind0 += 1
               vecindcoeff0[ind0,1:4] = ind1
               vecindcoeff1[ind0,end] = element

               vecmbind0 = vecmbindnn[pp]

               break

            end

        end
    end

    return vecindcoeff0, vecindcoeff1, ind0

end

function coefficientpair3(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},vecmbindnn3::Vector{Int64},vecmbindmm3::Vector{Int64},vecindcoeff0::Matrix{Int64},vecindcoeff1::Vector{Float64},Np::Int64)

    # down down up for g

    # vecmbindnn3 = zeros(Int64,2)
    # vecmbindmm3 = zeros(Int64,2)

    # a|n>= sqrt(n)|n-1>
    # a^+|n>= sqrt(n+1)|n+1>

    vecindcoeff0 .= 0 #zeros(Float64,3,2)
    vecindcoeff1 .= 0.0 #zeros(Float64,3,2)
    ind0 = 0
    # vecmbind0 = 0
    # common = 0

    if vecmbindnn[end] != vecmbindmm[end]
       return vecindcoeff0, vecindcoeff1, ind0
    end

    vecmbindnn3 .= vecmbindnn[1:2] # ll kk
    vecmbindmm3 .= vecmbindmm[1:2] # jj ii

    element = 1.0

    # note that the order of the indices matters if vecmbindnn3[1] != vecmbindnn3[2]
    if vecmbindnn3[1] == vecmbindnn3[2]
       element = sqrt(2)*element
    end
    if vecmbindmm3[1] == vecmbindmm3[2]
       element = sqrt(2)*element
    end

    # a^+_{ii} a^+_{jj} a_{kk} a_{ll}
    ind1 = [vecmbindmm3[2], vecmbindmm3[1], vecmbindnn3[2], vecmbindnn3[1]] # ii jj kk ll
    ind0 += 1
    vecindcoeff0[ind0,1:4] = ind1 #vecmbindnn3
    vecindcoeff1[ind0,end] = element

    # a^+_{ii} a^+_{jj} a_{ll} a_{kk}
    if vecmbindnn3[1] != vecmbindnn3[2]
       ind1 .= [vecmbindmm3[2], vecmbindmm3[1], vecmbindnn3[1], vecmbindnn3[2]]
       ind0 += 1
       vecindcoeff0[ind0,1:4] = ind1
       vecindcoeff1[ind0,end] = element
    end

    # a^+_{jj} a^+_{ii} a_{kk} a_{ll}
    if vecmbindmm3[1] != vecmbindmm3[2]
       ind1 .= [vecmbindmm3[1], vecmbindmm3[2], vecmbindnn3[2], vecmbindnn3[1]]
       ind0 += 1
       vecindcoeff0[ind0,1:4] = ind1
       vecindcoeff1[ind0,end] = element
    end

    # a^+_{jj} a^+_{ii} a_{ll} a_{kk}
    if vecmbindnn3[1] != vecmbindnn3[2] && vecmbindmm3[1] != vecmbindmm3[2]
       ind1 .= [vecmbindmm3[1], vecmbindmm3[2], vecmbindnn3[1], vecmbindnn3[2]]
       ind0 += 1
       vecindcoeff0[ind0,1:4] = ind1
       vecindcoeff1[ind0,end] = element
    end

    return vecindcoeff0, vecindcoeff1, ind0

end

function paircorrelation_fun(indvec::Vector{Int64}, indvec2::Vector{Int64}, Msize0::Int64, Np::Int64, matp::Matrix{Int64}, matp20::Matrix{Int64}, matp21::Matrix{Int64}, psi::Vector{ComplexF64}, xrange::LinRange{Float64}, yrange::LinRange{Float64})

    maxmatpcut = length(indvec)
    maxmatpcut2 = length(indvec2)

    # defines vectors and matrices
    # Mpairdown3int .= 0
    # Mpairdown2up1int .= 0
    # Mpairdown1up2int .= 0
    # Mpairup3int .= 0
    # Mpairdown3float .= 0.0
    # Mpairdown2up1float .= 0.0
    # Mpairdown1up2float .= 0.0
    # Mpairup3float .= 0.0
    vecmbindnn = zeros(Int64,Np)
    vecmbindmm = zeros(Int64,Np)
    vecmbindnn3 = zeros(Int64,2)
    vecmbindmm3 = zeros(Int64,2)
    vecindcoeff0 = zeros(Int64,3*4,4) #zeros(Int64,3,4)
    vecindcoeff1 = zeros(Float64,3*4) #zeros(Float64,3)

    # pair correlation
    Nx = length(xrange)
    Ny = length(yrange)
    fun_nudown = zeros(Float64,Nx,Ny)
    fun_nudu = zeros(Float64,Nx,Ny)
    fun_nuup = zeros(Float64,Nx,Ny)
    phiHO = setfunHO(xrange,Msize0)

    fun_nudown1 = zeros(Float64,Nx,Ny)
    fun_nuup1 = zeros(Float64,Nx,Ny)

    # define a pair correlation for down down down
    for nn = 1:maxmatpcut # parfor

        # ket
        in2bind!(indvec[nn],Msize0,Np,matp,vecmbindnn)

        for mm = 1:nn

            # bra
            in2bind!(indvec[mm],Msize0,Np,matp,vecmbindmm)

            vecindcoeff0, vecindcoeff1, ind0 = coefficientpair(vecmbindnn,vecmbindmm,vecmbindnn3,vecmbindmm3,vecindcoeff0,vecindcoeff1,Np)

            if nn == mm

               for jjx = 1:Nx
                   for jjy = 1:Ny
                       for jj = 1:ind0
                           fun_nudown[jjx,jjy] += abs(conj(psi[mm])*psi[nn])*vecindcoeff1[jj]*phiHO[jjx,vecindcoeff0[jj,1]]*
                                                                                              phiHO[jjy,vecindcoeff0[jj,2]]*
                                                                                              phiHO[jjx,vecindcoeff0[jj,3]]*
                                                                                              phiHO[jjy,vecindcoeff0[jj,4]]

                           fun_nuup[jjx,jjy] += abs(conj(psi[maxmatpcut+maxmatpcut2*2+mm])*psi[maxmatpcut+maxmatpcut2*2+nn])*vecindcoeff1[jj]*phiHO[jjx,vecindcoeff0[jj,1]]*
                                                                                                                                              phiHO[jjy,vecindcoeff0[jj,2]]*
                                                                                                                                              phiHO[jjx,vecindcoeff0[jj,3]]*
                                                                                                                                              phiHO[jjy,vecindcoeff0[jj,4]]

                       end
                   end
               end

            else # mm < nn

               for jjx = 1:Nx
                   for jjy = 1:Ny
                       for jj = 1:ind0
                           fun_nudown[jjx,jjy] += 2*real(conj(psi[mm])*psi[nn])*vecindcoeff1[jj]*phiHO[jjx,vecindcoeff0[jj,1]]*
                                                                                                 phiHO[jjy,vecindcoeff0[jj,2]]*
                                                                                                 phiHO[jjx,vecindcoeff0[jj,3]]*
                                                                                                 phiHO[jjy,vecindcoeff0[jj,4]]

                           fun_nuup[jjx,jjy] += 2*real(conj(psi[maxmatpcut+maxmatpcut2*2+mm])*psi[maxmatpcut+maxmatpcut2*2+nn])*vecindcoeff1[jj]*phiHO[jjx,vecindcoeff0[jj,1]]*
                                                                                                                                                 phiHO[jjy,vecindcoeff0[jj,2]]*
                                                                                                                                                 phiHO[jjx,vecindcoeff0[jj,3]]*
                                                                                                                                                 phiHO[jjy,vecindcoeff0[jj,4]]
                       end
                   end
               end

            end

            # for jj = 1:ind0
            #     # Mpairdown3int[mm,nn,4*(jj-1)+1:4*(jj-1)+4] = vecindcoeff0[jj,1:4]
            #     # Mpairdown3float[mm,nn,jj] = vecindcoeff1[jj]
            #     Mpairdown3int[(nn-1)*maxmatpcut+mm,4*(jj-1)+1:4*(jj-1)+4] = vecindcoeff0[jj,1:4]
            #     Mpairdown3float[(nn-1)*maxmatpcut+mm,jj] = vecindcoeff1[jj]
            # end

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

            vecindcoeff0, vecindcoeff1, ind0 = coefficientpair2(vecmbindnn,vecmbindmm,vecmbindnn3,vecmbindmm3,vecindcoeff0,vecindcoeff1,Np)

            if nn == mm
               for jjx = 1:Nx
                   for jjy = 1:Ny
                       for jj = 1:ind0
                           fun_nudu[jjx,jjy] += abs(conj(psi[maxmatpcut+mm])*psi[maxmatpcut+nn])*vecindcoeff1[jj]*phiHO[jjx,vecindcoeff0[jj,1]]*
                                                                                                                  phiHO[jjy,vecindcoeff0[jj,2]]*
                                                                                                                  phiHO[jjx,vecindcoeff0[jj,3]]*
                                                                                                                  phiHO[jjy,vecindcoeff0[jj,4]]

                           fun_nudu[jjx,jjy] += abs(conj(psi[maxmatpcut+maxmatpcut2+mm])*psi[maxmatpcut+maxmatpcut2+nn])*vecindcoeff1[jj]*phiHO[jjx,vecindcoeff0[jj,1]]*
                                                                                                                                          phiHO[jjy,vecindcoeff0[jj,2]]*
                                                                                                                                          phiHO[jjx,vecindcoeff0[jj,3]]*
                                                                                                                                          phiHO[jjy,vecindcoeff0[jj,4]]
                       end
                   end
               end

            else # mm < nn
               for jjx = 1:Nx
                   for jjy = 1:Ny
                       for jj = 1:ind0
                           fun_nudu[jjx,jjy] += 2*real(conj(psi[maxmatpcut+mm])*psi[maxmatpcut+nn])*vecindcoeff1[jj]*phiHO[jjx,vecindcoeff0[jj,1]]*
                                                                                                                     phiHO[jjy,vecindcoeff0[jj,2]]*
                                                                                                                     phiHO[jjx,vecindcoeff0[jj,3]]*
                                                                                                                     phiHO[jjy,vecindcoeff0[jj,4]]

                           fun_nudu[jjx,jjy] += 2*real(conj(psi[maxmatpcut+maxmatpcut2+mm])*psi[maxmatpcut+maxmatpcut2+nn])*vecindcoeff1[jj]*phiHO[jjx,vecindcoeff0[jj,1]]*
                                                                                                                                             phiHO[jjy,vecindcoeff0[jj,2]]*
                                                                                                                                             phiHO[jjx,vecindcoeff0[jj,3]]*
                                                                                                                                             phiHO[jjy,vecindcoeff0[jj,4]]
                       end
                   end
               end
            end

            vecindcoeff0, vecindcoeff1, ind0 = coefficientpair3(vecmbindnn,vecmbindmm,vecmbindnn3,vecmbindmm3,vecindcoeff0,vecindcoeff1,Np)

            if nn == mm
               for jjx = 1:Nx
                   for jjy = 1:Ny
                       for jj = 1:ind0
                           fun_nudown1[jjx,jjy] += abs(conj(psi[maxmatpcut+mm])*psi[maxmatpcut+nn])*vecindcoeff1[jj]*phiHO[jjx,vecindcoeff0[jj,1]]*
                                                                                                                     phiHO[jjy,vecindcoeff0[jj,2]]*
                                                                                                                     phiHO[jjx,vecindcoeff0[jj,3]]*
                                                                                                                     phiHO[jjy,vecindcoeff0[jj,4]]

                           fun_nuup1[jjx,jjy] += abs(conj(psi[maxmatpcut+maxmatpcut2+mm])*psi[maxmatpcut+maxmatpcut2+nn])*vecindcoeff1[jj]*phiHO[jjx,vecindcoeff0[jj,1]]*
                                                                                                                                           phiHO[jjy,vecindcoeff0[jj,2]]*
                                                                                                                                           phiHO[jjx,vecindcoeff0[jj,3]]*
                                                                                                                                           phiHO[jjy,vecindcoeff0[jj,4]]
                       end
                   end
               end

            else # mm < nn
               for jjx = 1:Nx
                   for jjy = 1:Ny
                       for jj = 1:ind0
                           fun_nudown1[jjx,jjy] += 2*real(conj(psi[maxmatpcut+mm])*psi[maxmatpcut+nn])*vecindcoeff1[jj]*phiHO[jjx,vecindcoeff0[jj,1]]*
                                                                                                                        phiHO[jjy,vecindcoeff0[jj,2]]*
                                                                                                                        phiHO[jjx,vecindcoeff0[jj,3]]*
                                                                                                                        phiHO[jjy,vecindcoeff0[jj,4]]

                           fun_nuup1[jjx,jjy] += 2*real(conj(psi[maxmatpcut+maxmatpcut2+mm])*psi[maxmatpcut+maxmatpcut2+nn])*vecindcoeff1[jj]*phiHO[jjx,vecindcoeff0[jj,1]]*
                                                                                                                                              phiHO[jjy,vecindcoeff0[jj,2]]*
                                                                                                                                              phiHO[jjx,vecindcoeff0[jj,3]]*
                                                                                                                                              phiHO[jjy,vecindcoeff0[jj,4]]
                       end
                   end
               end
            end

        end

    end

    #
    # for kk0 = 1:maxmatpcut
    #     for kk1 = kk0+1:maxmatpcut
    #         # Mpairdown3int[kk1,kk0,:] = Mpairdown3int[kk0,kk1,:]
    #         # Mpairdown3float[kk1,kk0,:] = Mpairdown3float[kk0,kk1,:]
    #         Mpairdown3int[(kk0-1)*maxmatpcut+kk1,:] = Mpairdown3int[(kk1-1)*maxmatpcut+kk0,:]
    #         Mpairdown3float[(kk0-1)*maxmatpcut+kk1,:] = Mpairdown3float[(kk1-1)*maxmatpcut+kk0,:]
    #     end
    # end

    # for kk0 = 1:maxmatpcut2 #maxmatpcut+1:maxmatpcut+maxmatpcut2
    #     for kk1 = kk0+1:maxmatpcut2
    #         # Mpairdown3int[maxmatpcut+kk1,maxmatpcut+kk0,:] = Mpairdown3int[maxmatpcut+kk0,maxmatpcut+kk1,:]
    #         # Mpairdown3float[maxmatpcut+kk1,maxmatpcut+kk0,:] = Mpairdown3float[maxmatpcut+kk0,maxmatpcut+kk1,:]
    #         # Mpairdown2up1int[kk1,kk0,:] = Mpairdown2up1int[kk0,kk1,:]
    #         # Mpairdown2up1float[kk1,kk0,:] = Mpairdown2up1float[kk0,kk1,:]
    #         Mpairdown3int[(maxmatpcut+kk0-1)*(maxmatpcut+maxmatpcut2)+maxmatpcut+kk1,:] = Mpairdown3int[(maxmatpcut+kk1-1)*(maxmatpcut+maxmatpcut2)+maxmatpcut+kk0,:]
    #         Mpairdown3float[(maxmatpcut+kk0-1)*(maxmatpcut+maxmatpcut2)+maxmatpcut+kk1,:] = Mpairdown3float[(maxmatpcut+kk1-1)*(maxmatpcut+maxmatpcut2)+maxmatpcut+kk0,:]
    #         Mpairdown2up1int[(kk0-1)*maxmatpcut2+kk1,:] = Mpairdown2up1int[(kk1-1)*maxmatpcut2+kk0,:]
    #         Mpairdown2up1float[(kk0-1)*maxmatpcut2+kk1,:] = Mpairdown2up1float[(kk1-1)*maxmatpcut2+kk0,:]
    #     end
    # end

    # Mpairup3int[maxmatpcut2+1:maxmatpcut2+maxmatpcut,maxmatpcut2+1:maxmatpcut2+maxmatpcut,:] = Mpairdown3int[1:maxmatpcut,1:maxmatpcut,:]
    # Mpairup3float[maxmatpcut2+1:maxmatpcut2+maxmatpcut,maxmatpcut2+1:maxmatpcut2+maxmatpcut,:] = Mpairdown3float[1:maxmatpcut,1:maxmatpcut,:]
    # Mpairup3int[1:maxmatpcut2,1:maxmatpcut2,:] = Mpairdown3int[1+maxmatpcut:maxmatpcut2+maxmatpcut,1+maxmatpcut:maxmatpcut2+maxmatpcut,:]
    # Mpairup3float[1:maxmatpcut2,1:maxmatpcut2,:] = Mpairdown3float[1+maxmatpcut:maxmatpcut2+maxmatpcut,1+maxmatpcut:maxmatpcut2+maxmatpcut,:]

    # Mpairdown1up2int .= Mpairdown2up1int
    # Mpairdown1up2float .= Mpairdown2up1float

    # # pair correlation
    # Nx = length(xrange)
    # Ny = length(yrange)
    #
    # fun_nudown = zeros(ComplexF64,Nx,Ny)
    # fun_nudu = zeros(ComplexF64,Nx,Ny)
    # fun_nuup = zeros(ComplexF64,Nx,Ny)
    #
    # phiHO = setfunHO(xrange,Msize0)

    # for jjy = 1:Ny
    #
    #     println("jjy=",jjy)
    #
    #     for jjx = 1:Nx
    #
    #         for kk1 = 1:maxmatpcut
    #             for kk0 = 1:maxmatpcut
    #                 for jj = 1:3
    #
    #                     # if Mpairdown3int[kk1,kk0,(jj-1)*4+1] != 0
    #                     #    fun_nudown[jjx,jjy] += conj(psi[kk1])*Mpairdown3float[kk1,kk0,jj]*psi[kk0]*phiHO[jjx,Mpairdown3int[kk1,kk0,(jj-1)*4+1]]*
    #                     #                                                                               phiHO[jjx,Mpairdown3int[kk1,kk0,(jj-1)*4+3]]*
    #                     #                                                                               phiHO[jjy,Mpairdown3int[kk1,kk0,(jj-1)*4+2]]*
    #                     #                                                                               phiHO[jjy,Mpairdown3int[kk1,kk0,(jj-1)*4+4]]
    #                     # end
    #                     #
    #                     # if Mpairup3int[maxmatpcut2+kk1,maxmatpcut2+kk0,(jj-1)*4+1] != 0
    #                     #    fun_nuup[jjx,jjy] += conj(psi[maxmatpcut+maxmatpcut2*2+kk1])*Mpairup3float[maxmatpcut2+kk1,maxmatpcut2+kk0,jj]*psi[maxmatpcut+maxmatpcut2*2+kk0]*phiHO[jjx,Mpairup3int[maxmatpcut2+kk1,maxmatpcut2+kk0,(jj-1)*4+1]]*
    #                     #                                                                                                                                                     phiHO[jjx,Mpairup3int[maxmatpcut2+kk1,maxmatpcut2+kk0,(jj-1)*4+3]]*
    #                     #                                                                                                                                                     phiHO[jjy,Mpairup3int[maxmatpcut2+kk1,maxmatpcut2+kk0,(jj-1)*4+2]]*
    #                     #                                                                                                                                                     phiHO[jjy,Mpairup3int[maxmatpcut2+kk1,maxmatpcut2+kk0,(jj-1)*4+4]]
    #                     # end
    #
    #                     indkk0 = kk0
    #                     indkk1 = kk1
    #
    #                     if kk1 >= kk0
    #                        indkk1 = kk0
    #                        indkk0 = kk1
    #                     end
    #
    #                     if Mpairdown3int[(indkk0-1)*maxmatpcut+indkk1,(jj-1)*4+1] != 0
    #                        fun_nudown[jjx,jjy] += conj(psi[kk1])*Mpairdown3float[(indkk0-1)*maxmatpcut+indkk1,jj]*psi[kk0]*phiHO[jjx,Mpairdown3int[(indkk0-1)*maxmatpcut+indkk1,(jj-1)*4+1]]*
    #                                                                                                                        phiHO[jjx,Mpairdown3int[(indkk0-1)*maxmatpcut+indkk1,(jj-1)*4+3]]*
    #                                                                                                                        phiHO[jjy,Mpairdown3int[(indkk0-1)*maxmatpcut+indkk1,(jj-1)*4+2]]*
    #                                                                                                                        phiHO[jjy,Mpairdown3int[(indkk0-1)*maxmatpcut+indkk1,(jj-1)*4+4]]
    #
    #                        fun_nuup[jjx,jjy] += conj(psi[maxmatpcut+maxmatpcut2*2+kk1])*Mpairdown3float[(indkk0-1)*maxmatpcut+indkk1,jj]*psi[maxmatpcut+maxmatpcut2*2+kk0]*phiHO[jjx,Mpairdown3int[(indkk0-1)*maxmatpcut+indkk1,(jj-1)*4+1]]*
    #                                                                                                                                                                        phiHO[jjx,Mpairdown3int[(indkk0-1)*maxmatpcut+indkk1,(jj-1)*4+3]]*
    #                                                                                                                                                                        phiHO[jjy,Mpairdown3int[(indkk0-1)*maxmatpcut+indkk1,(jj-1)*4+2]]*
    #                                                                                                                                                                        phiHO[jjy,Mpairdown3int[(indkk0-1)*maxmatpcut+indkk1,(jj-1)*4+4]]
    #                     end
    #
    #                 end
    #             end
    #         end
    #
    #         for kk1 = 1:maxmatpcut2 #maxmatpcut+1:maxmatpcut+maxmatpcut2
    #             for kk0 = 1:maxmatpcut2 #maxmatpcut+1:maxmatpcut+maxmatpcut2
    #                 for jj = 1:3
    #
    #                     # if Mpairdown3int[maxmatpcut+kk1,maxmatpcut+kk0,(jj-1)*4+1] != 0
    #                     #    fun_nudown[jjx,jjy] += conj(psi[maxmatpcut+kk1])*Mpairdown3float[maxmatpcut+kk1,maxmatpcut+kk0,jj]*psi[maxmatpcut+kk0]*phiHO[jjx,Mpairdown3int[maxmatpcut+kk1,maxmatpcut+kk0,(jj-1)*4+1]]*
    #                     #                                                                                                                           phiHO[jjx,Mpairdown3int[maxmatpcut+kk1,maxmatpcut+kk0,(jj-1)*4+3]]*
    #                     #                                                                                                                           phiHO[jjy,Mpairdown3int[maxmatpcut+kk1,maxmatpcut+kk0,(jj-1)*4+2]]*
    #                     #                                                                                                                           phiHO[jjy,Mpairdown3int[maxmatpcut+kk1,maxmatpcut+kk0,(jj-1)*4+4]]
    #                     # end
    #                     #
    #                     # if Mpairdown2up1int[kk1,kk0,(jj-1)*4+1] != 0
    #                     #    fun_nudu[jjx,jjy] += conj(psi[maxmatpcut+kk1])*Mpairdown2up1float[kk1,kk0,jj]*psi[maxmatpcut+kk0]*phiHO[jjx,Mpairdown2up1int[kk1,kk0,(jj-1)*4+1]]*
    #                     #                                                                                                      phiHO[jjx,Mpairdown2up1int[kk1,kk0,(jj-1)*4+3]]*
    #                     #                                                                                                      phiHO[jjy,Mpairdown2up1int[kk1,kk0,(jj-1)*4+2]]*
    #                     #                                                                                                      phiHO[jjy,Mpairdown2up1int[kk1,kk0,(jj-1)*4+4]]
    #                     # end
    #                     #
    #                     # if Mpairdown1up2int[kk1,kk0,(jj-1)*4+1] != 0
    #                     #    fun_nudu[jjx,jjy] += conj(psi[maxmatpcut+maxmatpcut2+kk1])*Mpairdown1up2float[kk1,kk0,jj]*psi[maxmatpcut+maxmatpcut2+kk0]*phiHO[jjx,Mpairdown1up2int[kk1,kk0,(jj-1)*4+1]]*
    #                     #                                                                                                                              phiHO[jjx,Mpairdown1up2int[kk1,kk0,(jj-1)*4+3]]*
    #                     #                                                                                                                              phiHO[jjy,Mpairdown1up2int[kk1,kk0,(jj-1)*4+2]]*
    #                     #                                                                                                                              phiHO[jjy,Mpairdown1up2int[kk1,kk0,(jj-1)*4+4]]
    #                     # end
    #                     #
    #                     # if Mpairup3int[kk1,kk0,(jj-1)*4+1] != 0
    #                     #    fun_nuup[jjx,jjy] += conj(psi[maxmatpcut+maxmatpcut2+kk1])*Mpairup3float[kk1,kk0,jj]*psi[maxmatpcut+maxmatpcut2+kk0]*phiHO[jjx,Mpairup3int[kk1,kk0,(jj-1)*4+1]]*
    #                     #                                                                                                                         phiHO[jjx,Mpairup3int[kk1,kk0,(jj-1)*4+3]]*
    #                     #                                                                                                                         phiHO[jjy,Mpairup3int[kk1,kk0,(jj-1)*4+2]]*
    #                     #                                                                                                                         phiHO[jjy,Mpairup3int[kk1,kk0,(jj-1)*4+4]]
    #                     # end
    #
    #                     indkk0 = kk0
    #                     indkk1 = kk1
    #
    #                     if kk1 >= kk0
    #                        indkk1 = kk0
    #                        indkk0 = kk1
    #                     end
    #
    #                     if Mpairdown3int[(maxmatpcut+indkk0-1)*(maxmatpcut+maxmatpcut2)+maxmatpcut+indkk1,(jj-1)*4+1] != 0
    #                        fun_nudown[jjx,jjy] += conj(psi[maxmatpcut+kk1])*Mpairdown3float[(maxmatpcut+indkk0-1)*(maxmatpcut+maxmatpcut2)+maxmatpcut+indkk1,jj]*psi[maxmatpcut+kk0]*phiHO[jjx,Mpairdown3int[(maxmatpcut+indkk0-1)*(maxmatpcut+maxmatpcut2)+maxmatpcut+indkk1,(jj-1)*4+1]]*
    #                                                                                                                                                                                  phiHO[jjx,Mpairdown3int[(maxmatpcut+indkk0-1)*(maxmatpcut+maxmatpcut2)+maxmatpcut+indkk1,(jj-1)*4+3]]*
    #                                                                                                                                                                                  phiHO[jjy,Mpairdown3int[(maxmatpcut+indkk0-1)*(maxmatpcut+maxmatpcut2)+maxmatpcut+indkk1,(jj-1)*4+2]]*
    #                                                                                                                                                                                  phiHO[jjy,Mpairdown3int[(maxmatpcut+indkk0-1)*(maxmatpcut+maxmatpcut2)+maxmatpcut+indkk1,(jj-1)*4+4]]
    #
    #                        fun_nuup[jjx,jjy] += conj(psi[maxmatpcut+maxmatpcut2+kk1])*Mpairdown3float[(maxmatpcut+indkk0-1)*(maxmatpcut+maxmatpcut2)+maxmatpcut+indkk1,jj]*psi[maxmatpcut+maxmatpcut2+kk0]*phiHO[jjx,Mpairdown3int[(maxmatpcut+indkk0-1)*(maxmatpcut+maxmatpcut2)+maxmatpcut+indkk1,(jj-1)*4+1]]*
    #                                                                                                                                                                                                        phiHO[jjx,Mpairdown3int[(maxmatpcut+indkk0-1)*(maxmatpcut+maxmatpcut2)+maxmatpcut+indkk1,(jj-1)*4+3]]*
    #                                                                                                                                                                                                        phiHO[jjy,Mpairdown3int[(maxmatpcut+indkk0-1)*(maxmatpcut+maxmatpcut2)+maxmatpcut+indkk1,(jj-1)*4+2]]*
    #                                                                                                                                                                                                        phiHO[jjy,Mpairdown3int[(maxmatpcut+indkk0-1)*(maxmatpcut+maxmatpcut2)+maxmatpcut+indkk1,(jj-1)*4+4]]
    #                     end
    #
    #                     if Mpairdown2up1int[kk1,kk0,(jj-1)*4+1] != 0
    #                        fun_nudu[jjx,jjy] += conj(psi[maxmatpcut+kk1])*Mpairdown2up1float[(kk0-1)*maxmatpcut2+kk1,jj]*psi[maxmatpcut+kk0]*phiHO[jjx,Mpairdown2up1int[(kk0-1)*maxmatpcut2+kk1,(jj-1)*4+1]]*
    #                                                                                                                                          phiHO[jjx,Mpairdown2up1int[(kk0-1)*maxmatpcut2+kk1,(jj-1)*4+3]]*
    #                                                                                                                                          phiHO[jjy,Mpairdown2up1int[(kk0-1)*maxmatpcut2+kk1,(jj-1)*4+2]]*
    #                                                                                                                                          phiHO[jjy,Mpairdown2up1int[(kk0-1)*maxmatpcut2+kk1,(jj-1)*4+4]]
    #
    #                        fun_nudu[jjx,jjy] += conj(psi[maxmatpcut+maxmatpcut2+kk1])*Mpairdown2up1float[(kk0-1)*maxmatpcut2+kk1,jj]*psi[maxmatpcut+maxmatpcut2+kk0]*phiHO[jjx,Mpairdown2up1int[(kk0-1)*maxmatpcut2+kk1,(jj-1)*4+1]]*
    #                                                                                                                                                                  phiHO[jjx,Mpairdown2up1int[(kk0-1)*maxmatpcut2+kk1,(jj-1)*4+3]]*
    #                                                                                                                                                                  phiHO[jjy,Mpairdown2up1int[(kk0-1)*maxmatpcut2+kk1,(jj-1)*4+2]]*
    #                                                                                                                                                                  phiHO[jjy,Mpairdown2up1int[(kk0-1)*maxmatpcut2+kk1,(jj-1)*4+4]]
    #                     end
    #
    #                 end
    #             end
    #         end
    #
    #         # for kk1 = 1:maxmatpcut2 #maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2*2
    #         #     for kk0 = 1:maxmatpcut2 #maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2*2
    #         #         for jj = 1:3
    #         #
    #         #             if Mpairdown1up2int[kk1,kk0,(jj-1)*4+1] != 0
    #         #                fun_nudu[jjx,jjy] += conj(psi[maxmatpcut+maxmatpcut2+kk1])*Mpairdown1up2float[kk1,kk0,jj]*psi[maxmatpcut+maxmatpcut2+kk0]*phiHO[jjx,Mpairdown1up2int[kk1,kk0,(jj-1)*4+1]]*
    #         #                                                                                                                                          phiHO[jjx,Mpairdown1up2int[kk1,kk0,(jj-1)*4+3]]*
    #         #                                                                                                                                          phiHO[jjy,Mpairdown1up2int[kk1,kk0,(jj-1)*4+2]]*
    #         #                                                                                                                                          phiHO[jjy,Mpairdown1up2int[kk1,kk0,(jj-1)*4+4]]
    #         #             end
    #         #
    #         #             if Mpairup3int[kk1,kk0,(jj-1)*4+1] != 0
    #         #                fun_nuup[jjx,jjy] += conj(psi[maxmatpcut+maxmatpcut2+kk1])*Mpairup3float[kk1,kk0,jj]*psi[maxmatpcut+maxmatpcut2+kk0]*phiHO[jjx,Mpairup3int[kk1,kk0,(jj-1)*4+1]]*
    #         #                                                                                                                                     phiHO[jjx,Mpairup3int[kk1,kk0,(jj-1)*4+3]]*
    #         #                                                                                                                                     phiHO[jjy,Mpairup3int[kk1,kk0,(jj-1)*4+2]]*
    #         #                                                                                                                                     phiHO[jjy,Mpairup3int[kk1,kk0,(jj-1)*4+4]]
    #         #             end
    #         #
    #         #         end
    #         #     end
    #         # end
    #
    #         # for kk1 = 1:maxmatpcut #maxmatpcut+maxmatpcut2*2+1:maxmatpcut+maxmatpcut2*2+maxmatpcut
    #         #     for kk0 = 1:maxmatpcut #maxmatpcut+maxmatpcut2*2+1:maxmatpcut+maxmatpcut2*2+maxmatpcut
    #         #         for jj = 1:3
    #         #             if Mpairup3int[maxmatpcut2+kk1,maxmatpcut2+kk0,(jj-1)*4+1] != 0
    #         #                fun_nuup[jjx,jjy] += conj(psi[maxmatpcut+maxmatpcut2*2+kk1])*Mpairup3float[maxmatpcut2+kk1,maxmatpcut2+kk0,jj]*psi[maxmatpcut+maxmatpcut2*2+kk0]*phiHO[jjx,Mpairup3int[maxmatpcut2+kk1,maxmatpcut2+kk0,(jj-1)*4+1]]*
    #         #                                                                                                                                                                 phiHO[jjx,Mpairup3int[maxmatpcut2+kk1,maxmatpcut2+kk0,(jj-1)*4+3]]*
    #         #                                                                                                                                                                 phiHO[jjy,Mpairup3int[maxmatpcut2+kk1,maxmatpcut2+kk0,(jj-1)*4+2]]*
    #         #                                                                                                                                                                 phiHO[jjy,Mpairup3int[maxmatpcut2+kk1,maxmatpcut2+kk0,(jj-1)*4+4]]
    #         #             end
    #         #         end
    #         #     end
    #         # end
    #
    #     end
    # end

    fun_nudown .= fun_nudown/Np/(Np-1)
    fun_nudu .= fun_nudu/Np/(Np-1)
    fun_nuup .= fun_nuup/Np/(Np-1)

    fun_nudown1 .= fun_nudown1/Np/(Np-1)
    fun_nuup1 .= fun_nuup1/Np/(Np-1)

    return fun_nudown, fun_nudu, fun_nuup, fun_nudown1, fun_nuup1

end
