include("in2bind.jl")
include("delta.jl")

# function check_duplicates(vecmbindnn::Vector{Int64},Np::Int64,hist::Vector{Int64},Msize0::Int64)
#
#     if Np != 3
#        error("This function is for Np=3 but can be extended for any particle number.")
#     end
#
#     hist .= 0 # hist = zeros(Int64,Msize0)
#
#     for jj = 1:Np
#         hist[vecmbindnn[jj]] += 1
#     end
#
#     for jj = 1:Msize0
#         if hist[jj] == Np
#            return -1 #1
#         elseif hist[jj] == Np-1
#            return jj #2
#         end
#     end
#
#     return 0 # if hist[jj] == Np-2 = 1
#
# end

# function epsilonHO(vecmbindnn::Vector{Int64},Np::Int64,hist::Vector{Int64},Msize0::Int64)
#
#     # This function is for Np=3 but can be extended for any particle number
#
#     case = check_duplicates(vecmbindnn,Np,hist,Msize0)
#
#     if case == Msize0+1
#        return (vecmbindnn[1]-1+1/2)*Np
#     elseif case != 0 # case <= Msize0 && case > 0
#        return (vecmbindnn[case]-1+1/2)*(Np-1) + (vecmbindnn[2]-1+1/2)
#     else
#        return (vecmbindnn[1]-1+1/2) + (vecmbindnn[2]-1+1/2) + (vecmbindnn[3]-1+1/2)
#     end
#
#     # if vecmbindnn[1] == vecmbindnn[2]
#     #    Hho[nn,nn] = (vecmbindnn[1]-1+1/2)*Np
#     # else
#     #    Hho[nn,nn] = (vecmbindnn[1]-1+1/2) + (vecmbindnn[2]-1+1/2)
#     # end
#
# end

# function check_duplicates_soc(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},Np::Int64,hist::Vector{Int64},Msize0::Int64)
#
#     if Np != 3
#        error("This function is for Np=3 but can be extended for any particle number.")
#     end
#
#     hist .= 0 # hist = zeros(Int64,Msize0)
#
#     for jj = 1:Np
#         hist[vecmbindnn[jj]] += 1
#     end
#
#     for jj = 1:Msize0
#         if hist[jj] == Np
#            return -1 #1
#         elseif hist[jj] == Np-1
#            return jj #2
#         end
#     end
#
#     return 0 # if hist[jj] == Np-2 = 1
#
# end

function epsilonsoc(ii::Int64,jj::Int64,commonii::Int64,commonjj::Int64,Np::Int64,ksoc::Float64)

    # if abs(jj-ii) != 1
       # return 0.0
    # end

    # Note jj != ii
    epsilpn = sqrt(commonjj+1)
    epsilpn = epsilpn*sqrt(commonii+1)

    njj = jj - 1
    nii = ii - 1
    epsilpn = epsilpn*(-1im)*ksoc/sqrt(2)*(sqrt(njj+1)*delta(nii-njj-1) - sqrt(njj)*delta(nii-njj+1)) # for down down

    return epsilpn

end

function epsilonsoc2(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},common::Vector{Int64},Np::Int64,ksoc::Float64)

   # bra and ket have to be the same except for one element, i.e. two elements have to be the same

   # vecmbindnn_tot = sum(vecmbindnn)
   # vecmbindmm_tot = sum(vecmbindmm)
   # if abs(vecmbindnn_tot - vecmbindmm_tot) == 1 # avoid some states very diffrent from each other (but cannot be all of them)
   #    return 0.0
   # end

   ind0 = 0
   ind1 = 0
   ind2 = 0
   # common .= 0

   for kk = 1:Np
       if vecmbindnn[kk] == vecmbindmm[kk]
          ind0 += 1
       elseif abs(vecmbindnn[kk] - vecmbindmm[kk]) == 1
          ind1 = kk
          ind2 += 1
       end
   end

   if ind0 == Np-1 && ind2 == 1
      common .= vecmbindnn[findall(x->x!=ind1,1:Np)]
      jj = vecmbindnn[ind1]
      ii = vecmbindmm[ind1]
   else
      return 0.0
   end

   a = findlast(x->x==jj,common)
   if !isa(a,Nothing)
      commonjj = a-findfirst(x->x==jj,common) + 1
   else
      commonjj = 0
   end

   a = findlast(x->x==ii,common)
   if !isa(a,Nothing)
      commonii = a-findfirst(x->x==ii,common) + 1
   else
      commonii = 0
   end

   return epsilonsoc(ii,jj,commonii,commonjj,Np,1.0)

end

function epsilonW(ii::Int64,common::Int64,Np::Int64,Omega::Float64)

    # a|1> = |0>

    if ii == common # a^+|1> = sqrt(2)|2>
       epsilpn = sqrt(Np)
    else
       epsilpn = 1.0
    end

    epsilpn = epsilpn*Omega/2

    return epsilpn

end

function Hsocfunccutoffk1W1!(indvec::Vector{Int64}, indvec2::Vector{Int64}, Msize0::Int64, Np::Int64, matp::Matrix{Int64}, matp20::Matrix{Int64}, matp21::Matrix{Int64}, Hho::SparseMatrixCSC{Float64}, Hsoc::SparseMatrixCSC{ComplexF64}, HW::SparseMatrixCSC{Float64})

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
    common = zeros(Int64,Np-1)

    # define a matrix for the Hamiltonian for down down
    for nn = 1:maxmatpcut # parfor

        # ket
        in2bind!(indvec[nn],Msize0,Np,matp,vecmbindnn)

        # bra
        # in2bind!(indvec[mm],Msize0,Np,matp,vecmbindmm)
        # vecmbindmm .= vecmbindnn

        # Hho for Np=3
        Hho[nn,nn] = sum(vecmbindnn[:]) - Np/2
        # (vecmbindnn[1]-1+1/2) + (vecmbindnn[2]-1+1/2) + (vecmbindnn[3]-1+1/2)
        # epsilonHO(vecmbindnn,Np,hist,Msize0)

        for mm = 1:nn-1

            # bra
            in2bind!(indvec[mm],Msize0,Np,matp,vecmbindmm)

            # Hsoc
            Hsoc[mm,nn] = epsilonsoc2(vecmbindnn,vecmbindmm,common,Np,1.0)

        end

    end

    # define the Hamiltonian for up up
    Hho[end-(maxmatpcut-1):end,end-(maxmatpcut-1):end] = Hho[1:maxmatpcut,1:maxmatpcut]
    Hsoc[end-(maxmatpcut-1):end,end-(maxmatpcut-1):end] = Hsoc[1:maxmatpcut,1:maxmatpcut]*(-1) # due to up up

    # define the Hamiltonian for down up
    vecmbindnn2 = zeros(Int64,Np)

    for nn = 1:maxmatpcut2

        # ket
        # indvec2[nn] = (nndown-1)*maxmatp21 + nnup
        indketup = mod(indvec2[nn],maxmatp21)
        if indketup != 0
           indketdown = div(indvec2[nn],maxmatp21) + 1
        else # indketup == 0
           indketdown = div(indvec2[nn],maxmatp21)
           indketup = maxmatp21
        end

        # indket2 = mod(indvec2[nn],Msize0)
        # if indket2 != 0
        #    indket1 = div(indvec2[nn],Msize0)+1
        # else
        #    indket1 = div(indvec2[nn],Msize0)
        #    indket2 = Msize0
        # end

        in2bind!(indketdown,Msize0,Np-1,matp20,vecmbindnn2)
        vecmbindnn[1:Np-1] = vecmbindnn2[1:Np-1]
        in2bind!(indketup,Msize0,1,matp20,vecmbindnn2)
        vecmbindnn[Np:Np] = vecmbindnn2[1:1]

        # Hho
        Hho[maxmatpcut+nn,maxmatpcut+nn] = sum(vecmbindnn[:]) - Np/2

        ##########

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

    # HW for <down,down|mixed>
    for nn = 1:maxmatpcut2

        # ket having mixed spins
        indket2 = mod(indvec2[nn],Msize0)
        if indket2 != 0
           indket1 = div(indvec2[nn],Msize0)+1
        else
           indket1 = div(indvec2[nn],Msize0)
           indket2 = Msize0
        end
        in2bind!(indket1,Msize0,Np-1,matp2,vecmbindnn2)
        vecmbindnn[1] = vecmbindnn2[1] # down
        in2bind!(indket2,Msize0,Np-1,matp2,vecmbindnn2)
        vecmbindnn[2] = vecmbindnn2[1] # up

        for mm = 1:maxmatpcut

            # bra having down down
            in2bind!(indvec[mm],Msize0,Np,matp,vecmbindmm)

            # HW
            if vecmbindnn[1] == vecmbindmm[1] && vecmbindnn[2] == vecmbindmm[2]

               common = vecmbindnn[1]
               # jj = vecmbindnn[2]
               ii = vecmbindmm[2]
               HW[mm,maxmatpcut+nn] = epsilonW(ii,common,Np,1.0)
               break

            elseif vecmbindnn[1] == vecmbindmm[2] && vecmbindnn[2] == vecmbindmm[1]

               common = vecmbindnn[1]
               # jj = vecmbindnn[2]
               ii = vecmbindmm[1]
               HW[mm,maxmatpcut+nn] = epsilonW(ii,common,Np,1.0)
               break

            end

        end

    end

    # HW for <up,up|mixed>

    # But the construction given below is the same as <down,down|mixed>, and so skip it and copy HW for <down,down|mixed>.
    HW[maxmatpcut+maxmatpcut2+1:end,maxmatpcut+1:maxmatpcut+maxmatpcut2] = HW[1:maxmatpcut,maxmatpcut+1:maxmatpcut+maxmatpcut2]

    # for nn = 1:maxmatpcut2
    #
    #     # ket having mixed spins
    #     indket2 = mod(indvec2[nn],Msize0)
    #     if indket2 != 0
    #        indket1 = div(indvec2[nn],Msize0)+1
    #     else
    #        indket1 = div(indvec2[nn],Msize0)
    #        indket2 = Msize0
    #     end
    #     in2bind!(indket1,Msize0,Np-1,matp2,vecmbindnn2)
    #     vecmbindnn[1] = vecmbindnn2[1] # down
    #     in2bind!(indket2,Msize0,Np-1,matp2,vecmbindnn2)
    #     vecmbindnn[2] = vecmbindnn2[1] # up
    #
    #     for mm = 1:maxmatpcut
    #
    #         # bra having up up
    #         in2bind!(indvec[mm],Msize0,Np,matp,vecmbindmm)
    #
    #         # HW
    #         if vecmbindnn[1] == vecmbindmm[1] && vecmbindnn[2] == vecmbindmm[2]
    #
    #            common = vecmbindnn[2]
    #            # jj = vecmbindnn[1]
    #            ii = vecmbindmm[1]
    #            HW[maxmatpcut+maxmatpcut2+mm,maxmatpcut+nn] = epsilonW(ii,common,Np,Omega)
    #
    #         elseif vecmbindnn[1] == vecmbindmm[2] && vecmbindnn[2] == vecmbindmm[1]
    #
    #            common = vecmbindnn[2]
    #            # jj = vecmbindnn[1]
    #            ii = vecmbindmm[2]
    #            HW[maxmatpcut+maxmatpcut2+mm,maxmatpcut+nn] = epsilonW(ii,common,Np,Omega)
    #
    #         end
    #
    #     end
    #
    # end

    # since the Hamiltonian is helmitian
    Hsoc .= Hsoc + Hsoc' - spdiagm(diag(Hsoc))
    HW .= HW + HW' - spdiagm(diag(HW))

end
