include("in2bind.jl")
include("delta.jl")
include("coefficientonebodysummary.jl")

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

function epsilonsoc(ii::Int64,jj::Int64,commonii::Int64,commonjj::Int64,ksoc::Float64)

    # if abs(jj-ii) != 1
       # return 0.0
    # end

    # Note jj != ii
    epsilpn = sqrt(commonjj+1)
    epsilpn = epsilpn*sqrt(commonii+1)

    njj = jj - 1
    nii = ii - 1
    epsilpn = epsilpn*(-1)*ksoc/sqrt(2)*(sqrt(njj+1)*delta(nii-njj-1) - sqrt(njj)*delta(nii-njj+1)) # for down down
    # epsilpn = epsilpn*(1im)

    return epsilpn

end

function epsilonsoc2(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},common::Vector{Int64},ksoc::Float64)

   # bra and ket have to be the same except for one element, i.e. two elements have to be the same

   indNp = length(vecmbindnn)

   ind0 = 0
   ind1 = 0
   ind2 = 0
   common .= 0

   # for kk = 1:Np
   for kk = 1:indNp
       if vecmbindnn[kk] == vecmbindmm[kk]
          ind0 += 1
       elseif abs(vecmbindnn[kk] - vecmbindmm[kk]) == 1
          ind1 = kk
          ind2 += 1
       end
   end

   # if ind0 == Np-1 && ind2 == 1
   if ind0 == indNp-1 && ind2 == 1
      # common .= vecmbindnn[findall(x->x!=ind1,1:Np)]
      if indNp > 1
         common[1:indNp-1] .= vecmbindnn[findall(x->x!=ind1,1:indNp)]
      end
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

   return epsilonsoc(ii,jj,commonii,commonjj,1.0)

end

function epsilonW(ii::Int64,common::Vector{Int64},Omega::Float64)

    # a^+|n> = sqrt(n+1)|n+1>

    ind0 = 0

    for kk = 1:length(common)
        if ii == common[kk]
           ind0 += 1
        end
    end

    epsilpn = sqrt(ind0+1)*Omega/2

    return epsilpn

    # if ii == common # a^+|1> = sqrt(2)|2>
    #    epsilpn = sqrt(Np)
    # else
    #    epsilpn = 1.0
    # end
    #
    # epsilpn = epsilpn*Omega/2
    #
    # return epsilpn

end

function epsilonW0_3(ii::Int64,common::Vector{Int64},Omega::Float64)

   # ket nn is up up down
   # bra mm is down down up

    # a|n> = sqrt(n)|n-1>
    # a^+|n> = sqrt(n+1)|n+1>

    epsilpn = 1.0

    for kk = 1:length(common)
       ind0 = 0
       if ii == common[kk]
          ind0 += 1
       end
       epsilpn = sqrt(ind0+1)*epsilpn
    end

    # ind0 = 0
    # if ii == common[1]
    #    ind0 += 1
    # end
    # epsilpn = sqrt(ind0+1)
    #
    # ind0 = 0
    # if ii == common[2]
    #    ind0 += 1
    # end
    # epsilpn = sqrt(ind0+1)*epsilpn

    epsilpn = Omega/2*epsilpn

    return epsilpn

end

# function epsilonW0_3(ii::Int64,common::Vector{Int64},Np::Int64,Omega::Float64)
#
#     # a|n> = sqrt(n)|n-1>
#     # a^+|n> = sqrt(n+1)|n+1>
#
#     epsilpn = 1.0
#
#     for kk = 1:length(common)
#        ind0 = 0
#        if ii == common[kk]
#           ind0 += 1
#        end
#        epsilpn = sqrt(ind0+1)*epsilpn
#     end
#
#     # ind0 = 0
#     # if ii == common[1]
#     #    ind0 += 1
#     # end
#     # epsilpn = sqrt(ind0+1)
#     #
#     # ind0 = 0
#     # if ii == common[2]
#     #    ind0 += 1
#     # end
#     # epsilpn = sqrt(ind0+1)*epsilpn
#
#     epsilpn = Omega/2*epsilpn
#
#     return epsilpn
#
# end

function epsilonW2(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},common::Vector{Int64},Np::Int64,Omega::Float64)

   # ket nn is down down up
   # bra mm is down down down

   # vecmbindnn2 is down down up, down up down, or up down down

   # if vecmbindnn2 != vecmbindmm
   #    return 0.0
   # end

   ind0 = 0
   common .= 0 # zeros(Int64,Np-1)

   for kk = 1:Np-1
       for ll = ind0+1:Np
           if vecmbindnn[kk] == vecmbindmm[ll]
              common[kk] = vecmbindnn[kk]
              ind0 = ll
              break
           end
       end
   end

   ii = vecmbindnn[end]

   return epsilonW(ii,common,1.0)

   # if vecmbindnn[1] == vecmbindmm[1] && vecmbindnn[2] == vecmbindmm[2]
   #
   #    common = vecmbindnn[1]
   #    # jj = vecmbindnn[2]
   #    ii = vecmbindmm[2]
   #    HW[mm,maxmatpcut+nn] = epsilonW(ii,common,Np,1.0)
   #    break
   #
   # elseif vecmbindnn[1] == vecmbindmm[2] && vecmbindnn[2] == vecmbindmm[1]
   #
   #    common = vecmbindnn[1]
   #    # jj = vecmbindnn[2]
   #    ii = vecmbindmm[1]
   #    HW[mm,maxmatpcut+nn] = epsilonW(ii,common,Np,1.0)
   #    break
   #
   # end

end

function epsilonW3(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},common::Vector{Int64},Omega::Float64)

   # ket nn is up up down
   # bra mm is down down up

   ind0 = 0
   common .= 0 # zeros(Int64,Np-1)

   if vecmbindnn[1] == vecmbindmm[3]

      if vecmbindnn[3] == vecmbindmm[1]

         if vecmbindnn[2] == vecmbindmm[2]
            ii = vecmbindnn[2]
            common[1] = vecmbindnn[3] # down
            common[2] = vecmbindnn[1] # up
            return epsilonW0_3(ii,common,1.0)
         end

      elseif vecmbindnn[3] == vecmbindmm[2]

         if vecmbindnn[2] == vecmbindmm[1]
            ii = vecmbindnn[2]
            common[1] = vecmbindnn[3] # down
            common[2] = vecmbindnn[1] # up
            return epsilonW0_3(ii,common,1.0)
         end

      end

   elseif vecmbindnn[2] == vecmbindmm[3]

      if vecmbindnn[3] == vecmbindmm[1]

         if vecmbindnn[1] == vecmbindmm[2]
            ii = vecmbindnn[1]
            common[1] = vecmbindnn[3] # down
            common[2] = vecmbindnn[2] # up
            return epsilonW0_3(ii,common,1.0)
         end

      elseif vecmbindnn[3] == vecmbindmm[2]

         if vecmbindnn[1] == vecmbindmm[1]
            ii = vecmbindnn[1]
            common[1] = vecmbindnn[3] # down
            common[1] = vecmbindnn[2] # up
            return epsilonW0_3(ii,common,1.0)
         end

      end

   end

   return 0.0

end

# function epsilonW3(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},common::Vector{Int64},Np::Int64,Omega::Float64)
#
#    # ket nn is down up up
#    # bra mm is down down up
#
#    ind0 = 0
#    common .= 0 # zeros(Int64,Np-1)
#
#    # for ll = 1:2
#    #     if vecmbindnn[1] == vecmbindmm[ll]
#    #        for qq = 2:3
#    #            if vecmbindnn[qq] == vecmbindmm[3]
#    #               if vecmbindnn[3] == vecmbindmm[2]
#    #               end
#    #               break
#    #            end
#    #            break
#    #        end
#    #        break
#    #     end
#    # end
#
#    if vecmbindnn[1] == vecmbindmm[1]
#
#       if vecmbindnn[2] == vecmbindmm[3]
#
#          if vecmbindnn[3] == vecmbindmm[2]
#             ii = vecmbindnn[3]
#             common[1] = vecmbindnn[2]
#             common[2] = vecmbindnn[1]
#             return epsilonW0_3(ii,common,Np,1.0)
#          end
#
#       elseif vecmbindnn[3] == vecmbindmm[3]
#
#          if vecmbindnn[2] == vecmbindmm[2]
#             ii = vecmbindnn[2]
#             common[1] = vecmbindnn[3]
#             common[2] = vecmbindnn[1]
#             return epsilonW0_3(ii,common,Np,1.0)
#          end
#
#       end
#
#    elseif vecmbindnn[1] == vecmbindmm[2]
#
#       if vecmbindnn[2] == vecmbindmm[3]
#
#          if vecmbindnn[3] == vecmbindmm[1]
#             ii = vecmbindnn[3]
#             common[1] = vecmbindnn[2]
#             common[2] = vecmbindnn[1]
#             return epsilonW0_3(ii,common,Np,1.0)
#          end
#
#       elseif vecmbindnn[3] == vecmbindmm[3]
#
#          if vecmbindnn[2] == vecmbindmm[1]
#             ii = vecmbindnn[2]
#             common[1] = vecmbindnn[3]
#             common[1] = vecmbindnn[2]
#             return epsilonW0_3(ii,common,Np,1.0)
#          end
#
#       end
#
#    end
#
#    return 0.0
#
# end

# function epsilonW4(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},common::Vector{Int64},Np::Int64,Omega::Float64)
#
#    # ket nn is down up up
#    # bra mm is up up up
#
#    ind0 = 0
#    common .= 0 # zeros(Int64,Np-1)
#
#    for kk = 2:Np
#        for ll = ind0+1:Np
#            if vecmbindnn[kk] == vecmbindmm[ll]
#               common[kk-1] = ll
#               ind0 = ll
#               break
#            end
#        end
#    end
#
#    ii = vecmbindnn[1]
#
#    return epsilonW(ii,common,Np,1.0)
#
# end

function Hsocfunccutoffk1W1!(indvec::Vector{Int64}, indvec2::Vector{Int64}, Msize0::Int64, Np::Int64, matp::Matrix{Int64}, matp20::Matrix{Int64}, matp21::Matrix{Int64}, Hho::SparseMatrixCSC{Float64}, Hsoc::SparseMatrixCSC{Float64}, HW::SparseMatrixCSC{Float64}, matonebody::SparseMatrixCSC{Float64})

    maxmatpcut = length(indvec)
    maxmatpcut2 = length(indvec2)
    # maxmatpcut3 = length(indvec3) # maxmatpcut3 == maxmatpcut2

    # reset matrices
    Hho .= 0
    Hsoc .= 0
    HW .= 0

    # define vectors
    vecmbindnn = zeros(Int64,Np)
    vecmbindmm = zeros(Int64,Np)
    common = zeros(Int64,Np-1)
    vecmbindnn3 = zeros(Int64,2)
    vecmbindmm3 = zeros(Int64,2)
    vecindcoeff = zeros(Float64,3,2)

    # define a matrix for the Hamiltonian for down down down
    for nn = 1:maxmatpcut # parfor

        # ket
        in2bind!(indvec[nn],Msize0,Np,matp,vecmbindnn)

        # Hho
        Hho[nn,nn] = sum(vecmbindnn[:]) - Np/2
        # (vecmbindnn[1]-1+1/2) + (vecmbindnn[2]-1+1/2) + (vecmbindnn[3]-1+1/2)

        # one body
        matonebody[nn,nn] = sum(vecmbindnn[:]) - Np/2

        for mm = 1:nn-1

            # bra
            in2bind!(indvec[mm],Msize0,Np,matp,vecmbindmm)

            # Hsoc
            Hsoc[mm,nn] = epsilonsoc2(vecmbindnn,vecmbindmm,common,1.0)

            # one body
            matonebody[mm,nn] = coefficientonebody(vecmbindnn,vecmbindmm,vecmbindnn3,vecmbindmm3,vecindcoeff,Np)

        end

    end

    # define the Hamiltonian for up up up
    Hho[end-(maxmatpcut-1):end,end-(maxmatpcut-1):end] = Hho[1:maxmatpcut,1:maxmatpcut]
    Hsoc[end-(maxmatpcut-1):end,end-(maxmatpcut-1):end] = Hsoc[1:maxmatpcut,1:maxmatpcut]*(-1) # due to up up

    # one body
    matonebody[end-(maxmatpcut-1):end,end-(maxmatpcut-1):end] = matonebody[1:maxmatpcut,1:maxmatpcut]

    # define the Hamiltonian for down down up
    vecmbindnn2 = zeros(Int64,Np)
    vecmbindmm2 = zeros(Int64,Np)
    maxmatp21 = matp21[Msize0+1,1+1]

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
        in2bind!(indketdown,Msize0,Np-1,matp20,vecmbindnn2)
        vecmbindnn[1:Np-1] = vecmbindnn2[1:Np-1]
        in2bind!(indketup,Msize0,1,matp21,vecmbindnn2)
        vecmbindnn[Np:Np] = vecmbindnn2[1:1]

        # Hho
        Hho[maxmatpcut+nn,maxmatpcut+nn] = sum(vecmbindnn[:]) - Np/2

        # one body
        matonebody[nn,nn] = sum(vecmbindnn[:]) - Np/2

        for mm = 1:nn-1

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

            # Hsoc
            if vecmbindnn[Np:Np] == vecmbindmm[Np:Np]
               Hsoc[maxmatpcut+mm,maxmatpcut+nn] = epsilonsoc2(vecmbindnn[1:Np-1],vecmbindmm[1:Np-1],common,1.0) # down down
            elseif vecmbindnn[1:Np-1] == vecmbindmm[1:Np-1]
               Hsoc[maxmatpcut+mm,maxmatpcut+nn] = epsilonsoc2(vecmbindnn[Np:Np],vecmbindmm[Np:Np],common,1.0)*(-1) # up up
            end

            # one body
            matonebody[mm,nn] = coefficientonebody2(vecmbindnn,vecmbindmm,vecmbindnn3,vecmbindmm3,vecindcoeff,Np)

        end

    end

    # define the Hamiltonian for up up down
    Hho[maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2,maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2] = Hho[maxmatpcut+1:maxmatpcut+maxmatpcut2,maxmatpcut+1:maxmatpcut+maxmatpcut2]
    Hsoc[maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2,maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2] = Hsoc[maxmatpcut+1:maxmatpcut+maxmatpcut2,maxmatpcut+1:maxmatpcut+maxmatpcut2]*(-1) # due to up up down

    # one body
    matonebody[maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2,maxmatpcut+maxmatpcut2+1:maxmatpcut+maxmatpcut2+maxmatpcut2] = matonebody[maxmatpcut+1:maxmatpcut+maxmatpcut2,maxmatpcut+1:maxmatpcut+maxmatpcut2]

    # # define the Hamiltonian for down up up
    # maxmatp31 = matp31[Msize0+1,2+1]
    # indNp = 2
    #
    # for nn = 1:maxmatpcut2
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
    #     # Hho
    #     Hho[maxmatpcut+maxmatpcut2+nn,maxmatpcut+maxmatpcut2+nn] = sum(vecmbindnn[:]) - Np/2
    #
    #     for mm = 1:nn-1
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
    #         # Hsoc
    #         if vecmbindnn[Np-indNp+1:Np] == vecmbindmm[Np-indNp+1:Np]
    #            Hsoc[maxmatpcut+mm,maxmatpcut+nn] = epsilonsoc2(vecmbindnn[1:Np-indNp],vecmbindmm[1:Np-indNp],common[1:Np-indNp],Np-indNp,1.0) # down
    #         elseif vecmbindnn[1:Np-indNp] == vecmbindmm[1:Np-indNp]
    #            Hsoc[maxmatpcut+mm,maxmatpcut+nn] = epsilonsoc2(vecmbindnn[Np-indNp+1:Np],vecmbindmm[Np-indNp+1:Np],common[Np-indNp+1:Np],indNp,1.0)*(-1) # up
    #         end
    #
    #     end
    #
    # end

    # HW for <down,down,down|down,down,up>
    for nn = 1:maxmatpcut2

        # ket having mixed spins
        # indvec2[nn] = (nndown-1)*maxmatp21 + nnup
        indketup = mod(indvec2[nn],maxmatp21)
        if indketup != 0
           indketdown = div(indvec2[nn],maxmatp21) + 1
        else # indketup == 0
           indketdown = div(indvec2[nn],maxmatp21)
           indketup = maxmatp21
        end
        in2bind!(indketdown,Msize0,Np-1,matp20,vecmbindnn2)
        vecmbindnn[1:Np-1] = vecmbindnn2[1:Np-1] # down
        in2bind!(indketup,Msize0,1,matp21,vecmbindnn2)
        vecmbindnn[Np:Np] = vecmbindnn2[1:1] # up
        vecmbindnn2 .= vecmbindnn

        # sort vecmbindnn2 so note vecmbindnn can be not only down down up, down up down, up down down
        sort!(vecmbindnn2)

        for mm = 1:maxmatpcut

            # bra having down down down
            in2bind!(indvec[mm],Msize0,Np,matp,vecmbindmm)

            # HW
            if vecmbindnn2 == vecmbindmm
               HW[mm,maxmatpcut+nn] = epsilonW2(vecmbindnn,vecmbindmm,common,Np,1.0)
            end

            # one body
            matonebody[mm,nn] = coefficientonebody3(vecmbindnn,vecmbindmm,vecmbindnn3,vecmbindmm3,vecindcoeff,Np)

        end

    end

    # HW for <up,up,up|up,up,down>
    HW[end-(maxmatpcut-1):end,end-(maxmatpcut+maxmatpcut2-1):end-maxmatpcut] = HW[1:maxmatpcut,maxmatpcut+1:maxmatpcut+maxmatpcut2]

    # one body
    matonebody[end-(maxmatpcut-1):end,end-(maxmatpcut+maxmatpcut2-1):end-maxmatpcut] = matonebody[1:maxmatpcut,maxmatpcut+1:maxmatpcut+maxmatpcut2]

    # HW for <down,down,up|up,up,down>
    indNp = 1
    for nn = 1:maxmatpcut2

        # ket having up up down
        # indvec2[nn] = (nndown-1)*maxmatp21 + nnup
        indketup = mod(indvec2[nn],maxmatp21)
        if indketup != 0
           indketdown = div(indvec2[nn],maxmatp21) + 1
        else # indketup == 0
           indketdown = div(indvec2[nn],maxmatp21)
           indketup = maxmatp21
        end
        in2bind!(indketdown,Msize0,Np-indNp,matp20,vecmbindnn2)
        vecmbindnn[1:Np-indNp] = vecmbindnn2[1:Np-indNp] # up
        in2bind!(indketup,Msize0,indNp,matp21,vecmbindnn2)
        vecmbindnn[Np-indNp+1:Np] = vecmbindnn2[1:indNp] # down
        # vecmbindnn2 .= vecmbindnn
        # sort!(vecmbindnn2)

        for mm = 1:maxmatpcut2

            # bra having down down up
            # indvec2[nn] = (nndown-1)*maxmatp21 + nnup
            indbraup = mod(indvec2[mm],maxmatp21)
            if indbraup != 0
               indbradown = div(indvec2[mm],maxmatp21) + 1
            else # indketup == 0
               indbradown = div(indvec2[mm],maxmatp21)
               indbraup = maxmatp21
            end
            in2bind!(indbradown,Msize0,Np-indNp,matp20,vecmbindmm2)
            vecmbindmm[1:Np-indNp] = vecmbindmm2[1:Np-indNp] # down
            in2bind!(indbraup,Msize0,indNp,matp21,vecmbindmm2)
            vecmbindmm[Np-indNp+1:Np] = vecmbindmm2[1:indNp] # up
            # vecmbindmm2 .= vecmbindmm
            # sort!(vecmbindmm2)

            # HW
            HW[maxmatpcut+mm,maxmatpcut+maxmatpcut2+nn] = epsilonW3(vecmbindnn,vecmbindmm,common,1.0)
            # if vecmbindnn2 == vecmbindmm2
            #    HW[maxmatpcut+mm,maxmatpcut+maxmatpcut2+nn] = epsilonW3(vecmbindnn,vecmbindmm,common,Np,1.0)
            # end

            # one body
            matonebody[mm,nn] = coefficientonebody4(vecmbindnn,vecmbindmm,vecmbindnn3,vecmbindmm3,vecindcoeff,Np)

        end

    end

    # # HW for <up,up,up|down,up,up>
    # indNp = 2
    #
    # for nn = 1:maxmatpcut2
    #
    #     # ket having mixed spins
    #     # indvec3[nn] = (nndown-1)*maxmatp31 + nnup
    #     indketup = mod(indvec3[nn],maxmatp31)
    #     if indketup != 0
    #        indketdown = div(indvec3[nn],maxmatp31) + 1
    #     else # indketup == 0
    #        indketdown = div(indvec3[nn],maxmatp31)
    #        indketup = maxmatp31
    #     end
    #
    #     in2bind!(indketdown,Msize0,Np-indNp,matp30,vecmbindnn2)
    #     vecmbindnn[1:Np-indNp] = vecmbindnn2[1:Np-indNp] # down
    #     in2bind!(indketup,Msize0,indNp,matp31,vecmbindnn2)
    #     vecmbindnn[Np-indNp+1:Np] = vecmbindnn2[Np-indNp+1:Np] # up
    #     vecmbindnn2 .= vecmbindnn
    #
    #     # sort vecmbindnn2
    #     sort!(vecmbindnn2)
    #
    #     for mm = 1:maxmatpcut
    #
    #         # bra having up up up
    #         in2bind!(indvec[mm],Msize0,Np,matp,vecmbindmm)
    #
    #         # HW
    #         if vecmbindnn2 == vecmbindmm
    #            HW[end-(maxmatpcut-mm),maxmatpcut+maxmatpcut2+nn] = epsilonW4(vecmbindnn,vecmbindmm,common[1:Np-1],Np,1.0)
    #         end
    #
    #     end
    #
    # end

    # since the Hamiltonian is helmitian
    # Hsoc .= Hsoc + Hsoc' - spdiagm(diag(Hsoc))
    Hsoc .= Hsoc + (-1)*Hsoc' # Hsco is an imaginary off-diagonal matrix for the reality
    # HW .= HW + HW' - spdiagm(diag(HW))
    HW .= HW + HW' # HW is an off-diagonal matrix

    # one body
    matonebody .= matonebody + matonebody' - spdiagm(diag(matonebody))

end
