
function coefficientonebody1!(ii::Int64, jj::Int64, psi::Vector{ComplexF64}, Msize0::Int64, maxmatpcut::Int64)

   for nn = 1:maxmatpcut # down down down

       # ket
       in2bind!(indvec[nn],Msize0,Np,matp,vecmbindnn)

       check0 = findall(x->x==jj,vecmbindnn)
       if length(check0) == 1
          psi[nn] = 1.0
       elseif length(check0) == 2
          psi[nn] = sqrt(2)
       elseif length(check0) == 3
          psi[nn] = sqrt(3)
       # else # length(check0) == 0
          # psi[nn] = 0.0
       end

   end

   for nn = 1:maxmatpcut2 # down down up

       # ket
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

       check0 = findall(x->x==jj,vecmbindnn[1:Np-1])
       if length(check0) == 1
          psi[maxmatpcut+nn] = 1.0
       elseif length(check0) == 2
          psi[maxmatpcut+nn] = sqrt(2)
       # else # length(check0) == 0
          # psi[maxmatpcut+nn] = 0.0
       end

   end

   for nn = 1:maxmatpcut2 # up up down

       # ket
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

       if vecmbindnn[Np:Np] == jj
          maxmatpcut[maxmatpcut+maxmatpcut2+nn] = 1.0
       # else # vecmbindnn[Np:Np] != jj
          # maxmatpcut[maxmatpcut+maxmatpcut2+nn] = 0.0
       end

   end

   # maxmatpcut[maxmatpcut+maxmatpcut2++maxmatpcut2+1:end] = 0.0 # up up up

end

function coefficientonebody(ii::Int64, jj::Int64, psi::Vector{Float64}, Msize0::Int64, maxmatpcut::Int64)

   if jj <= Msize0
      coefficientonebody1!(jj,psi)
   else # jj >= Msize0+1
      coefficientonebody2!(jj-Msize0,psi)
   end




   if jj <= Msize0
      if ii <= Msize0 # down down
         output = coefficientonebody1(ii,jj,psi)
      else # ii >= Msize0+1 # up down
         output = coefficientonebody2(ii-Msize0,jj,psi)
      end
   else # jj >= Msize0+1
      if ii <= Msize0 # down up
         output = coefficientonebody3(ii,jj-Msize0,psi)
      else # ii >= Msize0+1 # up up
         output = coefficientonebody4(ii-Msize0,jj-Msize0,psi)
      end
   end

end



function coefficientonebody0(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},vecmbindnn3::Vector{Int64},vecmbindmm3::Vector{Int64},Np::Int64)

    # vecmbindnn3 = zeros(Int64,2)
    # vecmbindmm3 = zeros(Int64,2)

    # a|n>= sqrt(n)|n-1>
    # a^+|n>= sqrt(n+1)|n+1>

    # initilise
    vecmbindnn3 .= 0
    vecmbindmm3 .= 0

    for pp = 1:Np

        vecmbindnn3 .= vecmbindnn[findall(x->x!=pp,1:Np)]

        for qq = 1:Np

            vecmbindmm3 .= vecmbindmm[findall(x->x!=qq,1:Np)]

            if vecmbindnn3 == vecmbindmm3

               element = 1.0

               check0 = findall(x->x==vecmbindnn[pp],vecmbindnn3)
               if length(check0) == 0
                  element = element*1.0
               elseif length(check0) == 1
                  element = element*sqrt(2)
               elseif length(check0) == 2
                  element = element*sqrt(3)
               end

               check1 = findall(x->x==vecmbindmm[qq],vecmbindmm3)
               if length(check1) == 0
                  element = element*1.0
               elseif length(check1) == 1
                  element = element*sqrt(2)
               elseif length(check1) == 2
                  element = element*sqrt(3)
               end

               return element

            end

        end
    end

    return 0.0

end

function coefficientonebody2(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},vecmbindnn3::Vector{Int64},vecmbindmm3::Vector{Int64},Np::Int64)

    # vecmbindnn3 = zeros(Int64,2)
    # vecmbindmm3 = zeros(Int64,2)

    # a|n>= sqrt(n)|n-1>
    # a^+|n>= sqrt(n+1)|n+1>

    # initilise
    vecmbindnn3 .= 0
    vecmbindmm3 .= 0

    if vecmbindnn[1:Np-1] == vecmbindmm[1:Np-1]
       element = 1.0
    end


       vecmbindnn3 .= vecmbindnn[1:Np-1]
       vecmbindmm3 .= vecmbindmm[1:Np-1]

       for pp = 1:Np-1

           for qq = 1:Np-1

               if vecmbindnn3[pp] == vecmbindmm3[qq]

                  check0 = findall(x->x==vecmbindnn3[pp],vecmbindnn3)
          if length(check0) == 0
             element = element*1.0
          elseif length(check0) == 1
             element = element*sqrt(2)
          elseif length(check0) == 2
             element = element*sqrt(3)
         end
      end

    end







    for pp = 1:Np-1

        vecmbindnn3 .= vecmbindnn[findall(x->x!=pp,1:Np)]

        for qq = 1:Np-1

            vecmbindmm3 .= vecmbindmm[findall(x->x!=qq,1:Np)]

            if vecmbindnn3 == vecmbindmm3

               element = 1.0

               check0 = findall(x->x==vecmbindnn[pp],vecmbindnn3)
               if length(check0) == 0
                  element = element*1.0
               elseif length(check0) == 1
                  element = element*sqrt(2)
               elseif length(check0) == 2
                  element = element*sqrt(3)
               end

               check1 = findall(x->x==vecmbindmm[qq],vecmbindmm3)
               if length(check1) == 0
                  element = element*1.0
               elseif length(check1) == 1
                  element = element*sqrt(2)
               elseif length(check1) == 2
                  element = element*sqrt(3)
               end

               return element

            end

        end
    end

    return 0.0

end
