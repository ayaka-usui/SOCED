
function coefficientonebody1!(jj::Int64, psi::Vector{Float64}, Msize0::Int64, Np::Int64, maxmatpcut::Int64, maxmatpcut2::Int64, indvec::Vector{Int64}, indvec2::Vector{Int64}, matp::Matrix{Int64}, matp20::Matrix{Int64}, matp21::Matrix{Int64})

   psi .= 0.0

   vecmbindnn = zeros(Int64,Np)
   maxmatp21 = matp21[Msize0+1,1+1]
   vecmbindnn2 = zeros(Int64,Np)

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

   for nn = 1:maxmatpcut2

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

       # down down up
       check0 = findall(x->x==jj,vecmbindnn[1:Np-1])
       if length(check0) == 1
          psi[maxmatpcut+nn] = 1.0
       elseif length(check0) == 2
          psi[maxmatpcut+nn] = sqrt(2)
       # else # length(check0) == 0
          # psi[maxmatpcut+nn] = 0.0
       end

       # up up down
       if vecmbindnn[Np] == jj
          psi[maxmatpcut+maxmatpcut2+nn] = 1.0
       # else # vecmbindnn[Np:Np] != jj
          # psi[maxmatpcut+maxmatpcut2+nn] = 0.0
       end

   end

   # for nn = 1:maxmatpcut2 # up up down
   #
   #     # ket
   #     indketup = mod(indvec2[nn],maxmatp21)
   #     if indketup != 0
   #        indketdown = div(indvec2[nn],maxmatp21) + 1
   #     else # indketup == 0
   #        indketdown = div(indvec2[nn],maxmatp21)
   #        indketup = maxmatp21
   #     end
   #     in2bind!(indketdown,Msize0,Np-1,matp20,vecmbindnn2)
   #     vecmbindnn[1:Np-1] = vecmbindnn2[1:Np-1]
   #     in2bind!(indketup,Msize0,1,matp21,vecmbindnn2)
   #     vecmbindnn[Np:Np] = vecmbindnn2[1:1]
   #
   #     if vecmbindnn[Np:Np] == jj
   #        psi[maxmatpcut+maxmatpcut2+nn] = 1.0
   #     # else # vecmbindnn[Np:Np] != jj
   #        # psi[maxmatpcut+maxmatpcut2+nn] = 0.0
   #     end
   #
   # end

   # psi[maxmatpcut+maxmatpcut2++maxmatpcut2+1:end] = 0.0 # up up up

end

function coefficientonebody2!(jj::Int64, psi::Vector{Float64}, Msize0::Int64, Np::Int64, maxmatpcut::Int64, maxmatpcut2::Int64, indvec::Vector{Int64}, indvec2::Vector{Int64}, matp::Matrix{Int64}, matp20::Matrix{Int64}, matp21::Matrix{Int64})

   psi .= 0.0

   vecmbindnn = zeros(Int64,Np)
   maxmatp21 = matp21[Msize0+1,1+1]
   vecmbindnn2 = zeros(Int64,Np)

   # psi[1:maxmatpcut] = 0.0 # down down down

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

       if vecmbindnn[Np] == jj
          psi[maxmatpcut+nn] = 1.0
       # else # vecmbindnn[Np:Np] != jj
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

       check0 = findall(x->x==jj,vecmbindnn[1:Np-1])
       if length(check0) == 1
          psi[maxmatpcut+maxmatpcut2+nn] = 1.0
       elseif length(check0) == 2
          psi[maxmatpcut+maxmatpcut2+nn] = sqrt(2)
       # else # length(check0) == 0
          # psi[maxmatpcut+nn] = 0.0
       end

   end

   for nn = 1:maxmatpcut # up up up

       # ket
       in2bind!(indvec[nn],Msize0,Np,matp,vecmbindnn)

       check0 = findall(x->x==jj,vecmbindnn)
       if length(check0) == 1
          psi[maxmatpcut+maxmatpcut2+maxmatpcut2+nn] = 1.0
       elseif length(check0) == 2
          psi[maxmatpcut+maxmatpcut2+maxmatpcut2+nn] = sqrt(2)
       elseif length(check0) == 3
          psi[maxmatpcut+maxmatpcut2+maxmatpcut2+nn] = sqrt(3)
       # else # length(check0) == 0
          # psi[nn] = 0.0
       end

   end

end

function coefficientonebody!(jj::Int64, psijj::Vector{Float64}, Msize0::Int64, Np::Int64, maxmatpcut::Int64, maxmatpcut2::Int64, indvec::Vector{Int64}, indvec2::Vector{Int64}, matp::Matrix{Int64}, matp20::Matrix{Int64}, matp21::Matrix{Int64})

   psijj .= 0.0

   if jj <= Msize0
      coefficientonebody1!(jj,psijj,Msize0,Np,maxmatpcut,maxmatpcut2,indvec,indvec2,matp,matp20,matp21)
   else # jj >= Msize0+1
      coefficientonebody2!(jj-Msize0,psijj,Msize0,Np,maxmatpcut,maxmatpcut2,indvec,indvec2,matp,matp20,matp21)
   end

end
