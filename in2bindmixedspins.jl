function in2bindmixedspins(Msize0::Int64,Np::Int64,nn::Int64,indvec2::Vector{Int64},maxmatp21::Int64,matp20::Matrix{Int64},matp21::Matrix{Int64},vecmbindnn2::Vector{Int64},vecmbindnn::Vector{Int64})
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
   return vecmbindnn
end
