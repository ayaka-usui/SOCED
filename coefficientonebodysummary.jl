
function coefficientonebody(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},vecmbindnn3::Vector{Int64},vecmbindmm3::Vector{Int64},vecindcoeff::Matrix{Float64},Np::Int64)

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

    ind0 = 0
    ind1 = 0
    ind2 = 0
    ind3 = 0

    for pp = 1:Np

        for qq = 1:Np

            if vecmbindnn[pp] == vecmbindmm[qq]
               ind0 += 1
               continue
            else
               ind1 = pp
               ind2 = qq
            end

        end
    end

    if ind0 != 5
       return 0.0
    end






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
