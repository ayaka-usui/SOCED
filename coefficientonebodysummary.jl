
function coefficientonebody(vecmbindnn::Vector{Int64},vecmbindmm::Vector{Int64},vecmbindnn3::Vector{Int64},vecmbindmm3::Vector{Int64},Np::Int64)

    # vecmbindnn3 = zeros(Int64,2)
    # vecmbindmm3 = zeros(Int64,2)

    # a|n>= sqrt(n)|n-1>
    # a^+|n>= sqrt(n+1)|n+1>

    # initilise
    vecmbindnn3 .= 0
    vecmbindmm3 .= 0

    check0 = 0

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
