

function changefrom2ndto1st(indvec::Vector{Int64}, Msize0::Int64, Np::Int64, matp::Matrix{Int64})

    maxmatpcut = length(indvec)
    vecmbindnn = zeros(Int64,Np)
    perms(a) = unique(collect(permutations(a))) # list of all permutations

    for nn = 1:maxmatpcut

        # 2nd quantisation label
        in2bind!(indvec[nn],Msize0,Np,matp,vecmbindnn)

        

    end

end