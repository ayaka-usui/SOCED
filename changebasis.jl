include("in2bind.jl")

function changefrom2ndto1st(indvec::Vector{Int64}, Msize0::Int64, Np::Int64, matp::Matrix{Int64})

    maxmatpcut = length(indvec)
    vecmbindnn = zeros(Int64,Np)
    # perms(a) = unique(collect(permutations(a))) # list of all permutations
    indvec_new = zeros(Int64,maxmatpcut)
    
    for nn = 1:maxmatpcut

        # 2nd quantisation label
        in2bind!(indvec[nn],Msize0,Np,matp,vecmbindnn)
        array_ind = length(unique(vecmbindnn))

        if array_ind == 3 #length(vecmbindnn)
           indvec_new[nn] = 6 #3!
        elseif array_ind == 2
           indvec_new[nn] = 3
        else # array_ind == 1
           indvec_new[nn] = 1
        end

        println(vecmbindnn)

    end

    maxmatpcut_new = sum(indvec_new)
    mat_from2ndto1st = zeros(Float64,maxmatpcut_new,maxmatpcut)
    ind0 = 1

    for nn = 1:maxmatpcut

        ind1 = ind0 + indvec_new[nn] - 1
        mat_from2ndto1st[ind0:ind1,nn] .= 1/sqrt(indvec_new[nn]) 
        ind0 = ind1 + 1

    end

    return mat_from2ndto1st

end

function ()
    
end








