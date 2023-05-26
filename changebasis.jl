include("in2bind.jl")

function perms(a::Vector{Int64})

    return unique(collect(permutations(a))) # list of all permutations

end

function changefrom2ndto1st_downdown_upup(indvec::Vector{Int64}, Msize0::Int64, Np::Int64, matp::Matrix{Int64})

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

function list_from_index_downup!(indvec2_nn::Int64, Msize0::Int64, Np::Int64, matp20::Matrix{Int64}, matp21::Matrix{Int64}, maxmatp21::Int64, vecmbindnn::Vector{Int64}, vecmbindnn2::Vector{Int64})

    # indvec2[nn] = (nndown-1)*maxmatp21 + nnup
    indketup = mod(indvec2_nn,maxmatp21)
    if indketup != 0
       indketdown = div(indvec2_nn,maxmatp21) + 1
    else # indketup == 0
       indketdown = div(indvec2_nn,maxmatp21)
       indketup = maxmatp21
    end
    in2bind!(indketdown,Msize0,Np-1,matp20,vecmbindnn2)
    vecmbindnn[1:Np-1] = vecmbindnn2[1:Np-1]
    in2bind!(indketup,Msize0,1,matp21,vecmbindnn2)
    vecmbindnn[Np:Np] = vecmbindnn2[1:1]

end

function changefrom2ndto1st_downup(indvec2::Vector{Int64}, Msize0::Int64, Np::Int64, matp20::Matrix{Int64}, matp21::Matrix{Int64})
    
    maxmatpcut2 = length(indvec2)

    vecmbindnn = zeros(Int64,Np)
    vecmbindnn2 = zeros(Int64,Np)
    maxmatp21 = matp21[Msize0+1,1+1]

    indvec2_new = zeros(Int64,maxmatpcut2)

    # compute size of basis_1st
    for nn = 1:maxmatpcut2

        list_from_index_downup!(indvec2[nn],Msize0,Np,matp20,matp21,maxmatp21,vecmbindnn,vecmbindnn2)

        # index for down is less than or equal to M, and index for up is more than or equal to M+1. 
        vecmbindnn[3] += Msize0

        # println(vecmbindnn)
        
        array_ind = length(unique(vecmbindnn))

        if array_ind == 3
           indvec2_new[nn] = 6
        elseif array_ind == 2
           indvec2_new[nn] = 3
        else # array_ind == 1
           indvec2_new[nn] = 1
        end

    end

    maxmatpcut2_new = sum(indvec2_new)
    basis_1st = zeros(Int64,maxmatpcut2_new,6)

    ind0 = 0
    # ind1 = ind0 + indvec2_new[nn] - 1
    # ind0 = ind1 + 1
    vec_space_jj = zeros(Int64,3)

    for nn = 1:maxmatpcut2

        list_from_index_downup!(indvec2[nn],Msize0,Np,matp20,matp21,maxmatp21,vecmbindnn,vecmbindnn2)

        # index for down is less than or equal to M, and index for up is more than or equal to M+1. 
        vecmbindnn[3] += Msize0    

        vec_space = perms(vecmbindnn)
        # println(vec_space)

        for jj = 1:indvec2_new[nn]

            vec_space_jj .= vec_space[jj]

            indup = findfirst(x -> x>Msize0, vec_space_jj)
            vec_space_jj[indup] = vec_space_jj[indup] - Msize0
            basis_1st[ind0+jj,1:3] .= vec_space_jj
            # basis_1st[ind0+jj,4:6] .= 0
            basis_1st[ind0+jj,Np+indup] = 1

        end
        ind0 += indvec2_new[nn]

    end

    return basis_1st

    # # maxmatpcut2_new = sum(indvec2_new)
    # mat_from2ndto1st = zeros(Float64,maxmatpcut2_new,maxmatpcut2) #spzeros(Float64,maxmatpcut2_new,maxmatpcut2)
    # ind0 = 1

    # for nn = 1:maxmatpcut2

    #     ind1 = ind0 + indvec2_new[nn] - 1
    #     mat_from2ndto1st[ind0:ind1,nn] .= 1/sqrt(indvec2_new[nn]) 
    #     ind0 = ind1 + 1

    # end

    # return mat_from2ndto1st
    
end














