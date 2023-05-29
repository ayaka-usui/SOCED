include("in2bind.jl")

function perms(a::Vector{Int64})

    return unique(collect(permutations(a))) # list of all permutations

end

function changefrom2ndto1st_downdown_upup(indvec::Vector{Int64}, Msize0::Int64, Np::Int64, matp::Matrix{Int64})

    maxmatpcut = length(indvec)
    vecmbindnn = zeros(Int64,Np)
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

        # index for down (up) is less than or equal to M, and index for up (down) is more than or equal to M+1. 
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
    vec_space_jj = zeros(Int64,3)

    mat_from2ndto1st = spzeros(Float64,maxmatpcut2_new,maxmatpcut2) #zeros(Float64,maxmatpcut2_new,maxmatpcut2)

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
        mat_from2ndto1st[ind0+1:ind0+indvec2_new[nn],nn] .= 1/sqrt(indvec2_new[nn])
        ind0 += indvec2_new[nn]

    end

    return basis_1st, mat_from2ndto1st
    
end

function distinguish_downdownup!(nn::Int64,basis_1st::Matrix{Int64},index_spin::Matrix{Int64},ind::Int64)

    if basis_1st[nn,6] == 1 # down down up  
       index_spin[ind,1] = nn
    elseif basis_1st[nn,5] == 1 # down up down
       index_spin[ind,2] = nn
    else # basis_1st[nn,4] == 1 # up down down
       index_spin[ind,3] = nn
    end
    
end

function check_nn_included(nn::Int64,index_spin_checked::Matrix{Int64})

    if length(index_spin_checked) == 0
       check = findfirst(x->x==1,0)
    else
       check = findfirst(x->x==nn, index_spin_checked)
    end
    
    return check

end

function changefrom1sttospin_downup(basis_1st::Matrix{Int64})

    basis_1st_edit = copy(basis_1st)
    maxmatpcut2_new = size(basis_1st)[1]
    mat_from1sttospin_downdownup = spzeros(Float64,maxmatpcut2_new,maxmatpcut2_new) #zeros(Float64,maxmatpcut2_new,maxmatpcut2_new)
    mat_from1sttospin_upupdown = spzeros(Float64,maxmatpcut2_new,maxmatpcut2_new) #zeros(Float64,maxmatpcut2_new,maxmatpcut2_new)
    index_spin = zeros(Int64,Int64(maxmatpcut2_new/3),3)
    ind = 0
    ind_avoid = zeros(Int64,maxmatpcut2_new)

    mat_S2 = spzeros(Int64,maxmatpcut2_new,maxmatpcut2_new) #zeros(Int64,maxmatpcut2_new,maxmatpcut2_new)
    mat_M2 = spzeros(Int64,maxmatpcut2_new,maxmatpcut2_new)
    mat_M4 = spzeros(Int64,maxmatpcut2_new,maxmatpcut2_new)

    mat_S3 = spzeros(Int64,maxmatpcut2_new,maxmatpcut2_new)
    mat_M1 = spzeros(Int64,maxmatpcut2_new,maxmatpcut2_new)
    mat_M3 = spzeros(Int64,maxmatpcut2_new,maxmatpcut2_new)

    for nn = 1:maxmatpcut2_new
        
        if ind_avoid[nn] == 1
           continue
        end
        # check = check_nn_included(nn,index_spin[1:ind,:])
        # # check = findfirst(x -> x==nn, index_spin[1:ind,1:3])
        # if !isa(check,Nothing) # if check is not Nothing
        #    continue
        # end

        ind += 1
        distinguish_downdownup!(nn,basis_1st,index_spin,ind)
        ind_avoid[nn] = 1
        count = 0

        for jj = 1:maxmatpcut2_new

            if jj == nn
               continue
            end

            if ind_avoid[jj] == 1
               continue
            end
            # check = check_nn_included(jj,index_spin[1:ind,:])
            # # check = findfirst(x -> x==jj, index_spin[1:ind,1:3])
            # if !isa(check,Nothing)
            #    continue
            # end

            if basis_1st[nn,1] == basis_1st[jj,1] && basis_1st[nn,2] == basis_1st[jj,2] && basis_1st[nn,3] == basis_1st[jj,3]
               distinguish_downdownup!(jj,basis_1st,index_spin,ind)
               ind_avoid[jj] = 1
               count += 1
            end
            if count == 2
               break
            end

        end

        # println("ind=",ind)
        # println("It ends at ind=",Int64(maxmatpcut2_new/3))

        # println(index_spin[ind,1:3])

        # down down up -> S2 M2 M4 
        mat_from1sttospin_downdownup[index_spin[ind,1],index_spin[ind,1]] = 1/sqrt(3)
        mat_from1sttospin_downdownup[index_spin[ind,1],index_spin[ind,2]] = 1/sqrt(3)
        mat_from1sttospin_downdownup[index_spin[ind,1],index_spin[ind,3]] = 1/sqrt(3)
        mat_from1sttospin_downdownup[index_spin[ind,2],index_spin[ind,1]] =-2/sqrt(6)
        mat_from1sttospin_downdownup[index_spin[ind,2],index_spin[ind,2]] = 1/sqrt(6)
        mat_from1sttospin_downdownup[index_spin[ind,2],index_spin[ind,3]] = 1/sqrt(6)
        mat_from1sttospin_downdownup[index_spin[ind,3],index_spin[ind,1]] = 0
        mat_from1sttospin_downdownup[index_spin[ind,3],index_spin[ind,2]] =-1/sqrt(2)
        mat_from1sttospin_downdownup[index_spin[ind,3],index_spin[ind,3]] = 1/sqrt(2)
        mat_S2[index_spin[ind,1],index_spin[ind,1]] = 1
        mat_M2[index_spin[ind,2],index_spin[ind,2]] = 1
        mat_M4[index_spin[ind,3],index_spin[ind,3]] = 1

        # up up down -> S3 M1 M3 
        mat_from1sttospin_upupdown[index_spin[ind,1],index_spin[ind,1]] = 1/sqrt(3)
        mat_from1sttospin_upupdown[index_spin[ind,1],index_spin[ind,2]] = 1/sqrt(3)
        mat_from1sttospin_upupdown[index_spin[ind,1],index_spin[ind,3]] = 1/sqrt(3)
        mat_from1sttospin_upupdown[index_spin[ind,2],index_spin[ind,1]] = 2/sqrt(6)
        mat_from1sttospin_upupdown[index_spin[ind,2],index_spin[ind,2]] =-1/sqrt(6)
        mat_from1sttospin_upupdown[index_spin[ind,2],index_spin[ind,3]] =-1/sqrt(6)
        mat_from1sttospin_upupdown[index_spin[ind,3],index_spin[ind,1]] = 0
        mat_from1sttospin_upupdown[index_spin[ind,3],index_spin[ind,2]] = 1/sqrt(2)
        mat_from1sttospin_upupdown[index_spin[ind,3],index_spin[ind,3]] =-1/sqrt(2)
        mat_S3[index_spin[ind,1],index_spin[ind,1]] = 1
        mat_M1[index_spin[ind,1],index_spin[ind,1]] = 1
        mat_M3[index_spin[ind,1],index_spin[ind,1]] = 1
        
    end

    # return index_spin, basis_1st

    # println(sum(ind_avoid))
    # println(maxmatpcut2_new)

    return mat_from1sttospin_downdownup, mat_from1sttospin_upupdown, mat_S2, mat_M2, mat_M4, mat_S3, mat_M1, mat_M3

end

function changefrom2ndtospin_downup(indvec2::Vector{Int64}, Msize0::Int64, Np::Int64, matp20::Matrix{Int64}, matp21::Matrix{Int64})

    basis_1st, mat_from2ndto1st_downup = changefrom2ndto1st_downup(indvec2,Msize0,Np,matp20,matp21)
    mat_from1sttospin_downdownup, mat_from1sttospin_upupdown, mat_S2, mat_M2, mat_M4, mat_S3, mat_M1, mat_M3 = changefrom1sttospin_downup(basis_1st)
    # index_spin, basis_1st = changefrom1sttospin_downup(basis_1st)
    
    return mat_from2ndto1st_downup, mat_from1sttospin_downdownup, mat_from1sttospin_upupdown, mat_S2, mat_M2, mat_M4, mat_S3, mat_M1, mat_M3
    # return index_spin, basis_1st

end










