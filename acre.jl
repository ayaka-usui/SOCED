function acre!(jj::Int64, input::Union{SparseVector{Int64},Vector{Int64}}, fock::Union{SparseVector{Int64},Array{Int64}}, tid)

    # creation operator, a_{jj}^{\dagger}

    # fock is a fock state, a vector having integers
    # its last element is square of coefficient
    fock[:,tid] .= input
    focksize = size(fock)[1]

    fock[jj,tid] = fock[jj,tid]+1
    fock[focksize,tid] = fock[focksize,tid]*fock[jj,tid]

    # return fock

end
