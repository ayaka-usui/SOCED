function acre_spinless!(jj::Int64, input::Union{SparseVector{Int64},Vector{Int64}}, fock::Union{SparseVector{Int64},Array{Int64}})

    # creation operator, a_{jj}^{\dagger}

    # fock is a fock state, a vector having integers
    # its last element is square of coefficient
    fock[:] .= input
    focksize = size(fock)[1]

    fock[jj] = fock[jj]+1
    fock[focksize] = fock[focksize]*fock[jj]

    # return fock

end
