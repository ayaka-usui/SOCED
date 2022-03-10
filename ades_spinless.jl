function ades_spinless!(jj::Int64, input::Union{SparseVector{Int64},Vector{Int64}}, fock::Union{SparseVector{Int64},Array{Int64}})

    # anihilation operator, a_{jj}

    # fock is a fock state, a vector having integers
    # its last element is square of coefficient
    fock[:] .= input
    focksize = size(fock)[1]

    fock[jj] = fock[jj]-1
    fock[focksize] = input[focksize]*(fock[jj]+1)

    if fock[jj] < 0
       fock[focksize] = 0
    end

    # return fock

end
