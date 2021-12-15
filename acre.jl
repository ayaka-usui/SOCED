function acre(jj::Int,fock::Vector{Int})

    # fock is a fock state, a vector having integers
    # its last element is square of coefficient

    focksize = length(fock) 

    fock[jj] = fock[jj]+1
    fock[focksize] = fock[focksize]*fock[jj]

    return fock

end
