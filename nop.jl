function nop(jj::Int,fock::Vector{Int})

    # fock is a fock state, a vector having integers
    # its last element is square of coefficient

    focksize = length(fock)

    fock[focksize] = fock[focksize]*fock[jj]^2

    return fock

end
