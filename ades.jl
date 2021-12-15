function ades(jj::Int,fock::Vector{Int})

    # fock is a fock state, a vector having integers
    # its last element is square of coefficient

    focksize = length(fock) 

    fock[jj] = fock[jj]-1
    fock[focksize] = fock[focksize]*(fock[jj]+1)

    if fock[jj] < 0
       fock[focksize] = 0
    end

    return fock

end
