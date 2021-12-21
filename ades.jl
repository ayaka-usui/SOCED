function ades!(jj::Int64, input::SparseVector{Int64}, fock::SparseVector{Int64})

    # fock is a fock state, a vector having integers
    # its last element is square of coefficient
    fock .= input
    focksize = length(fock)

    fock[jj] = fock[jj]-1
    fock[focksize] = input[focksize]*(fock[jj]+1)

    if fock[jj] < 0
       fock[focksize] = 0
    end

    # return fock

end
