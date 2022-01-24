function ades!(jj::Int64, input::Union{SparseVector{Int64},Vector{Int64}}, fock::Union{SparseVector{Int64},Array{Int64}}, tid::Int64)

    # anihilation operator, a_{jj}

    # fock is a fock state, a vector having integers
    # its last element is square of coefficient
    fock[:,tid] .= input
    focksize = size(fock)[1]

    fock[jj,tid] = fock[jj]-1
    fock[focksize,tid] = input[focksize]*(fock[jj,tid]+1)

    if fock[jj,tid] < 0
       fock[focksize,tid] = 0
    end

    # return fock

end
