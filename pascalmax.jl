function pascalmax(Msize::Int,Np::Int)

    include("pascaltriangle.jl")
    matp = pascaltriangle(Msize,Np)

    return matp[Msize+1,Np+1]

end
