include("delta.jl")

function epsilon2(ii::Int64, jj::Int64)

    # snigle-particle energy, epsilon

    energy0 = 0 # Int64
    # njj = ceil(Int64,jj/2) - 1
    # nii = ceil(Int64,ii/2) - 1

    if ii == jj

       # energy0 = njj + 1/2 # Float64
       energy0 = (ii-1) + 1/2 # Float64

    end

    return energy0

end
