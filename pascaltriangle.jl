function pascaltriangle!(Msize0::Int64, Np::Int64, matp::Matrix{Int64})

    Lm = Msize0+1
    Ln = Np + 1

    matp .= 0
    matp[:,1] .= 1

    for jj = 2:Lm
        for ii = 2:Ln
            matp[jj,ii] = sum(matp[jj-1,1:ii])
        end
    end

end
