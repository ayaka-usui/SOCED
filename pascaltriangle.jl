function pascaltriangle(Msize::Int,Np::Int)

    Lm = Msize+1
    Ln = Np + 1

    matp = zeros(Int,Lm,Ln)
    matp[:,1] .= 1

    for jj = 2:Lm
        for ii = 2:Ln
            matp[jj,ii] = sum(matp[jj-1,1:ii])
        end
    end

    return matp

end
