function Vijkl(ii::Int64,jj::Int64,kk::Int64,ll::Int64)

    n1 = ii
    n2 = jj
    n3 = kk
    n4 = ll

    A1234 = 1

    for m1 = 0:ceil(Int64,n1/2)
        for m2 = 0:ceil(Int64,n2/2)
            for m3 = 0:ceil(Int64,n3/2)
                for m4 = 0:ceil(Int64,n4/2)

                    if m1==m2==m3==m4
                       continue
                    end

                    A1234 = A1234*(n4-2*m4)*(n4-2*m4-1)/(m4+1)*2*((n1+n2+n3+n4)/2-m1-m2-m3-m4)/(n1+n2+n3+n4-2m1-2m2-2m3-2m4)/(n1+n2+n3+n4-2m1-2m2-2m3-2m4-1)

                end
            end
        end
    end

end
