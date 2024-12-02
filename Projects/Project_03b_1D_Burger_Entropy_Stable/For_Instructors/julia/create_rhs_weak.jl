#---------------------------------------------------------------------#
#This function computes the RHS vector R from the Global Matrices.
#Written by F.X. Giraldo on April 19 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function create_rhs_weak(q,M,Dwe,intma,periodicity,DFloat)

    #Get lengths of arrays
    (Np,nelem)=size(intma)
    Npoin=size(q)

    #Store local arrays
    inode=zeros(Int64,Np)
    rhs=zeros(DFloat,Npoin)

    #store flux term
    f=0.5*q.^2

    #build differentiation matrix contribution
    fe=zeros(DFloat,Np)
    for e=1:nelem
        for i=1:Np
            I=intma[i,e]
            fe[i]=f[I]
        end
        for i=1:Np
            I=intma[i,e]
            for j=1:Np
                rhs[periodicity[I]]+=Dwe[i,j]*fe[j]
            end
        end
    end

    #build flux matrix contribution (numerical flux)
    for e=1:nelem
        L=intma[Np,e]
        ep1=e+1
        if (e == nelem)
            ep1=1
        end
        R=intma[1,ep1]
        umax=max(abs(q[L]),abs(q[R]))
        flux=0.5*( f[L] + f[R] - umax*(q[R] - q[L]) )
        rhs[periodicity[L]]-=flux
        rhs[periodicity[R]]+=flux
    end

    #Build RHS vector
    R=M\rhs

    #Return arrays
    return(R)
end
