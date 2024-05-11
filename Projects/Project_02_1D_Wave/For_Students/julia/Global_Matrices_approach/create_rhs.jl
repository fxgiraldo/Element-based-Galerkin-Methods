#---------------------------------------------------------------------#
#This function computes the global Mass and Differentiation Matrices.
#Written by F.X. Giraldo on April 19 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function create_rhs(q,u,Dhat,Fhat,intma,DFloat)

    #Get lengths of arrays
    (Np,nelem)=size(intma)
    Npoin=size(q)

    #store flux term
    f=q*u

    #build numerical flux: fstar
    fstar=zeros(DFloat,Npoin)
    for e=1:nelem
        L=intma[Np,e]
        ep1=e+1
        if (e == nelem)
            ep1=1
        end
        R=intma[1,ep1]
        flux=0.5*( f[L] + f[R] - u*(q[R] - q[L]) )
        fstar[L]=flux
        fstar[R]=flux
    end

    #Build RHS vector
    R=Dhat*f - Fhat*fstar

    #Return arrays
    return(R)
end
