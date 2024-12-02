#---------------------------------------------------------------------#
#This function computes the RHS vector R from the Element Matrices.
#It can do either standard fluxes (NES) or Entropy-Stable (ES).
#Written by F.X. Giraldo on December 2, 2024.
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function create_rhs_strong(q,M,De,intma,periodicity,flux_method,DFloat)

    #Get lengths of arrays
    (Np,nelem)=size(intma)
    Npoin=size(q)

    #Store local arrays
    rhs=zeros(DFloat,Npoin)
    qe=zeros(DFloat,Np)
    fe=zeros(DFloat,Np,Np)

    #build differentiation matrix contribution
    for e=1:nelem
        #store local solution
        for i=1:Np
            I=intma[i,e]
            qe[i]=q[I]
        end
        #build flux matrix
        if (flux_method == "ES")
            for i=1:Np
                for j=1:Np
                    fe[i,j]=( qe[i]*qe[i] + qe[i]*qe[j] + qe[j]*qe[j] )/6.0
                end
            end
        elseif (flux_method == "NES")
            for i=1:Np
                for j=1:Np
                    fe[i,j]=( 0.5*qe[i]*qe[i] + 0.5*qe[j]*qe[j] )/2.0
                end
            end
        end

        #Augment RHS vector with Volume Flux terms
        for i=1:Np
            I=intma[i,e]
            for j=1:Np
                rhs[periodicity[I]]-=2.0*De[i,j]*fe[i,j]
            end
        end
    end

    #build flux matrix contribution (numerical flux)
    for e=0:nelem
        L=e
        R=e+1

        #Store Left State
        if (L>=1)
            IL=intma[Np,L]
            q_L=q[IL]
        else
            IL=0
            q_L=0 #homogeneous Dirichlet BC
        end
        
        #Store Right State
        if (R<=nelem)
            IR=intma[1,R]
            q_R=q[IR]
        else
            IR=0
            q_R=0 #homogeneous Dirichlet BC
        end
    
        #Build flux*
        umax=max(abs(q_L),abs(q_R))
        f_L=0.5*q_L^2
        f_R=0.5*q_R^2
        if flux_method == "ES"
            flux=( q_L^2 + q_L*q_R + q_R^2)/6.0
        elseif flux_method == "NES"
            flux=( 0.5*q_L^2 + 0.5*q_R^2)/2.0
        end
        flux_star=flux - 0.5*umax*(q_R - q_L) 

        #Augment RHS vector with Surface Flux terms
        if (IL > 0)
            rhs[periodicity[IL]]-=(flux_star-f_L)
        end
        if (IR > 0)
            rhs[periodicity[IR]]+=(flux_star-f_R)
        end
    end

    #Build RHS vector
    R=M\rhs

    #Return arrays
    return(R)
end
