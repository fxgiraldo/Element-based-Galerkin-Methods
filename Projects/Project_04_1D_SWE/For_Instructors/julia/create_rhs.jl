#---------------------------------------------------------------------#
#This function computes the global Mass and Differentiation Matrices.
#Written by F.X. Giraldo on April 19 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
function create_rhs(q,qb,M,Dwe,intma,periodicity,gravity,Δ_NL,DFloat)

    #Get lengths of arrays
    (Np,nelem)=size(intma)
    Npoin=size(q,1)

    #Store local arrays
    rhs=zeros(DFloat,Npoin,2)

    #build differentiation matrix contribution
    fe=zeros(DFloat,Np,2)
    for e=1:nelem
        for i=1:Np
            I=intma[i,e]
            hs=q[I,1]
            hb=qb[I,1]
            h=hs+hb
            u=q[I,2]/h
            fe[i,1]=h*u
            fe[i,2]=(h*u*u + 0.5*(gravity*hs^2))*Δ_NL + gravity*hs*hb
        end
        for i=1:Np
            I=intma[i,e]
            for j=1:Np
                rhs[periodicity[I],:]+=Dwe[i,j]*fe[j,:]
            end
        end
    end

    #build flux matrix contribution (numerical flux)
    flux=zeros(DFloat,2)
    feL=zeros(DFloat,2)
    feR=zeros(DFloat,2)
    for e=1:nelem
        #Left State
        eL=e
        L=intma[Np,eL]
        hsL=q[L,1]
        hbL=qb[L]
        hL=hsL + hbL
        UL=q[L,2]
        uL=UL/hL
        feL[1]=hL*uL
        feL[2]=( hL*uL*uL + 0.5*gravity*hsL^2 )*Δ_NL + gravity*hsL*hbL

        #Right State
        eR=e+1
        if (e == nelem)
            eR=1
        end
        R=intma[1,eR]
        hsR=q[R,1]
        hbR=qb[R]
        hR=hsR + hbR
        UR=q[R,2]
        uR=UR/hR
        feR[1]=hR*uR
        feR[2]=( hR*uR*uR + 0.5*gravity*hsR^2 )*Δ_NL + gravity*hsR*hbR

        #Jump Speed
        λL=( abs(uL) + sqrt(gravity*hL) )*Δ_NL + sqrt(gravity*hbL)*(1.0-Δ_NL);
        λR=( abs(uR) + sqrt(gravity*hR) )*Δ_NL + sqrt(gravity*hbR)*(1.0-Δ_NL);
        λ=max( λL, λR )
        flux[1]=0.5*( feL[1] + feR[1] - λ*(hR - hL) )
        flux[2]=0.5*( feL[2] + feR[2] - λ*(UR - UL) )
        rhs[periodicity[L],:]-=flux[:]
        rhs[periodicity[R],:]+=flux[:]
    end

    #Build RHS vector
    rhs[:,1]=M\rhs[:,1]
    rhs[:,2]=M\rhs[:,2]

    #Return arrays
    return(rhs)
end
