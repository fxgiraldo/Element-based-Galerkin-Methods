#---------------------------------------------------------------------#
#This function computes the RK Coefficients.
#Written by F.X. Giraldo on April 19, 2019
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey; CA 93943-5216
#---------------------------------------------------------------------#
include("create_rhs.jl")

function ti_SSP!(q0,u,Dhat,Fhat,intma,periodicity,time,ntime,dt,stages,DFloat)

    #Initialize RK coefficients
    a0=zeros(DFloat,stages)
    a1=zeros(DFloat,stages)
    beta=zeros(DFloat,stages)
    if stages == 1 #RK1-SSP
        a0[1]=1; a1[1]=0; beta[1]=1
    elseif stages == 2 #RK2-SSP
        a0[1]=1; a1[1]=0; beta[1]=1; #SSP
        a0[2]=1.0/2.0; a1[2]=1.0/2.0; beta[2]=1.0/2.0
    elseif stages == 3 #RK3-SSP
        a0[1]=1; a1[1]=0; beta[1]=1;
        a0[2]=3.0/4.0; a1[2]=1.0/4.0; beta[2]=1.0/4.0
        a0[3]=1.0/3.0; a1[3]=2.0/3.0; beta[3]=2.0/3.0
    end

    #Initialize arrays
    q1=copy(q0); qp=copy(q0)
    Npoin=length(q0)

    #Time Integration
    for itime=1:ntime
        time=time + dt
        #RK Stages
        for s=1:stages
            #Create RHS Matrix
#            R=-Dhat*qp
            R=create_rhs(qp,u,Dhat,Fhat,intma,DFloat)

            #Solve System
            @inbounds for I=1:Npoin
                qp[I]=a0[s]*q0[I] + a1[s]*q1[I] + dt*beta[s]*R[I]
            end #s
            if periodicity[Npoin] == periodicity[1]
                qp[Npoin]=qp[1] #periodicity
            end
            #Update
            q1 .= qp
        end #ik
#        println(" itime = ",itime," time = ",time," qmax = ",maximum(qp)," qmin = ",minimum(qp))

        #Update Q
        q0 .= qp
    end #itime

    return (q0)
end
