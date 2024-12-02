#=
---------------------------------------------------------------------
This function advances the solution in time using the LSRK45 Method by Carpenter-Kennedy 1994.

Written by F.X. Giraldo on July 12, 2021
           Department of Applied Mathematics
           Naval Postgraduate School
           Monterey; CA 93943-5216
---------------------------------------------------------------------
=#

include("compute_courant.jl")
include("create_rhs.jl")
include("plot_solution.jl")

function ti_LSRK!(q0,coord,M,De,Dwe,intma,periodicity,time,time_final,Courant_max,space_method,cgdg_method,flux_method,plot_movie,DFloat)

    #Initialize RK coefficients
    RKA = (DFloat(0),
           DFloat(-567301805773)  / DFloat(1357537059087),
           DFloat(-2404267990393) / DFloat(2016746695238),
           DFloat(-3550918686646) / DFloat(2091501179385),
           DFloat(-1275806237668) / DFloat(842570457699 ))

    RKB = (DFloat(1432997174477) / DFloat(9575080441755 ),
           DFloat(5161836677717) / DFloat(13612068292357),
           DFloat(1720146321549) / DFloat(2090206949498 ),
           DFloat(3134564353537) / DFloat(4481467310338 ),
           DFloat(2277821191437) / DFloat(14882151754819))

    RKC = (DFloat(0),
           DFloat(1432997174477) / DFloat(9575080441755),
           DFloat(2526269341429) / DFloat(6820363962896),
           DFloat(2006345519317) / DFloat(3224310063776),
           DFloat(2802321613138) / DFloat(2924317926251))

    Npoin=length(q0)
    dq=zeros(DFloat,Npoin)
    qp=copy(q0)
    stages=length(RKA)
    dx_min=coord[2]-coord[1]
    itime=0

    #Time Integration
    while time < time_final
        itime=itime + 1

        #Compute new dt and courant
        (dt,courant)=compute_courant(qp,Courant_max,dx_min,DFloat)       
        time=time+dt
        
        #Make sure we finish at Time_Final
        if time>time_final
            time=time - dt;
            dt=time_final - time;
            time=time + dt;
        end

        #RK Stages
        for s = 1:stages

            #Create RHS Matrix
            #= -----------Students add your CREATE_RHS routine here-----------=#
            #= -----------Students add your CREATE_RHS routine here-----------=#
            R=create_rhs(qp,M,De,Dwe,intma,periodicity,cgdg_method,flux_method,DFloat)
            #= -----------Students add your CREATE_RHS routine here-----------=#
            #= -----------Students add your CREATE_RHS routine here-----------=#

            #Solve System
            @inbounds for I=1:Npoin
                dq[I] = RKA[s]*dq[I] + dt*R[I]
                qp[I] = qp[I] + RKB[s]*dq[I]
            end
            if periodicity[Npoin] == periodicity[1]
                qp[Npoin]=qp[1] #periodicity
            end
        end #s

        #Update Q
        q0 .= qp

        #Plot Solution
        if (plot_movie)
            plot_solution(q0,coord,space_method,time)
        end

         #Print to Screen
         if ( mod(itime,10) == 0 )
             println(" itime = ",itime," time = ",time," Courant = ",courant," dt = ",dt," qmax = ",maximum(abs,q0))
end

    end #itime

    return(q0,time,itime)

end
