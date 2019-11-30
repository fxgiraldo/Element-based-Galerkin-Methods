#---------------------------------------------------------------------#
#This routine advances the solution in time using Split-Explicit using
#a General Order RK method in Butcher tableau form following Schlegel et al. GMD 2012.
#Written by F.X. Giraldo / P.R. Mugg on 9/22/19
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey, CA 93943-5216
#---------------------------------------------------------------------#
function multirate_erk_butcher_substepping(q0,c,αI,βI,I,αJ,βJ,J,M,Δt,time)

    #Compute Coefficients
    #Outer Loop
    α̃I = zeros(I+1,I)
    c̃I = zeros(I+1,1)
    cI = zeros(I,1)
    for i = 2:I
        cI[i] = sum(αI[i,:])
        α̃I[i,:] = αI[i,:] - αI[i-1,:]
        c̃I[i] = cI[i] - cI[i-1]
    end
    for j = 1:I
        α̃I[I+1,j] = βI[j] - αI[I,j]
    end
    c̃I[I+1] = 1.0 - cI[I]

    #Inner Loop
    α̃J = zeros(J+1,J)
    c̃J = zeros(J+1,1)
    cJ = zeros(J,1)
    for i = 2:J
        cJ[i] = sum(αJ[i,:])
        α̃J[i,:] = αJ[i,:] - αJ[i-1,:]
        c̃J[i] = cJ[i] - cJ[i-1]
    end
    for j = 1:J
        α̃J[J+1,j] = βJ[j] - αJ[J,j]
    end
    c̃J[J+1] = 1.0 - cJ[J]

    #Initialize Arrays
    Q = zeros(I+1,1)
    V = zeros(I+1,J+1)

    #Outer Loop
    Q[1] = q0
    for i = 2:I+1
        r_i = 0
        for j = 1:i-1
            t = time + Δt*cI[j]
            (rhs,f𝑗,g𝑗) = rhs_function(Q[j],c,t)
            r_i += α̃I[i,j]*g𝑗
        end

        #Substepping Loop
        qm = Q[i-1]
        for m = 1:M
            V[i,1] = qm

            #Inner Loop
            for j = 2:J+1
                R_sum = 0
                for k = 1:j-1
                    t = time + Δt*( cI[i-1] + c̃I[i]*cJ[k]/M*(m-1) ) #not correct so am losing accuracy
                    (rhs,f𝑘,g𝑘) = rhs_function(V[i,k],c,t)
                    R_sum += α̃J[j,k]*( r_i + c̃I[i]*f𝑘 )
                end #k
                V[i,j] = V[i,j-1] + Δt/M*R_sum
            end #j
            qm = V[i,J+1]
        end #m
        Q[i] = qm
    end #i

    return qp = Q[I+1]
end #end function
