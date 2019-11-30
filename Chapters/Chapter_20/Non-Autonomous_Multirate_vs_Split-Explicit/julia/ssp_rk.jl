#---------------------------------------------------------------------#
#This routine advances the solution in time using SSP Osher-Shu RK2 or RK3
#Written by F.X. Giraldo / P.R. Mugg on 9/20/19
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey, CA 93943-5216
#---------------------------------------------------------------------#
function ssp_rk(q0,c,kstages,Δt,time)

   qp = q0
   q1 = q0
   for ik = 1:kstages
      if (kstages == 2)     #RK2
         if (ik == 1)
            α0 = 1.0; α1 = 0.0; β = 1.0
            t = time
         elseif (ik == 2)
            α0 = 0.5; α1 = 0.5; β = 0.5
            t = time + Δt
         end # ik
      elseif (kstages == 3)  #RK3
         if (ik == 1)
            α0 = 1.0; α1 = 0.0; β = 1.0
            t = time
         elseif (ik == 2)
            α0 = 3.0/4.0; α1 = 1.0/4.0; β = 1.0/4.0
            t = time + Δt
         elseif (ik == 3)
            α0 = 1.0/3.0; α1 = 2.0/3.0; β = 2.0/3.0
            t = time + (Δt/2)
         end # ik
      end #kstages

      #Construct RHS
      (rhs,f,g) = rhs_function(qp,c,t)

      #Update
      qp = α0*q0 + α1*q1 + Δt*β*rhs
      q1 = qp
   end #for ik

   return qp
end #end function
