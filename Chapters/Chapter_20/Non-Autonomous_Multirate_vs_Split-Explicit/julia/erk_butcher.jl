#---------------------------------------------------------------------#
#This routine advances the solution in time using Explicit RK methods
#written in Butcher tableau form.
#Written by F.X. Giraldo / P.R. Mugg on 9/22/19
#           Department of Applied Mathematics
#           Naval Postgraduate School
#           Monterey, CA 93943-5216
#---------------------------------------------------------------------#
function erk_butcher(q0,c,αI,βI,I,Δt,time)

   Q = zeros(I,1)
   R = zeros(I,1)

   c𝑖 = zeros(I)
   for i = 1:I
      c𝑖[i] = sum(αI[i,:])
   end

   Q[1] = q0
   t = time
   (rhs,f,g) = rhs_function(Q[1],c,t)
   R[1] = rhs

   for i = 2:I
      R_sum = 0
      for j = 1:i-1
         R_sum += αI[i,j]*R[j]
      end
      Q[i] = q0 + Δt*R_sum
      t = time + Δt*c𝑖[i]
      (rhs,f,g) = rhs_function(Q[i],c,t)
      R[i] = rhs
   end

   R_sum = 0
   for i = 1:I
      R_sum += βI[i]*R[i]
   end

   return qp = q0 + Δt*R_sum
end #end function
