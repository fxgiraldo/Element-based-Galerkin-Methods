%---------------------------------------------------------------------%
%This routine advances the solution in time using SSP Osher-Shu RK2 or RK3
%Written by F.X. Giraldo on 9/6/19
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qp] = ssp_rk(q0,c,kstages,dt,time)

qp=q0;
q1=q0;
for ik=1:kstages
  switch kstages
      case 2  %RK2
          switch ik
              case 1
                 a0=1;
                 a1=0;
                 beta=1;
                 t=time;
              case (2)
                 a0=0.5;
                 a1=0.5;
                 beta=0.5;
                 t=time+dt;
          end %ik
      case 3 %RK3
          switch ik
              case 1
                 a0=1;
                 a1=0;
                 beta=1;
                 t=time;
              case (2)
                 a0=3.0/4.0;
                 a1=1.0/4.0;
                 beta=1.0/4.0;
                 t=time+dt;
              case (3)
                 a0=1.0/3.0;
                 a1=2.0/3.0;
                 beta=2.0/3.0;
                 t=time+0.5*dt;
          end %ik
  end %switch kstages

  %Construct RHS
  [rhs,f,g] = rhs_function(qp,c,t);
  
  %Update
  qp=a0*q0 + a1*q1 + dt*beta*rhs;
  q1=qp;
end %for ik
   
   
