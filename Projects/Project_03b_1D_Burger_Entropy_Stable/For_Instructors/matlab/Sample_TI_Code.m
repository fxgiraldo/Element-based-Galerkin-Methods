%---------------------------------------------------------------------%
%This code computes the 1D Shallow Water Equations using either CG or 
%DG method with either LG or LGL points for interpolation and integration.
%The time-integration is accomplished via 2nd, 3rd Order, or 3rd Order 4-stage
%RK.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%

%----------Snippet of Code and how to use COMPUTE_TI_Coefficients-------%

%Compute RK Time-Integration Coefficients
kstages = 4; %=1 is RK1; =2 is RK2; = 3 is RK3; and = 4 is RK4
[a0,a1,beta] = compute_ti_coefficients(kstages);

%Time Integration
for itime=1:ntime
   time=time + dt;
   
   for ik=1:kstages %RK Stages
      
      %Create RHS Matrix
      rhs = create_rhs();
      
      %Solve System
      qp=a0(ik)*q0 + a1(ik)*q1 + dt*beta(ik)*rhs;
                  
      %Update
      rhs_rk4(:,:,:,ik)=rhs(:,:,:);
      q1=qp;
   end %ik
   
   %Do RK4 Solution
   if (kstages == 4)
       qp(:,:,:)=q0(:,:,:) + dt/6.0*( rhs_rk4(:,:,:,1) + 2*rhs_rk4(:,:,:,2)...
                      + 2*rhs_rk4(:,:,:,3) + rhs_rk4(:,:,:,4) );
   end
   
   %Update Q
   q0=qp;
   
end