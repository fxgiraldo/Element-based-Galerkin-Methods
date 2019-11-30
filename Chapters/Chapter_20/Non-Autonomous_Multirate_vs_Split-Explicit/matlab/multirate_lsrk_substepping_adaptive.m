function [q,Msteps_max] = multirate_lsrk_substepping_adaptive(q,c,alphaI,betaI,cI,alphaJ,betaJ,cJ,dt,t0,DT_Ratio,Msteps_max)

I = length(alphaI);
J = length(alphaJ);
M_local=0;

% Outer Loop 
r_slow = 0;
r_fast = 0;
cI = [cI;1];
for i = 1:I
  % slow mode
  t=t0 + cI(i) * dt;
  [~, ~, K_slow] = rhs_function(q,c,t);
  r_slow = alphaI(i) * r_slow + K_slow;

  % fast mode
  c_slow = (cI(i+1) - cI(i));
  d_slow = betaI(i) / c_slow;
  
%   dt_fast_real=c_slow*dt;
%   M=ceil(dt_fast_real/dt_fast_ideal);
  M=ceil(DT_Ratio*c_slow);
  M_local=M_local + M;
    
  %Substepping Loop
  for m=1:M
      for j = 1:J
          t=t0 + cI(i)*dt + c_slow*cJ(j)*dt*m/M;
          [~, K_fast, ~] = rhs_function(q,c,t);
          r_fast = alphaJ(j) * r_fast + K_fast + d_slow * r_slow;
          q = q + dt/M * c_slow * betaJ(j) * r_fast;
      end %j
  end %m
end %i

Msteps_max=max(Msteps_max,M_local);
