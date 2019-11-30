function [q] = multirate_lsrk(q,c,alphaI,betaI,cI,alphaJ,betaJ,cJ,dt,t0)

I = length(alphaI);
J = length(alphaJ);
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
  for j = 1:J
    [~, K_fast, ~] = rhs_function(q,c,t + c_slow * cJ(j) * dt);
    r_fast = alphaJ(j) * r_fast + K_fast + d_slow * r_slow;
    q = q + dt * c_slow * betaJ(j) * r_fast;
  end
end
