function [q] = lsrk(q,c,alphaI,betaI,cI,dt,t0)

I = length(alphaI);
r = 0;
for i = 1:I
  [K, ~, ~] = rhs_function(q, c, t0 + dt * cI(i));
  r = alphaI(i) * r + K;
  q = q + betaI(i) * dt * r;
end
