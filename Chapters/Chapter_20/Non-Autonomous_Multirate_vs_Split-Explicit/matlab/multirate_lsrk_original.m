function [q] = multirate_lsrk_original(q, c, RKA, RKB, RKC, dt, t0)

I = length(RKA);
% Outer Loop 
r = 0;
r_slow = 0;
RKC = [RKC;1];
for i = 1:I
  %{
  [K ~, ~] = rhs_function(q, c, t0 + dt * RKC(i));
  r = RKA(i) * r + K;
  q = q + RKB(i) * dt * r;
  continue
  %}

  % slow mode
  [~, ~, K_slow] = rhs_function(q, c, t0 + RKC(i) * dt);
  r_slow = RKA(i) * r_slow + K_slow;

  % fast mode
  c_slow = (RKC(i+1) - RKC(i));
  d_slow = RKB(i) / c_slow;

  r_fast = 0;
  t = t0 + RKC(i) * dt;
  for j = 1:I
    [~, K_fast, ~] = rhs_function(q, c, t + c_slow * RKC(j) * dt);
    r_fast = RKA(j) * r_fast + K_fast + d_slow * r_slow;
    q = q + dt * c_slow * RKB(j) * r_fast;
  end
end
