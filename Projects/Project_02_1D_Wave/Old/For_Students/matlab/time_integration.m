%---------------------------------------------------------------------%
%This function performs the time-integration
%Written by F.X. Giraldo on April 19, 2019
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey; CA 93943-5216
%---------------------------------------------------------------------%
function [q0,time] = time_integration(q0,Dhat,periodicity,time,ntime,dt)
    
[q0,time] = ti_LSRK(q0,Dhat,periodicity,time,ntime,dt);
