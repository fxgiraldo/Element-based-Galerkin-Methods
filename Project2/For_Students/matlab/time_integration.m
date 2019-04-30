%---------------------------------------------------------------------%
%This function performs the time-integration
%Written by F.X. Giraldo on April 19, 2019
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey; CA 93943-5216
%---------------------------------------------------------------------%
function [q0] = time_integration(q0,u,Dhat,Fhat,intma,periodicity,time,ntime,dt,stages,ti_type)

if (ti_type == 1)
    q0 = ti_SSP(q0,u,Dhat,Fhat,intma,periodicity,time,ntime,dt,stages);
elseif (ti_type == 2)
    q0 = ti_LSRK(q0,u,Dhat,Fhat,intma,periodicity,time,ntime,dt);
end