%---------------------------------------------------------------------%
%This function computes the Courant Number on Quadrilaterals.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qi_movie,time_movie,iframe] = update_movie(qi_movie,time_movie,q0,coord,npoin,nelem,ngl,ntime,itime,iplot,iframe,xi,yi)


if (mod(itime,iplot) == 0 || itime==ntime)
    iframe=iframe + 1;

   %Compute gridpoint solution
    q_sol=zeros(npoin,1);
    lhowm=zeros(npoin,1);
    for ie=1:nelem
        for j=1:ngl
        for i=1:ngl
          ip=intma(ie,i,j);
          lhowm(ip)=lhowm(ip)+1;
          q_sol(ip)=q_sol(ip) + q0(ie,i,j);
        end %i
        end %j
    end
    for i=1:npoin
       q_sol(i)=q_sol(i)/lhowm(i);
    end
   qi_movie(:,:,iframe)=griddata(coord(:,1),coord(:,2),q_sol,xi,yi,'cubic');
   time_movie(iframe)=timec;
end

