%---------------------------------------------------------------------%
%This function updates/stores the frames of the movie for a CG/DG code.
%Written by F.X. Giraldo on September 16, 2014
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%Input Parameters: q0(nelem,ngl,ngl) is the solution
%                  coord(npoin,2) are the coordinates from CREATE_GRID_2D
%                  intma(nelem,ngl,ngl) is the element connectivity
%                  npoin is the number of points
%                  nelem is the number of elements
%                  ngl is the number of interpolation points in an element
%                  in one direction
%                  ntime is the number of time-steps to reach the final
%                  solution (e.g., ntime=time_final/dt)
%                  itime is the time-loop index (e.g., itime=1:ntime)
%                  nframes are the number of frames you want to store
%                  (user's choice)
%                  iframe is the current frame number.
%                  xi yi are the coordinates reshaped in INITIALIZE_MOVIE
%                  timec is the simulation time of the frame.
%Output Paramters: qi_movie(201,201,nframes+1) is the solution interpolated
%                  for plotting the movie
%                  time_movie(nframes+1,1) is the time for each frame
%                  iframe is the new frame number
%---------------------------------------------------------------------%
function [qi_movie,time_movie,iframe] = update_movie(qi_movie,time_movie,q0,coord,intma,npoin,nelem,ngl,ntime,itime,nframes,iframe,xi,yi,timec)

iplot=round(ntime/nframes);
if (mod(itime,iplot) == 0 || itime==ntime)
% if (mod(itime,iplot) == 0)
    iframe=iframe + 1;

   %Compute gridpoint solution
    q_sol=zeros(npoin,1);
    lhowm=zeros(npoin,1);
    for e=1:nelem
        for j=1:ngl
        for i=1:ngl
          I=intma(e,i,j);
          lhowm(I)=lhowm(I)+1;
          q_sol(I)=q_sol(I) + q0(e,i,j);
        end %i
        end %j
    end
    for i=1:npoin
       q_sol(i)=q_sol(i)/lhowm(i);
    end
   qi_movie(:,:,iframe)=griddata(coord(:,1),coord(:,2),q_sol,xi,yi,'cubic');
   time_movie(iframe)=timec;
end

