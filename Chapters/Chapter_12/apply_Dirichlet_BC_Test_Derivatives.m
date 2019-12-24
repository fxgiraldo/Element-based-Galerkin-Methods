%----------------------------------------------------------------------%
%This subroutine builds the FLUX vector for the Strong Form DGM-SEM
%on Quadrilateral Elements for the 2D Euler Equations.
%Written by Francis X. Giraldo on 1/2001
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%----------------------------------------------------------------------%
function [Lmatrix] = apply_Dirichlet_BC_Test_Derivatives(Lmatrix,psideh,...
               nside,ngl,imapl,intma)

%Construct FVM-type Operators
for is=1:nside

   %Store Left Side Variables
   el=psideh(is,3);
   er=psideh(is,4);
   if (er == -4) %Dirichlet bc
      ilocl=psideh(is,1);
      for l=1:ngl
          %Get Pointers
          il=imapl(ilocl,1,l);
          jl=imapl(ilocl,2,l);
          I=intma(el,il,jl);
          
          %Left Element
          Lmatrix(I,:)=0;
          Lmatrix(I,I)=1;
      end %l  
   end %if er      
end %is

