%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on April 22, 2021
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function Fmatrix = Fmatrix_upwind_flux(intma,nelem,npoin,ngl,u)

Fmatrix=zeros(npoin,npoin);

for e=1:nelem
      
    %Visit Left-most DOF of each element
    i=1;
    I=intma(i,e);
    %shift left
    Im=I-1;
    if (Im < 1) 
        Im=npoin; %periodicity
    end
    Fmatrix(I,Im)=-1;
    
    %Visit Right-most DOF of each element
    i=ngl;
    I=intma(i,e);
    %shift left
    Ip=I+1;
    if (Ip > npoin) 
        Ip=1; %periodicity
    end
    Fmatrix(I,I)=1;
   
end %e

Fmatrix=Fmatrix*u;


      