%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function Fmatrix = Fmatrix_centered_flux(Fmatrix,intma,nelem,npoin,ngl,u)

Fmatrix=0;
for e=1:nelem
      
    %Visit Left-most DOF of each element
    i=1;
    I=intma(i,e);
    %shift left
    Im=I-1;
    if (Im < 1) 
        Im=npoin; %periodicity
    end
    Fmatrix(I,Im)=-0.5;
    Fmatrix(I,I)=-0.5;
    
    %Visit Right-most DOF of each element
    i=ngl;
    I=intma(i,e);
    %shift left
    Ip=I+1;
    if (Ip > npoin) 
        Ip=1; %periodicity
    end
    Fmatrix(I,Ip)=0.5;
    Fmatrix(I,I)=0.5;
   
end %e

Fmatrix=Fmatrix*u;


      