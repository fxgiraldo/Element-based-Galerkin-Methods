%---------------------------------------------------------------------%
%This function computes the Laplacian Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_Flux_q(rhs,intma,nelem,ngl,q,alpha,beta)

%Build Q: FLUX
for e=1:nelem
   eL=e;
   eR=e+1;
   if (eR > nelem) 
       eR=1;
   end
   
   %LGL Integration
   IL=intma(ngl,eL);
   IR=intma(1,eR);
   q_l=q(IL);
   q_r=q(IR);
   nvector=+1; %(pointing from Left -> Right)
   qstar=alpha*q_l + beta*q_r;
      
   %Add to Left
   rhs(IL)=rhs(IL) + nvector*qstar;
   
   %Add to Right
   rhs(IR)=rhs(IR) - nvector*qstar;
end %ie