%---------------------------------------------------------------------%
%This function computes the Laplacian Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_Flux_q(rhs,intma,nelem,ngl,q,qe,alpha,beta)

%Build Q: FLUX
for e=1:nelem-1
   el=e;
   er=e+1;
   
   %LGL Integration
   IL=intma(ngl,el);
   IR=intma(1,er);
   q_l=q(IL);
   q_r=q(IR);
   nvector=+1; %(pointing from Left -> Right)
   qstar=alpha*q_l + beta*q_r;
      
   %Add to Left
   rhs(IL)=rhs(IL) + nvector*qstar;
   
   %Add to Right
   rhs(IR)=rhs(IR) - nvector*qstar;
end %ie

%Right Lateral Boundary
IL=intma(ngl,nelem);
q_l=q(IL);
q_r=qe(IL);
nvector=+1; %(pointing from Left -> Right)
qstar=alpha*q_l + beta*q_r;
rhs(IL)=rhs(IL) + nvector*qstar;

%Left Lateral Boundary
IL=intma(1,1);
q_l=q(IL);
q_r=qe(IL);
nvector=+1; %(pointing from Left -> Right)
qstar=alpha*q_l + beta*q_r;
rhs(IL)=rhs(IL) - nvector*qstar;
