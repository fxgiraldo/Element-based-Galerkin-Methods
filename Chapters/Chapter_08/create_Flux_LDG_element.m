%---------------------------------------------------------------------%
%This function computes the Laplacian Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_Flux_LDG_element(rhs,intma,nelem,ngl,q,qe,alpha,beta)

%Build Q: FLUX
for e=1:nelem
   
   %-------I=0 Edge
   %Left Edge
   IL=intma(1,e);
   q_L=q(IL);
   n_L=-1; %(pointing from Left -> Right)
   
   %Right Edge
   if (e>1)
       IR=intma(ngl,e-1);
       q_R=q(IR);
   else
       q_R=qe(intma(1,1));
   end
   qstar=alpha*q_L + beta*q_R;
      
   %Add to Element
   rhs(IL)=rhs(IL) + n_L*qstar;
   
   %-------I=N Edge
   %Left Edge
   IL=intma(ngl,e);
   q_L=q(IL);
   n_L=+1; %(pointing from Left -> Right)
   
   %Right Edge
   if (e<nelem)
       IR=intma(1,e+1);
       q_R=q(IR);
   else
       q_R=qe(intma(ngl,nelem));
   end
   qstar=alpha*q_L + beta*q_R;
      
   %Add to Element
   rhs(IL)=rhs(IL) + n_L*qstar;
   
end %ie