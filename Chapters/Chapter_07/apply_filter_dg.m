%---------------------------------------------------------------------%
%This function applies the Legendre Filter.
%Written by F.X. Giraldo on 7/2007
%           Department of Applied Maths
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = apply_filter_dg(qp,f,nelem,ngl)

%Initialize
rhs=zeros(ngl,nelem);
q_e=zeros(ngl);
qf=zeros(ngl);

%Integrate Divergence of Flux
for ie=1:nelem
   
   %Store Local Solution
   for i=1:ngl
      q_e(i)=qp(i,ie);
   end

   %Form Filtered Variable
   qf=f*q_e;

   %Store It
   for i=1:ngl
      rhs(i,ie)=qf(i);
   end
end %ie
