%---------------------------------------------------------------------%
%This function applies the Legendre Filter.
%Written by F.X. Giraldo on 7/2007
%           Department of Applied Maths
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function qp = apply_filter_dg(qp,f,nelem,ngl)

%Initialize
p_e=zeros(ngl,1);
pu_e=zeros(ngl,1);
pf=zeros(ngl,1);
puf=zeros(ngl,1);

%Integrate Divergence of Flux
for ie=1:nelem
   
   %Store Local Solution
   for i=1:ngl
      p_e(i)=qp(1,i,ie);
      pu_e(i)=qp(2,i,ie);
   end

   %Form Filtered Variable
   pf=f*p_e;
   puf=f*pu_e;

   %Store It
   for i=1:ngl
      qp(1,i,ie)=pf(i);
      qp(2,i,ie)=puf(i);
   end
end %ie
