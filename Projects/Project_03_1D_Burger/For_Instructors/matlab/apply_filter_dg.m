%---------------------------------------------------------------------%
%This function applies the Legendre Filter.
%Written by F.X. Giraldo on 7/2007
%           Department of Applied Maths
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function qp = apply_filter_dg(qp,f,intma,nelem,ngl)

%Initialize
h_e=zeros(ngl,1);
U_e=zeros(ngl,1);
hf=zeros(ngl,1);
Uf=zeros(ngl,1);

%Integrate Divergence of Flux
for e=1:nelem
   
   %Store Local Solution
   for i=1:ngl
      ip=intma(i,e);
      h_e(i)=qp(ip,1);
      U_e(i)=qp(ip,2);
   end

   %Form Filtered Variable
   hf=f*h_e;
   Uf=f*U_e;

   %Store It
   for i=1:ngl
      ip=intma(i,e);
      qp(ip,1)=hf(i);
      qp(ip,2)=Uf(i);
   end
end %ie
