%---------------------------------------------------------------------%
%This function applies the Legendre Filter.
%Written by F.X. Giraldo on 7/2007
%           Department of Applied Maths
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function qp = apply_filter_dg(qp,f,nelem,ngl)

%Initialize
h_e=zeros(ngl,1);
U_e=zeros(ngl,1);
hf=zeros(ngl,1);
Uf=zeros(ngl,1);

%Integrate Divergence of Flux
for e=1:nelem
   
   %Store Local Solution
   h_e(:)=qp(1,:,e);
   U_e(:)=qp(2,:,e);
%    for i=1:ngl
%       h_e(i)=qp(1,i,e);
%       U_e(i)=qp(2,i,e);
%    end

   %Form Filtered Variable
   hf=f*h_e;
   Uf=f*U_e;

   qp(1,:,e)=hf(:);
   qp(2,:,e)=Uf(:);
%    %Store It
%    for i=1:ngl
%       qp(1,i,e)=hf(i);
%       qp(2,i,e)=Uf(i);
%    end
end %ie
