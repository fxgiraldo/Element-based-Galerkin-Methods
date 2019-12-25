%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [refel,refpt,coord,intma,face,facepa,iface,jeside,jesideh,iperiodic,tc,tp,tm,tl,new_el,nnew,nx,ny,jac_side] = refine_CG_DG(iref,nref,refel,refpt,coord,...
                                                intma,ngl,xgl,jeside,jesideh,iface,face,facepa,iperiodic,tc,tp,tm,tl,nq,wnq,psi,dpsi,nx,ny,jac_side)

k=0;
l=0;
m=0;
new_el = 0;
nnew = 0;

for i=1:nref
   parent = iref(i);
   if(tc(1,parent)==0)
        k=k+1;
        create(k) = parent;
   else
        l=l+1;
        refine(l) = parent;
   end
end
ncreate = k;
nrefine = l;

if ncreate>0
%     create
    iface_old = iface;
    [refel,refpt,coord,intma,face,facepa,iface,jeside,jesideh,iperiodic,tc,tp,tm,tl,ft,nft] = refine_elements_CG_DG_new(create,refel,refpt,coord,...
                                                intma,ngl,xgl,jeside,jesideh,iface,face,facepa,iperiodic,tc,tp,tm,tl);
    % compute normals          
    [nx,ny,jac_side]=compute_normals_local1(face,intma,coord,...
                          ngl,nq,wnq,psi,dpsi,nx,ny,jac_side,ft,nft);
    [nx,ny,jac_side]=compute_normals_local(face,intma,coord,...
                          ngl,nq,wnq,psi,dpsi,nx,ny,jac_side,iface,iface_old);
%     face = arrange_face(face,iface);
    
    for i=1:ncreate
        for j=1:4
            m=m+1;
            new_el(m) = tc(j,create(i));
            
%             %arrange face
%             for f=1:4
%                 face = arrange_face_local(face,jeside(new_el(m),f));
%             end
        end
    end
end
nnew = m;
if nrefine>0
%     refine
    [jeside,face,facepa,tm,ft,nft] = refine_elements_tree(refine,nrefine,jeside,jesideh,face,facepa,ngl,tc,tm);
    
    % compute normals          
    [nx,ny,jac_side]=compute_normals_local1(face,intma,coord,...
                          ngl,nq,wnq,psi,dpsi,nx,ny,jac_side,ft,nft);
end
end
