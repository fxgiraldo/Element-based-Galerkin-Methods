%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [icorr,ncorr] = balance_coarsen_ratio(icor,ncor,face,jeside,tl,tp,tc,ngl)

p=0;
r=0;
icorr = 0;
ncorr = 0;

[fe,~] = get_pointers(ngl);

for i=1:ncor
    parent = tp(icor(i));
    clev = tl(icor(i));
    conflict = 0;
    
    %go over external faces
    for f=1:4
        ch(1) = tc(fe(f,1),parent);
        ch(2) = tc(fe(f,2),parent);
        
        for j=1:2
            ifex = jeside(ch(j),f);
            %locate neighbours
            if face(ifex,5)==ch(j) || face(ifex,7)==ch(j)
                   ngb1 = face(ifex,6);
                   ngb2 = face(ifex,8);
            elseif face(ifex,6)==ch(j) || face(ifex,8)==ch(j)
                   ngb1 = face(ifex,5);
                   ngb2 = face(ifex,7);
            end
            
            if tl(ngb1)==clev || tl(ngb1)==clev-1
                continue
            else
                conflict = 1; 
            end            
        end
        
    end
    
    if conflict == 0
       r=r+1;
       icorr(r) = icor(i);
    end
    
end

ncorp=p;
ncorr=r;
end
