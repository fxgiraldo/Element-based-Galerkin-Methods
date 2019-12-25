%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [irefp,nrefp,irefr,nrefr] = balance_refinement_ratio(iref,nref,face,jeside,tl,irec,maxlev)

%balance refinement
p=0;
r=0;
irefr = 0;
irefp = 0;
nrefp = 0;
nrefr = 0;
irec=irec+1;

if irec > maxlev;
    return
end

for i=1:nref
    parent = iref(i);
    if parent==0
        continue
    end
    
    clev = tl(parent);
    m=0;
    for f=1:4 %go over all faces of the element
        ifex = jeside(parent,f);
        
        %locate neighbours
        if face(ifex,5)==parent || face(ifex,7)==parent
               ngb1 = face(ifex,6);
%                ngb2 = face(ifex,8);
        elseif face(ifex,6)==parent || face(ifex,8)==parent
               ngb1 = face(ifex,5);
%                ngb2 = face(ifex,6);
        end
        
        if tl(ngb1)==clev || tl(ngb1)==clev+1
            continue
        else
           
            m=m+1;
            isinp=0;
            for j=1:p
                if ngb1 == irefp(j)
                    isinp = 1;
                end
            end
            
            if isinp==0
                p=p+1;
                irefp(p) = ngb1;
            end
            
            for j=1:nref
                if ngb1==iref(j)
                   iref(j) =0; %delete from original refinement list, if it was there 
                end
            end
        end
    end
%     if m==0
        r=r+1;
        irefr(r) = parent;
%     end
    
end
nrefr = r;
nrefp = p;


%if any conflicts still present resolve them recursively
if nrefp>0
    disp('balancing')
    [irefp,nrefp,irefr1,nrefr1] = balance_refinement_ratio(irefp,nrefp,face,jeside,tl,irec,maxlev);
    [irefr,nrefr] = append_to_front(irefr1,nrefr1,irefr,nrefr); %append to previously obtained list
end
   
end
