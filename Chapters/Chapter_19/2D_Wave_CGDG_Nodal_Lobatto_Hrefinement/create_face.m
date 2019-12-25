%---------------------------------------------------------------------%
%Written by F.X. Giraldo and M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [face,facepa,facep] = create_face(nside,iside,psideh,ae)

%construct new face structure from iside and psideh
face = zeros(nside*ae,11);
facep = zeros(nside*ae,1);
facepa = zeros(nside*ae,1);


for is = 1:nside
    face(is,1) = iside(is,1);
    face(is,2) = iside(is,2);
    
    switch psideh(is,1)
        case 1
            face(is,3) = 1;
        case 2
            face(is,3) = 4;
        case 3
            face(is,3) = 2;
        case 4
            face(is,3) = 3;   
    end
    
    switch psideh(is,2)
        case 1
            face(is,4) = 1;
        case 2
            face(is,4) = 4;
        case 3
            face(is,4) = 2;
        case 4
            face(is,4) = 3;   
    end
%     if psideh(is,3)<0
%         face(is,5) = 0;
%         face(is,9) = psideh(is,3);
%     else
        face(is,5) = psideh(is,3);
%     end
%     if psideh(is,4)<0
%         face(is,7) = 0;
%         face(is,9) = psideh(is,4);
%     else
        face(is,6) = psideh(is,4); %changed here 7 to 6
%     end
    face(is,10) = psideh(is,5);
    facepa(is)=1;
    if(psideh(is,3)<0)
       facep(is)=psideh(is,3); 
       facepa(is)=psideh(is,3); 
    end
    face(is,11) = psideh(is,6);
end

end
