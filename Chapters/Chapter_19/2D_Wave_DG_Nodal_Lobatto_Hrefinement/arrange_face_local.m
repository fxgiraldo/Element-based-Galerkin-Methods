%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function face = arrange_face_local(face,f)

    if face(f,9)==1 && face(f,8)==0
           ch1 = face(f,3);
           face(f,3) = face(f,4);
           face(f,4) = ch1;
           ch1 = face(f,5);
           face(f,5) = face(f,6);
           face(f,8) = face(f,7);
           face(f,6) = ch1;
           face(f,7) = 0;
    end 
end
