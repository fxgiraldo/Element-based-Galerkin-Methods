%---------------------------------------------------------------------%
%This code builds the 2:1 Balanced H-refinement non-conforming AMR algorithm
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [refel,refpt,coord,intma,face,facepa,iface,jeside,jesideh,...
          iperiodic,tc,tp,tm,tl,nx,ny,jac_side,ksi_x,ksi_y,...
          eta_x,eta_y,jac,Mmatrix,Mmatrix_inv,qp,q0,q1,ue,ve,ffc,nffc,sfc,...
          nsfc,Courant,vel,ds,dt,nelemT,npoinT,nface] = ...
                    AMR(ref_threshold,cor_threshold,iadapt,icoarse,itime,ref_freq,...
                        sfc,nsfc,ffc,nffc,q0,q1,qp,ue,ve,refel,refpt,xgl,...
                        jeside,jesideh,iface,face,facepa,iperiodic,tc,tp,tm,tl,...
                        nq,wnq,psi,dpsi,nx,ny,jac_side,intma,coord,Courant_max,...
                        maxlev,nelem,nface,ngl,ksi_x,ksi_y,eta_x,eta_y,jac,...
                        Mmatrix,Mmatrix_inv,nelemT,npoinT,Courant,vel,ds,dt,...
                        P1s2d,P2s2d,P3s2d,P4s2d,P1g2d,P2g2d,P3g2d,P4g2d)

   %--- AMR ---
   
      %mark for mesh adaptation
   if(mod(itime,ref_freq) == 0 && iadapt == 1)
       iref = 0;
       ir = 0;
       icor = 0;
       ic = 0;  
       nrefr=0;
       
       
       for is=1:nsfc
           ie = sfc(is);
           clev = tl(ie);
           if clev<maxlev
               qoi = max(max(q0(ie,:,:)));
               if (qoi > ref_threshold(clev+1))
                   ir=ir+1;
                   iref(ir) = ie;
    %                disp(['ref ',int2str(ie),' qoi ',num2str(qoi)]);
                    continue;
               end
           end
           if clev>0
               if clev==maxlev
                   qoi = max(max(q0(ie,:,:)));
               end
               
               if (qoi < cor_threshold(clev) && icoarse==1)
                   ic=ic+1;
                   icor(ic) = ie;
    %                disp(['cor ',int2str(ie),' parent ',num2str(tp(ie)),' qoi ',num2str(qoi)]);
                   continue;
               end
           end
           
               
       end

       ncor = ic;
       nref = ir;
     
       
%             %make sure icor does not have any siblings
%             [icor,ncor] = correct_coarse_list_new(icor,ncor,tp,tc);
           
            
       %refinement first
       
       %go over levels in ascending order
       
%        [irefl,nrefl] = extract_level(iref,nref,tl,0); %zeroth-level first
%        if nrefl>0
%             [refel,refpt,coord,intma,face,facepa,iface,jeside,jesideh,iperiodic,tc,tp,tm,tl,new_el,nnew,nx,ny,jac_side] = refine_CG_DG(irefl,nrefl,...
%                                             refel,refpt,coord,intma,ngl,xgl,jeside,jesideh,iface,face,facepa,iperiodic,tc,tp,tm,tl,nq,wnq,psi,...
%                                             dpsi,nx,ny,jac_side);
%             nface = iface;
%             nelemT=refel-1;
%             npoinT=refpt-1;
% 
%             if nnew>0
%                 [ksi_x,ksi_y,eta_x,eta_y,jac] = metrics_local(coord,intma,psi,dpsi,wnq,ngl,nq,new_el,nnew,ksi_x,ksi_y,eta_x,eta_y,jac);
%                 [Mmatrix,Mmatrix_inv] = create_Mmatrix2d_dg_TensorProduct_inexact_local(Mmatrix,Mmatrix_inv,jac,ngl,new_el,nnew);
%             end
% 
%             [qp] = refine_project_data(qp,irefl,nrefl,P1s2d,P2s2d,P3s2d,P4s2d,ngl,tc);
%             [q1] = refine_project_data(q1,irefl,nrefl,P1s2d,P2s2d,P3s2d,P4s2d,ngl,tc);
%             [ue] = refine_project_data(ue,irefl,nrefl,P1s2d,P2s2d,P3s2d,P4s2d,ngl,tc);
%             [ve] = refine_project_data(ve,irefl,nrefl,P1s2d,P2s2d,P3s2d,P4s2d,ngl,tc);                            
%        end
       

       %now the rest
       if nref>0
           for lev=0:maxlev-1
                [irefl,nrefl] = extract_level(iref,nref,tl,lev);            

                if nrefl>0

                    
                    if lev==0
                        %skip balance
                        irefr=irefl;
                        nrefr=nrefl;
                        nrefp=0;
                        irefp=0;
                    else
                        %balance
                        [irefp,nrefp,irefr,nrefr] = balance_refinement_ratio(irefl,nrefl,face,jeside,tl,0,maxlev);
                    end

                    if nrefp>0
                       disp('unresolved refinement problem!'); 
                    end
                                        
                     for i=1:nrefr
                        [refel,refpt,coord,intma,face,facepa,iface,jeside,jesideh,iperiodic,tc,tp,tm,tl,new_el,nnew,nx,ny,jac_side] = refine_CG_DG(irefr(i),1,...
                                                        refel,refpt,coord,intma,ngl,xgl,jeside,jesideh,iface,face,facepa,iperiodic,tc,tp,tm,tl,nq,wnq,psi,...
                                                        dpsi,nx,ny,jac_side);
                                                 
                    nface = iface;
                    nelemT=refel-1;
                    npoinT=refpt-1;
                    

                    if nnew>0
                        [ksi_x,ksi_y,eta_x,eta_y,jac] = metrics_local(coord,intma,psi,dpsi,wnq,ngl,nq,new_el,nnew,ksi_x,ksi_y,eta_x,eta_y,jac);
                        [Mmatrix,Mmatrix_inv] = create_Mmatrix2d_dg_TensorProduct_inexact_local(Mmatrix,Mmatrix_inv,jac,ngl,new_el,nnew);
                    end
                     end
                    
                    [qp] = refine_project_data(qp,irefr,nrefr,P1s2d,P2s2d,P3s2d,P4s2d,ngl,tc);
                    [q0] = refine_project_data(q0,irefr,nrefr,P1s2d,P2s2d,P3s2d,P4s2d,ngl,tc);
                    [q1] = refine_project_data(q1,irefr,nrefr,P1s2d,P2s2d,P3s2d,P4s2d,ngl,tc);
                    [ue] = refine_project_data(ue,irefr,nrefr,P1s2d,P2s2d,P3s2d,P4s2d,ngl,tc);
                    [ve] = refine_project_data(ve,irefr,nrefr,P1s2d,P2s2d,P3s2d,P4s2d,ngl,tc);                            
                end

           end
       end
       
%        %de-refinement
% %               q0=qp;
%        
%        %Generate the Space Filling Curve - SFC
%        [sfc,nsfc] = space_filling_curve_tree(nelem,tc,tm);
%        [ffc,nffc] = face_filling_curve(nface,facepa);
%         for is=1:nsfc
%            ie = sfc(is);
%            clev = tl(ie);
%            
%            if clev>0
%                if clev==maxlev
%                    qoi = max(max(q0(ie,:,:)));
%                end
%                
%                if (qoi < cor_threshold(clev) && icoarse==1)
%                    ic=ic+1;
%                    icor(ic) = ie;
%     %                disp(['cor ',int2str(ie),' parent ',num2str(tp(ie)),' qoi ',num2str(qoi)]);
%                    continue;
%                end
%            end
%            
%                
%        end
% 
%        ncor = ic;
     
       
            %make sure icor does not have any siblings or inactive elements
            %(due to ripple propagation)
            [icora,ncora] = correct_coarse_list_new(icor,ncor,tp,tc,tm);
            
            %make sure icor does not contain any rippled or previously coelements
            if nrefr>nref
                k=0;
                icor = 0;
                for i=1:ncora
                   m=0;
                   for j=1:nrefr-nref
                       if tp(icora(i)) == tp(irefr(j))
                           m=1;
                       end
                   end
                   if  m==0
                       k=k+1;
                       icor(k)=icora(i);
                   end

                end
                ncor=k;
            else
                icor=icora;
                ncor=ncora;
            end
       
       if ncor>0
           
           for lev=maxlev:-1:1 %decreasing order
               [icorl,ncorl] = extract_level(icor,ncor,tl,lev);
               
                
               
               if ncorl>0
                   
                   if lev==maxlev
                       %skip balance
                       icorr=icorl;
                       ncorr=ncorl;
                   else
                       %balance
                       [icorr,ncorr] = balance_coarsen_ratio(icorl,ncorl,face,jeside,tl,tp,tc,ngl);
                   end
                   
                   [jeside,face,facepa,tm,ft,nft] = coarsen_elements_new(icorr,ncorr,jeside,jesideh,face,facepa,ngl,tc,tp,tm);
                   % compute normals          
                    [nx,ny,jac_side]=compute_normals_local1(face,intma,coord,...
                                  ngl,nq,wnq,psi,dpsi,nx,ny,jac_side,ft,nft);
                   
                   [qp] = coarse_project_data(qp,icorr,ncorr,P1g2d,P2g2d,P3g2d,P4g2d,ngl,tc,tp);
                   [q0] = coarse_project_data(q0,icorr,ncorr,P1g2d,P2g2d,P3g2d,P4g2d,ngl,tc,tp);
                   [q1] = coarse_project_data(q1,icorr,ncorr,P1g2d,P2g2d,P3g2d,P4g2d,ngl,tc,tp);
                   [ue] = coarse_project_data(ue,icorr,ncorr,P1g2d,P2g2d,P3g2d,P4g2d,ngl,tc,tp);
                   [ve] = coarse_project_data(ve,icorr,ncorr,P1g2d,P2g2d,P3g2d,P4g2d,ngl,tc,tp);
                   
               end
           end      
       end
            

%        if(nref>0)
%              iref
%             [refel,refpt,coord,intma,face,facepa,iface,jeside,jesideh,iperiodic,tc,tp,tm,tl,new_el,nnew,nx,ny,jac_side] = refine_CG_DG(iref,nref,...
%                                                 refel,refpt,coord,intma,ngl,xgl,jeside,jesideh,iface,face,facepa,iperiodic,tc,tp,tm,tl,nq,wnq,psi,...
%                                                 dpsi,nx,ny,jac_side);
% 
%             nface = iface;
%             nelemT=refel-1;
%             npoinT=refpt-1;
%             
%             if nnew>0
%                 [ksi_x,ksi_y,eta_x,eta_y,jac] = metrics_local(coord,intma,psi,dpsi,wnq,ngl,nq,new_el,nnew,ksi_x,ksi_y,eta_x,eta_y,jac);
%                 [Mmatrix,Mmatrix_inv] = create_Mmatrix2d_dg_TensorProduct_inexact_local(Mmatrix,Mmatrix_inv,jac,ngl,new_el,nnew);
%             end
%            
% 
%             [qp] = refine_project_data(qp,iref,nref,P1s2d,P2s2d,P3s2d,P4s2d,ngl,tc);
%             [q1] = refine_project_data(q1,iref,nref,P1s2d,P2s2d,P3s2d,P4s2d,ngl,tc);
%             [ue] = refine_project_data(ue,iref,nref,P1s2d,P2s2d,P3s2d,P4s2d,ngl,tc);
%             [ve] = refine_project_data(ve,iref,nref,P1s2d,P2s2d,P3s2d,P4s2d,ngl,tc);
%        end
%        
%        if nref1>0
%             [refel,refpt,coord,intma,face,facepa,iface,jeside,jesideh,iperiodic,tc,tp,tm,tl,new_el,nnew,nx,ny,jac_side] = refine_CG_DG(iref1,nref1,...
%                                                 refel,refpt,coord,intma,ngl,xgl,jeside,jesideh,iface,face,facepa,iperiodic,tc,tp,tm,tl,nq,wnq,psi,...
%                                                 dpsi,nx,ny,jac_side);
% 
%             nface = iface;
%             nelemT=refel-1;
%             npoinT=refpt-1;
% 
%             if nnew>0
%                 [ksi_x,ksi_y,eta_x,eta_y,jac] = metrics_local(coord,intma,psi,dpsi,wnq,ngl,nq,new_el,nnew,ksi_x,ksi_y,eta_x,eta_y,jac);
%                 [Mmatrix,Mmatrix_inv] = create_Mmatrix2d_dg_TensorProduct_inexact_local(Mmatrix,Mmatrix_inv,jac,ngl,new_el,nnew);
%             end
% 
% 
%             [qp] = refine_project_data(qp,iref1,nref1,P1s2d,P2s2d,P3s2d,P4s2d,ngl,tc);
%             [q1] = refine_project_data(q1,iref1,nref1,P1s2d,P2s2d,P3s2d,P4s2d,ngl,tc);
%             [ue] = refine_project_data(ue,iref1,nref1,P1s2d,P2s2d,P3s2d,P4s2d,ngl,tc);
%             [ve] = refine_project_data(ve,iref1,nref1,P1s2d,P2s2d,P3s2d,P4s2d,ngl,tc);
% 
%        end
% 
%        if(ncor1>0)
% %            icor
%            
%            [jeside,face,facepa,tm,ft,nft] = coarsen_elements_new(icor1,ncor1,jeside,jesideh,face,facepa,ngl,tc,tp,tm);
%            [qp] = coarse_project_data(qp,icor1,ncor1,P1g2d,P2g2d,P3g2d,P4g2d,ngl,tc,tp);
%            [q1] = coarse_project_data(q1,icor1,ncor1,P1g2d,P2g2d,P3g2d,P4g2d,ngl,tc,tp);
%            [ue] = coarse_project_data(ue,icor1,ncor1,P1g2d,P2g2d,P3g2d,P4g2d,ngl,tc,tp);
%            [ve] = coarse_project_data(ve,icor1,ncor1,P1g2d,P2g2d,P3g2d,P4g2d,ngl,tc,tp);
%            % compute normals          
%             [nx,ny,jac_side]=compute_normals_local1(face,intma,coord,...
%                           ngl,nq,wnq,psi,dpsi,nx,ny,jac_side,ft,nft);
%        end
%        
%        if(ncor>0)
% %            icor
%            
%            [jeside,face,facepa,tm,ft,nft] = coarsen_elements_new(icor,ncor,jeside,jesideh,face,facepa,ngl,tc,tp,tm);
%            [qp] = coarse_project_data(qp,icor,ncor,P1g2d,P2g2d,P3g2d,P4g2d,ngl,tc,tp);
%            [q1] = coarse_project_data(q1,icor,ncor,P1g2d,P2g2d,P3g2d,P4g2d,ngl,tc,tp);
%            [ue] = coarse_project_data(ue,icor,ncor,P1g2d,P2g2d,P3g2d,P4g2d,ngl,tc,tp);
%            [ve] = coarse_project_data(ve,icor,ncor,P1g2d,P2g2d,P3g2d,P4g2d,ngl,tc,tp);
%            % compute normals          
%             [nx,ny,jac_side]=compute_normals_local1(face,intma,coord,...
%                           ngl,nq,wnq,psi,dpsi,nx,ny,jac_side,ft,nft);
%        end
%        q0=qp;
       
       %Generate the Space Filling Curve - SFC
       [sfc,nsfc] = space_filling_curve_tree(nelem,tc,tm);
       [ffc,nffc] = face_filling_curve(nface,facepa);
       
       [Courant,vel,ds,dt] = compute_Courant(ue,ve,intma,coord,ngl,dt,Courant_max,sfc,nsfc);


   end
   
   


end
