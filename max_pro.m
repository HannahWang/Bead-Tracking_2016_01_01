% raw image z projection
% t=1;
% z_min=2;
% z_max=397;
% for k=z_min:z_max
% I=imread(sprintf('Copy_of_StrainEnergy3D_SD_2016-01-01/BFFRAME/BFframe_t%06i_%04i.tif',t,k));    
% comp_im(:,:,k-z_min+1)=I; 
% end
% max_im=max(comp_im,[],3);

% sampling=1;
% for thickness=5:5:30
% t=1;
% for z_min=200:200
%     z_max=z_min+thickness-1;
%     vector2contour(t,z_min,z_max,sampling,thickness);
% end
% end


t=1;
result=[];
layer=[2 397;2 99;100 199;200 299;300 397];
for i_layer=1:5
z_min=layer(i_layer,1);
z_max=layer(i_layer,2);
sampling=30;%30
padu = 1;
[x_padu,y_padu,z_padu,u_padu,v_padu,w_padu,count]=bemove_col4(t,z_min,z_max,sampling,padu,i_layer);
%[x,y,z,u,v,w,lar_x,lar_y,lar_z,lar_u,lar_v,lar_w,count]=bemove_col4(t,z_min,z_max,sampling,padu);

[a,b]=size(x);
num=b;
     if i_layer==1
     
         result(i_layer,:)=[mean(u) mean(v) mean(w) count];
     else
        m_u=sqrt((sum(u.*u)+sum(lar_u.*lar_u))/num);
        m_v=sqrt((sum(v.*v)+sum(lar_v.*lar_v))/num);
        m_w=sqrt((sum(w.*w)+sum(lar_w.*lar_w))/num);
        result(i_layer,:)=[m_u m_v m_w count];

     end
end
disp(result);
% figure,     
% imshow(max_im*30)
% hold on
% small_move=quiver3(x,y,z,u,v,w,'color','red');
% hold on
% Lar_move=quiver3(lar_x,lar_y,lar_z,lar_u,lar_v,lar_w,'color','green');
% hold on