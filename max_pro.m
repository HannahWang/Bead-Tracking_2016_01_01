% raw image z projection
% t=1;
% z_min=2;
% z_max=397;
% for k=z_min:z_max
% I=imread(sprintf('Copy_of_StrainEnergy3D_SD_2016-01-01/BFFRAME/BFframe_t%06i_%04i.tif',t,k));    
% comp_im(:,:,k-z_min+1)=I; 
% end
% max_im=max(comp_im,[],3);
t=1;
result=[];
layer=[2 99;100 199;200 299;300 397;2 397];
for i=1:4
z_min=layer(i,1);
z_max=layer(i,2);
sampling=10;
[x,y,z,u,v,w,lar_x,lar_y,lar_z,lar_u,lar_v,lar_w,count]=original_coor_ty(t,z_min,z_max,sampling);
[a,b]=size(x);
[c,d]=size(lar_x);
num=b+d;
     if i==5
         

         result(i,:)=[mean(u) mean(v) mean(w) count];
     else
        m_u=sqrt((sum(u.*u)+sum(lar_u.*lar_u))/num);
        m_v=sqrt((sum(v.*v)+sum(lar_v.*lar_v))/num);
        m_w=sqrt((sum(w.*w)+sum(lar_w.*lar_w))/num);
        result(i,:)=[m_u m_v m_w count];

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