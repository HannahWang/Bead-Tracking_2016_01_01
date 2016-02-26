% max intensity z projection
%function maxpro_im(t,z_min,z_max)
max_pro=zeros(1002,1004);
z_min=3;
z_max=4;
t=1;
for k=z_min:z_max
    for i=1:1002
        for j=1:1004
    %comp_image=imread(sprintf('Copy_of_StrainEnergy3D_SD_2015_12_31/BFFRAME/BFframe_t%06i_%04i.tif',t,k));
    comp_image=imread(getBFfileName(k,t));
         if max_pro(i,j)<comp_image(i,j)
            max_pro(i,j)=comp_image(i,j); 
         end
         %max_pro=max_pro*50;
    
   % imshowpair(BWdfill,I,'montage');
        end
    end
end
      imshow(max_pro)