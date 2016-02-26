function TEST_FILTER(t)
for z=0:169
I = imread(getBFfileName(z,t));
if t == 1
    I_BG=imread('BFframe_t000001_0001.tif');
elseif t == 2
    I_BG=imread('BFframe_t000002_00000004.tif');
end
BG_weight = sum(I(1:10,1:10))/sum(I_BG(1:10,1:10));
I_nor=1.7*(I-I_BG*BG_weight);
im_fil2=medfilt2(I_nor);
im_fil2 = im_fil2 - 155*ones(1002,1004,'uint16'); %necessary
BW1 = edge(im_fil2, 'sobel');

se90=strel('line',3,90);
se0=strel('line',3,0);
%sdisk = strel('disk',8,8);
BWdfill = imdilate(BW1,[se90 se0]);
%BWsdil = imdilate(BW1,sdisk);
%BWdfill = imfill(BWsdil,'holes');

eraseS = 8;
for i=1:1002-eraseS
   for j=1:1004-eraseS
      if sum(BWdfill(i,j:j+eraseS))==0 && sum(BWdfill(i+eraseS,j:j+eraseS))==0 && sum(BWdfill(i:i+eraseS,j))==0 && sum(BWdfill(i:i+eraseS,j+eraseS))==0
         BWdfill(i:i+eraseS-1,j:j+eraseS-1)=0;
      elseif i==1 && sum(BWdfill(i:i+eraseS,j))+sum(BWdfill(i:i+eraseS,j+eraseS))+sum(BWdfill(i+eraseS,j:j+eraseS))==0
         BWdfill(i:i+eraseS-1, j:j+eraseS-1)=0;
      elseif i==1002-eraseS && sum(BWdfill(i:i+eraseS,j))+sum(BWdfill(i:i+eraseS,j+eraseS))+sum(BWdfill(i,j:j+eraseS))==0
         BWdfill(i:i+eraseS-1, j:j+eraseS-1)=0;
      elseif j==1 && sum(BWdfill(i,j:j+eraseS))+sum(BWdfill(i+eraseS,j:j+eraseS))+sum(BWdfill(i:i+eraseS,j+eraseS))==0
         BWdfill(i:i+eraseS-1, j:j+eraseS-1)=0;
      elseif j==1004-eraseS && sum(BWdfill(i,j:j+eraseS))+sum(BWdfill(i+eraseS,j:j+eraseS))+sum(BWdfill(i:i+eraseS,j))==0
         BWdfill(i:i+eraseS-1, j:j+eraseS-1)=0;
      end
   end
end


imshowpair(BWdfill,I,'montage');
imwrite(BWdfill,sprintf('Copy_of_StrainEnergy3D_SD_2015_12_31/Cell_%d_%d.tif', t, z))
end
end