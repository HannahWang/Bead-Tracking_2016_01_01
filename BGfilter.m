function BGfilter(t,z, k1, nor)
I = imread(getBFfileName(z,t));
k2=0.2;
threshold = 170;
ll=4;
eraseS=9;
if t == 1
    I_BG3 = imread('CELL/Cell_BG.tif');
    I_BG3 = im2uint16(I_BG3)/300;
    I_BG168 = imread(getBFfileName(168,1));
    I_BG169_2 = imread(getBFfileName(169,1));
    I_BG169 = (I_BG168+I_BG169_2)/2;
    I_BG003 = imread(getBFfileName(3,1));
    I_BG004 = imread(getBFfileName(4,1));
    I_BG001 = (I_BG003+I_BG004)/2;
    I_BG=I_BG3*k2+I_BG169*k1+I_BG001*(1-k1);
elseif t == 2
    I_BG3 = imread('CELL/Cell_BG.tif');
    I_BG3 = im2uint16(I_BG3)/300;
    I_BG168=imread(getBFfileName(168,2));
    I_BG169_2=imread(getBFfileName(169,2));
    I_BG169 = (I_BG168+I_BG169_2)/2;
    I_BG003=imread(getBFfileName(3,2));
    I_BG004=imread(getBFfileName(4,2));
    I_BG001 = (I_BG003+I_BG004)/2;
    I_BG=I_BG3*k2+I_BG169*k1+I_BG001*(1-k1);
end
BG_weight = sum(I(1:10,1:10))/sum(I_BG(1:10,1:10));
I_nor=nor*(I-I_BG*BG_weight);
h = fspecial('gaussian');
I_blur = imfilter(I_nor, h);
im_fil2=medfilt2(I_blur); %necessary
im_threshold = im_fil2 - threshold*ones(1002,1004,'uint16'); %necessary
BW1 = edge(im_threshold, 'sobel');

se90=strel('line',ll,90);
se0=strel('line',ll,0);
%sdisk = strel('disk',8,8);
BWdfill = imdilate(BW1,[se90 se0]);
%BWsdil = imdilate(BW1,sdisk);
%BWdfill = imfill(BWsdil,'holes');

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
title(['fig',num2str(z)],'color','b');
%imwrite(BWdfill,sprintf('Copy_of_StrainEnergy3D_SD_2015_12_31/CELL/Cell_%d_%d.tif', t, z))

end
