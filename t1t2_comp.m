k1=0.245;
nor=2.45;

for t=1:2
    I = imread(sprintf('MAX_t%d_20_40.tif',t));
    %I2 = imread('MAX_t2_20_40.tif');
    k2=0.2;
    threshold = 170;
    ll=4;
    eraseS=9;
    if t == 1
        % set t=1 background
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
        % set t=2 background
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

%eval(['BWdfill_t',int2str(t),'=BWdfill;']);
    

% imshowpair(BWdfill_t1,BWdfill_t2,'diff');
% title(['fig',num2str(z)],'color','b');

    eval(['X',int2str(t),'=[];']);
    eval(['Y',int2str(t),'=[];']);
    eval(['Z',int2str(t),'=[];']);

 
  for i = 1:1002
     for j = 1:1004
        if BWdfill(i,j) == 1
            %for conti = 0:7
               
               eval(['X',int2str(t),'(end+1) = i/6;']);
               eval(['Y',int2str(t),'(end+1) = j/6;']); 
               %eval(['Z',int2str(t),'(end+1) = 19-',int2str(t),';']); 
              % Z(end+1) = ztmp;% + conti/8;
            %end
        end
     end
  end
end
figure
scatter_color_map = 'bk';
for t = 1:2
    eval(['scatter(X',int2str(t),', Y',int2str(t),' , 5, '' ',scatter_color_map(t),' ''); ']);
    hold on
end
small_move=quiver(x,y,u,v,'color','green');
hold on
Lar_move=quiver(lar_x,lar_y,lar_u,lar_v,'color','red');

%imwrite(BWdfill,sprintf('Copy_of_StrainEnergy3D_SD_2015_12_31/CELL/Cell_%d_%d.tif', t, z))



% threshold = 170;
% ll=4;
% se90=strel('line',ll,90);
% se0=strel('line',ll,0);
% t1_20_40 = imread('MAX_t1_20_40.tif');
% im_fil2_t1=medfilt2(t1_20_40); %necessary
% im_threshold_t1 = im_fil2_t1 - threshold*ones(1002,1004,'uint16'); %necessary
% BW1_t1 = edge(im_threshold_t1, 'sobel');
% BWdfill_t1= imdilate(BW1_t1,[se90 se0]);
% t2_20_40 = imread('MAX_t2_20_40.tif');
% im_fil2_t2=medfilt2(t2_20_40); %necessary
% im_threshold_t2 = im_fil2_t2 - threshold*ones(1002,1004,'uint16'); %necessary
% BW1_t2 = edge(im_threshold_t2, 'sobel');
% BWdfill_t2= imdilate(BW1_t2,[se90 se0]);
% figure,
% imshowpair(t1_20_40,t2_20_40,'diff');