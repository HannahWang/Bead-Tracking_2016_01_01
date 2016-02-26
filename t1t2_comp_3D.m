start_z=20;
max_z=40;
for t=1:2

    BG_matrix = imread(sprintf('CELL/Cell_%d_0.tif',t));
    for k=start_z:max_z
       BG_matrix(:,:,k+1) = imread(sprintf('CELL/Cell_%d_%d.tif', t, k));
    end
    %BG_matrix = padarray(BG_matrix,[6,6,2],0);

%erase1um = 1;
%erase1umij = erase1um * 6;
%for k=start_z+1+1:max_z+1+1
%    for i=1:1002+6-erase1umij*3+1
%        for j=1:1004+6-erase1umij*3+1
%           if (sum(sum(sum(BG_matrix(i:i+erase1umij*3-1,j:j+erase1umij*3-1,k-erase1um:k+erase1um))))-sum(sum(sum(BG_matrix(i+erase1umij:i+erase1umij*2-1,j+erase1umij:j+erase1umij*2-1,k-erase1um+1:k+erase1um-1)))))==0
%               if sum(sum(sum(BG_matrix(i:i+erase1umij*3-1,j:j+erase1umij*3-1,k-erase1um:k+erase1um)))) ~= 0
%                   [k i j]
%                   BG_matrix(i+erase1umij:i+erase1umij*2-1,j+erase1umij:j+erase1umij*2-1,k) = 0;
%               end             
%           end
%       end
%   end
%end


    eval(['X',int2str(t),'=[];']);
    eval(['Y',int2str(t),'=[];']);
    eval(['Z',int2str(t),'=[];']);
for ztmp=start_z:max_z
    for i = 1:1002
        for j = 1:1004
            if BG_matrix(i,j,ztmp+1) == 1
               eval(['X',int2str(t),'(end+1) = i/6;']);
               eval(['Y',int2str(t),'(end+1) = j/6;']); 
               eval(['Z',int2str(t),'(end+1) = ztmp;']); 

            end
        end
    end
end
end
% X = [];
% Y = [];
% Z = [];
% for ztmp = start_z:max_z
%   I = BG_matrix(:,:,ztmp+1);
%   for i = 1:1002
%      for j = 1:1004
%         if I(i,j) == 1
%             %for conti = 0:7
%                X(end+1) = i/6;
%                Y(end+1) = j/6;
%                Z(end+1) = ztmp;% + conti/8;
%             %end
%         end
%      end
%   end
% end

%des_color = [1 0 0;];

% plot vector map
figure
scatter_color_map = 'bk';
for t = 1:2
    eval(['scatter3(X',int2str(t),', Y',int2str(t),', Z',int2str(t),' , 5, '' ',scatter_color_map(t),' ''); ']);
    hold on
end
small_move=quiver3(x,y,z,u,v,w,'color','green');
hold on
Lar_move=quiver3(lar_x,lar_y,lar_z,lar_u,lar_v,lar_w,'color','red');
%hold on
%scatter3(X,Y,Z,5,'.','blue');
