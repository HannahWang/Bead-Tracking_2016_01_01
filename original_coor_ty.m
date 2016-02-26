function [x,y,z,u,v,w,lar_x,lar_y,lar_z,lar_u,lar_v,lar_w,count]=original_coor_ty(timestep,start_z,max_z,sampling)

%bead_tnxyz = evalin('base','bead_tnxyz');
load('Copy_of_StrainEnergy3D_SD_2016-01-01/bead_tnxyz.mat','bead_tnxyz');

%find start_num
n_row_idx = (bead_tnxyz(:,5) == max_z);
n_filtered = bead_tnxyz(n_row_idx,:);
n_filt_size = size(n_filtered);
while n_filt_size(1) == 0
    n_row_idx = (bead_tnxyz(:,5) == max_z-1);
    n_filtered = bead_tnxyz(n_row_idx,:);
end
start_num = min(n_filtered(:,2));
%find max_bead_num
n_row_idx = (bead_tnxyz(:,5) == start_z);
n_filtered = bead_tnxyz(n_row_idx,:);
n_filt_size = size(n_filtered);
while n_filt_size(1) == 0
    n_row_idx = (bead_tnxyz(:,5) == start_z+1);
    n_filtered = bead_tnxyz(n_row_idx,:);
end
max_bead_num = max(n_filtered(:,2));
% position 
number=((max_bead_num-start_num+1)-mod((max_bead_num-start_num+1),sampling))/sampling; 
Atemp=zeros(number,3);
Atemp2=zeros(number,3);
bead_num=1;
index_i=1;
for i=start_num:sampling:max_bead_num+1
%original
a1=get_bead_pos_ty(timestep, i);
Atemp(index_i, :) = a1;
%next
a2=get_bead_pos_ty(timestep+1,i);
Atemp2(index_i, :) = a2;

index_i=index_i+1;

end

ori=Atemp.';

next=Atemp2.';

%displacement
dis=next-ori;

dis2=dis;
x1=ori(1,:)/6;
y1=ori(2,:)/6;
z1=ori(3,:);
u1=dis2(1,:)/6;
v1=dis2(2,:)/6;
w1=dis2(3,:);

x = [];
y = [];
z = [];
u = [];
v = [];
w = [];
lar_x=[];
lar_y=[];
lar_z=[];

lar_u=[];
lar_v=[];
lar_w=[];
count = [0 0 0 0];
for i=1:number
    if  u1(i)^2+v1(i)^2+w1(i)^2 > 1
        lar_x=[lar_x x1(i)];
        lar_y=[lar_y y1(i)];
        lar_z=[lar_z z1(i)];
        lar_u=[lar_u u1(i)];
        lar_v=[lar_v v1(i)];
        lar_w=[lar_w w1(i)]; 
        if u1(i)^2+v1(i)^2+w1(i)^2 > 25
            count(4) = count(4) + 1;
        elseif u1(i)^2+v1(i)^2+w1(i)^2 > 4
            count(3) = count(3) + 1;
        else
            count(2) = count(2) + 1;
        end
        
    else
        x=[x x1(i)]; 
        y=[y y1(i)];
        z=[z z1(i)];
        u=[u u1(i)];
        v=[v v1(i)];
        w=[w w1(i)]; 
        count(1) = count(1) + 1;
    end
end

% for t=1:2
% 
%     BG_matrix = imread(sprintf('CELL/Cell_%d_0.tif',t));
%     for k=start_z:max_z
%        BG_matrix(:,:,k+1) = imread(sprintf('CELL/Cell_%d_%d.tif', t, k));
%     end
%     %BG_matrix = padarray(BG_matrix,[6,6,2],0);
% 
% %erase1um = 1;
% %erase1umij = erase1um * 6;
% %for k=start_z+1+1:max_z+1+1
% %    for i=1:1002+6-erase1umij*3+1
% %        for j=1:1004+6-erase1umij*3+1
% %           if (sum(sum(sum(BG_matrix(i:i+erase1umij*3-1,j:j+erase1umij*3-1,k-erase1um:k+erase1um))))-sum(sum(sum(BG_matrix(i+erase1umij:i+erase1umij*2-1,j+erase1umij:j+erase1umij*2-1,k-erase1um+1:k+erase1um-1)))))==0
% %               if sum(sum(sum(BG_matrix(i:i+erase1umij*3-1,j:j+erase1umij*3-1,k-erase1um:k+erase1um)))) ~= 0
% %                   [k i j]
% %                   BG_matrix(i+erase1umij:i+erase1umij*2-1,j+erase1umij:j+erase1umij*2-1,k) = 0;
% %               end             
% %           end
% %       end
% %   end
% %end
% 
% 
%     eval(['X',int2str(t),'=[];']);
%     eval(['Y',int2str(t),'=[];']);
%     eval(['Z',int2str(t),'=[];']);
% for ztmp=start_z:max_z
%     for i = 1:1002
%         for j = 1:1004
%             if BG_matrix(i,j,ztmp) == 1
%                eval(['X',int2str(t),'(end+1) = i/6;']);
%                eval(['Y',int2str(t),'(end+1) = j/6;']); 
%                eval(['Z',int2str(t),'(end+1) = ztmp;']); 
% 
%             end
%         end
%     end
% end
% end
% % X = [];
% % Y = [];
% % Z = [];
% % for ztmp = start_z:max_z
% %   I = BG_matrix(:,:,ztmp+1);
% %   for i = 1:1002
% %      for j = 1:1004
% %         if I(i,j) == 1
% %             %for conti = 0:7
% %                X(end+1) = i/6;
% %                Y(end+1) = j/6;
% %                Z(end+1) = ztmp;% + conti/8;
% %             %end
% %         end
% %      end
% %   end
% % end

%des_color = [1 0 0;];

% plot vector map
% figure
% scatter_color_map = 'bk';
% for t = 1:2
%     eval(['scatter3(X',int2str(t),', Y',int2str(t),', Z',int2str(t),' , 5, '' ',scatter_color_map(t),' ''); ']);
%     hold on
% end
figure,
small_move=quiver3(x,y,z,u,v,w,'color','green');
hold on
Lar_move=quiver3(lar_x,lar_y,lar_z,lar_u,lar_v,lar_w,'color','red');
%hold on
%scatter3(X,Y,Z,5,'.','blue');

end
