function [x_padu,y_padu,z_padu,u_padu,v_padu,w_padu,count]=bemove_col4(t,z_min,z_max,sampling,padu,i_layer);
% t=1;
% z_min=250;
% z_max=300;
% sampling=30;
load('cell_1.mat', 'out_matrix');
cell = out_matrix;
load('bead_tnxyz.mat','bead_tnxyz');
bead_tnxyz = evalin('base','bead_tnxyz');
%find start_num
n_row_idx = (bead_tnxyz(:,5) == z_max);
n_filtered = bead_tnxyz(n_row_idx,:);
n_filt_size = size(n_filtered);
while n_filt_size(1) == 0
    n_row_idx = (bead_tnxyz(:,5) == z_max-1);
    n_filtered = bead_tnxyz(n_row_idx,:);
end
start_num = min(n_filtered(:,2));
%find max_bead_num
n_row_idx = (bead_tnxyz(:,5) == z_min);
n_filtered = bead_tnxyz(n_row_idx,:);
n_filt_size = size(n_filtered);
while n_filt_size(1) == 0
    n_row_idx = (bead_tnxyz(:,5) == z_min+1);
    n_filtered = bead_tnxyz(n_row_idx,:);
end
max_bead_num = max(n_filtered(:,2));
% position 
number=((max_bead_num-start_num+1)-mod((max_bead_num-start_num+1),sampling))/sampling; 
Atemp=zeros(number,3);
Atemp2=zeros(number,3);
bead_num=1;
index_i=1;
for i=start_num:sampling:max_bead_num
%original
a1=get_bead_pos_ty(t, i);
Atemp(index_i, :) = a1;
%next
a2=get_bead_pos_ty(t+1,i);
Atemp2(index_i, :) = a2;

index_i=index_i+1;
end

ori=Atemp.';

next=Atemp2.';

%displacement
dis=next-ori;

dis2=dis;
x1=ori(1,:);
y1=ori(2,:);
z1=ori(3,:);
u1=dis2(1,:);
v1=dis2(2,:);
w1=dis2(3,:);

    x=0;
    y=0;
    z=0;
    u=0;
    v=0;
    w=0;

x_padu = [];
y_padu = [];
z_padu= [];
u_padu= [];
v_padu= [];
w_padu= [];

count = [0 0 0 0 0 0 0];

length=linspace(1,7,7);
colormap1=jet(7);
cmap_ref=[transpose(length) jet(7)];
cmap=[];
headWidth = 8;
headLength = 8;
LineLength = 0.08;


% figure,
% % draw xy vector map
% displacement=[];
% hold on
% 

% %cell = padarray(cell, [6*padu, 6*padu, padu]);
% for i=1:number
%     x_min_pad = round(x1(i)*6)-6*padu;
%     if x_min_pad < 1
%        x_min_pad = 1; 
%     end
%     x_max_pad = round(x1(i)*6)+6*padu;
%     if x_max_pad > 1002
%        x_max_pad = 1002; 
%     end
%     y_min_pad = round(y1(i)*6)-6*padu;
%     if y_min_pad < 1
%        y_min_pad = 1; 
%     end
%     y_max_pad = round(y1(i)*6)+6*padu;
%     if y_max_pad > 1004
%        y_max_pad = 1004; 
%     end
%     z_min_pad = round(z1(i))-padu;
%     if z_min_pad < 1
%        z_min_pad = 1; 
%     end
%     z_max_pad = round(z1(i))+padu;
%     if z_max_pad > 397
%        z_max_pad = 397; 
%     end
%     cumu_cellx = sum(sum(sum(cell(x_min_pad:x_max_pad, y_min_pad:y_max_pad, z_min_pad:z_max_pad))));
%     
%         x=[x x1(i)];
%         y=[y y1(i)];
%         z=[z z1(i)];
%         u=[u u1(i)];
%         v=[v v1(i)];
%         w=[w w1(i)];
%     if cumu_cellx > 0
%         x_padu=[x_padu x1(i)];
%         y_padu=[y_padu y1(i)];
%         z_padu=[z_padu z1(i)];
%         u_padu=[u_padu u1(i)];
%         v_padu=[v_padu v1(i)];
%         w_padu=[w_padu w1(i)];
%         
%        displacement(i)=sqrt(u1(i)^2+v1(i)^2+w1(i)^2);
%         ll=floor(5*displacement(i))+1;
% 
%             if ll>7
%                 count(7)=count(7)+1;
%                 cmap(i,:)=cmap_ref(7,[2 3 4]);
%             else
%                  count(ll)=count(ll)+1;
%                  cmap(i,:)=cmap_ref(ll,[2 3 4]);      
%             end
% 
%         %move_xy=quiver(x(i),y(i),u(i),v(i),'color',cmap(i,:),'MaxHeadSize',70,'AutoScaleFactor',0.89,'AutoScale','off');
% 
%         %small_move=quiver3(x(i),y(i),z(i)-z_min,u(i),v(i),w(i),'color',cmap(i,:),'MaxHeadSize',100,'AutoScaleFactor',0.89,'AutoScale','off');
%         alpha = 0.33; 
%         alpha = 0.33; 
%         beta = 0.23; 
%     end
% end
% 
%         
%             set(gca,'color',[0 0 0]);
%             title('0101_xy');
%             xlabel('x axis','fontsize',14);
%             ylabel('y axis','fontsize',14);
%            
%             axis([0 168 0 168]);
%             caxis([0 7]);
%             c=colorbar;
%             set(c,'YTick',[1,2,3,4,5,6,7]);
           
 
 figure,
% draw yz vector map
displacement=[];
hold on

for i=1:number
    x_min_pad = round(x1(i)*6)-6*padu;
    if x_min_pad < 1
       x_min_pad = 1; 
    end
    x_max_pad = round(x1(i)*6)+6*padu;
    if x_max_pad > 1002
       x_max_pad = 1002; 
    end
    y_min_pad = round(y1(i)*6)-6*padu;
    if y_min_pad < 1
       y_min_pad = 1; 
    end
    y_max_pad = round(y1(i)*6)+6*padu;
    if y_max_pad > 1004
       y_max_pad = 1004; 
    end
    z_min_pad = round(z1(i))-padu;
    if z_min_pad < 1
       z_min_pad = 1; 
    end
    z_max_pad = round(z1(i))+padu;
    if z_max_pad > 397
       z_max_pad = 397; 
    end
    cumu_cellx = sum(sum(sum(cell(x_min_pad:x_max_pad, y_min_pad:y_max_pad, z_min_pad:z_max_pad))));
    
    x=[x x1(i)];
    y=[y y1(i)];
    z=[z z1(i)];
    u=[u u1(i)];
    v=[v v1(i)];
    w=[w w1(i)];
    if cumu_cellx > 0
        
        x_padu=[x_padu x1(i)];
        y_padu=[y_padu y1(i)];
        z_padu=[z_padu z1(i)];
        u_padu=[u_padu u1(i)];
        v_padu=[v_padu v1(i)];
        w_padu=[w_padu w1(i)];
        
   displacement(i)=sqrt(u1(i)^2+v1(i)^2+w1(i)^2);
    ll=floor(5*displacement(i))+1;
   
        if ll>7
            count(7)=count(7)+1;
            cmap(i,:)=cmap_ref(7,[2 3 4]);
        else
             count(ll)=count(ll)+1;
             cmap(i,:)=cmap_ref(ll,[2 3 4]);      
        end
    
        
    if i_layer==1    
    move_yz=quiver(y(i),z(i),v(i),w(i),'color',cmap(i,:),'MaxHeadSize',70,'AutoScaleFactor',0.89,'AutoScale','off');
    end
    %small_move=quiver3(x(i),y(i),z(i)-z_min,u(i),v(i),w(i),'color',cmap(i,:),'MaxHeadSize',100,'AutoScaleFactor',0.89,'AutoScale','off');
    alpha = 0.33; 
    alpha = 0.33; 
    beta = 0.23; 
    end
end

            set(gca,'color',[0 0 0]);
            title('0101_yz');
            xlabel('y axis','fontsize',14);
            ylabel('z axis','fontsize',14);
            axis([0 168 0 399]);
            caxis([0 7]);
            %colorbar('Ticks',[0 1 2 3 4 5 6 7],'Ticklabels',{'0','1','2','3','4','5','6','7'},'fontsize',14)
            
%             c=colorbar;
%             set(c,'YTick',[1,2,3,4,5,6,7]);
% %             ax=findobj(gcf,'type','axes');
% %             c1=colorbar('peer',ax);
% %             cbfreeze(c1);
% 
%  figure,
% % draw xz vector map
% displacement=[];
% hold on
% 
% for i=1:number
%     x_min_pad = round(x1(i)*6)-6*padu;
%     if x_min_pad < 1
%        x_min_pad = 1; 
%     end
%     x_max_pad = round(x1(i)*6)+6*padu;
%     if x_max_pad > 1002
%        x_max_pad = 1002; 
%     end
%     y_min_pad = round(y1(i)*6)-6*padu;
%     if y_min_pad < 1
%        y_min_pad = 1; 
%     end
%     y_max_pad = round(y1(i)*6)+6*padu;
%     if y_max_pad > 1004
%        y_max_pad = 1004; 
%     end
%     z_min_pad = round(z1(i))-padu;
%     if z_min_pad < 1
%        z_min_pad = 1; 
%     end
%     z_max_pad = round(z1(i))+padu;
%     if z_max_pad > 397
%        z_max_pad = 397; 
%     end
%     cumu_cellx = sum(sum(sum(cell(x_min_pad:x_max_pad, y_min_pad:y_max_pad, z_min_pad:z_max_pad))));
%     
%         x=[x x1(i)];
%         y=[y y1(i)];
%         z=[z z1(i)];
%         u=[u u1(i)];
%         v=[v v1(i)];
%         w=[w w1(i)];
%     if cumu_cellx > 0
%         
%         x_padu=[x_padu x1(i)];
%         y_padu=[y_padu y1(i)];
%         z_padu=[z_padu z1(i)];
%         u_padu=[u_padu u1(i)];
%         v_padu=[v_padu v1(i)];
%         w_padu=[w_padu w1(i)];
%         
%         displacement(i)=sqrt(u1(i)^2+v1(i)^2+w1(i)^2);
%         ll=floor(5*displacement(i))+1;
% 
%             if ll>7
%                 count(7)=count(7)+1;
%                 cmap(i,:)=cmap_ref(7,[2 3 4]);
%             else
%                  count(ll)=count(ll)+1;
%                  cmap(i,:)=cmap_ref(ll,[2 3 4]);      
%             end
% 
%         %move_xz=quiver(x(i),z(i),u(i),w(i),'color',cmap(i,:),'MaxHeadSize',70,'AutoScaleFactor',0.89,'AutoScale','off');
% 
%         %small_move=quiver3(x(i),y(i),z(i)-z_min,u(i),v(i),w(i),'color',cmap(i,:),'MaxHeadSize',100,'AutoScaleFactor',0.89,'AutoScale','off');
%         alpha = 0.33; 
%         alpha = 0.33; 
%         beta = 0.23; 
%     end
% end
% 
%             set(gca,'color',[0 0 0]);
%             title('0101_xz');
%             xlabel('x axis','fontsize',14);
%             ylabel('z axis','fontsize',14);
%             axis([0 168 0 399]);
%             caxis([0 7]);
%             %colorbar('Ticks',[0 1 2 3 4 5 6 7],'Ticklabels',{'0','1','2','3','4','5','6','7'},'fontsize',14)
%             
%             c=colorbar;
%             set(c,'YTick',[1,2,3,4,5,6,7]);
% %set(small_move,'CurrentAxes',ax)
% 
% % c=colorbar;
% % colorbar('peer',ax)
% % cbfreeze(colorbar)
% 
% % set(c,'Tag','colorbars1','UserData',struct('associatedAxes',small_move));
% % h=get(c,'UserData');
% % set(h.associatedAxes,'Selected','on');
% % h=findobj(gcf,'Tag','colorbars1');
% % set(h,'Selected', 'on');
% 
% % imshow(max_t)
% % max_t.colorbar('off')
% 
% 
% 
% %max_t=imagesc([resol resol*1002],[resol resol*1004],max_im);
% 
% % for k=z_min:z_max
% % I=imread(sprintf('Copy_of_StrainEnergy3D_SD_2016-01-01/BFFRAME/BFframe_t%06i_%04i.tif',t,k));    
% % comp_im(:,:,k-z_min+1)=I; 
% % end
% 
% %z projection
% % max_im=5*max(comp_im,[],3);
% % resol=1/6;
% % max_t=imresize(max_im,resol);
% % im_t=imshow(max_t);
% 
% 
% %y projection
% 
% %x-z projection 
% % max_2=max(comp_im,[],2);
% % max_2t = squeeze(max_2);
% % max_2t=15*max_2t;
% % resol=1;
% % max_t2=imresize(max_2t,resol);
% % 
% % imshow(max_t2)
% % 
% % imagesc([resol resol*1002],[resol resol*1004],max_im);
% 
end
