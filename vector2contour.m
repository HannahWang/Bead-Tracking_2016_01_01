% function [x,y,z,u,v,w,count]=bemove_col4(t,z_min,z_max,sampling)
t=1;
z_min=240;
z_max=270;
sampling=10;

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

x = [];
y = [];
z = [];
u = [];
v = [];
w = [];

count = [0 0 0 0 0 0 0];

length=linspace(1,7,7);
colormap1=jet(7);
cmap_ref=[transpose(length) jet(7)];
cmap=[];
headWidth = 8;
headLength = 8;
LineLength = 0.08;


figure,


hold on

% for contour map
height=zeros(168,168);

for x_idx=1:168
    
    for y_idx=1:168
        unit_idx=0;
        sum_idx=0;
        for i=1:number
            
            while abs(x1(i)-x_idx)<=1 && abs(y1(i)-y_idx)<=1 && abs(z1(i)-(z_max-z_min/2))<= z_max-z_min/2
                unit_idx=sqrt(u1(i)^2+v1(i)^2+w1(i)^2);
                sum_idx=[unit_idx unit_idx];
                
                break
            end
        end
        height(x_idx,y_idx)=mean(sum_idx);
        
    end
    
end
 x_coor=linspace(1,168,168);
 y_coor=linspace(1,168,168);
 contourf(x_coor,y_coor,height);
 
% for vector map/quiver3
displacement=[];

for i=1:number
    x=[x x1(i)];
    y=[y y1(i)];
    z=[z z1(i)];
    u=[u u1(i)];
    v=[v v1(i)];
    w=[w w1(i)];
    
    
    displacement(i)=sqrt(u1(i)^2+v1(i)^2+w1(i)^2);
    ll=floor(displacement(i))+1;
    
        if ll>7
            count(7)=count(7)+1;
            cmap(i,:)=cmap_ref(7,[2 3 4]);
        else
             count(ll)=count(ll)+1;
             cmap(i,:)=cmap_ref(ll,[2 3 4]);      
        end   
    small_move=quiver3(x(i)*6,y(i)*6,z(i)-z_min,u(i)*6,v(i)*6,w(i),'color',cmap(i,:),'MaxHeadSize',70,'AutoScaleFactor',0.89,'AutoScale','on');
    alpha = 0.33; 
    beta = 0.23; 
    
end

        
            set(gca,'color',[0 0 0]);
            xlabel('x axis','fontsize',14);
            ylabel('y axis','fontsize',14);
            zlabel('z axis','fontsize',14);
            axis([0 168 0 168]);
            caxis([0 7]);
            colorbar('Ticks',[0 1 2 3 4 5 6 7],'TickLabels',{'0','1','2','3','4','5','6','7'},'fontsize',14)
            ylabel(colorbar,'Displacement(um)','fontsize',14);
            c=colorbar;
            ax=findobj(gcf,'type','axes');
            c1=colorbar('peer',ax);
            cbfreeze(c1);


%set(small_move,'CurrentAxes',ax)

% c=colorbar;
% colorbar('peer',ax)
% cbfreeze(colorbar)

% set(c,'Tag','colorbars1','UserData',struct('associatedAxes',small_move));
% h=get(c,'UserData');
% set(h.associatedAxes,'Selected','on');
% h=findobj(gcf,'Tag','colorbars1');
% set(h,'Selected', 'on');

% imshow(max_t)
% max_t.colorbar('off')



%max_t=imagesc([resol resol*1002],[resol resol*1004],max_im);

for k=z_min:z_max
I=imread(sprintf('Copy_of_StrainEnergy3D_SD_2016-01-01/BFFRAME/BFframe_t%06i_%04i.tif',t,k));    
comp_im(:,:,k-z_min+1)=I; 
end

%z projection
% max_im=5*max(comp_im,[],3);
% resol=1/6;
% max_t=imresize(max_im,resol);
% im_t=imshow(max_t);


%y projection

%x-z projection 
% max_2=max(comp_im,[],2);
% max_2t = squeeze(max_2);
% max_2t=15*max_2t;
% resol=1;
% max_t2=imresize(max_2t,resol);
% 
% imshow(max_t2)
% 
% imagesc([resol resol*1002],[resol resol*1004],max_im);

%end
