function [x,y,z,u,v,w,count]=bemove_col4(timestep,start_z,max_z,sampling)


load('bead_tnxyz.mat','bead_tnxyz');
bead_tnxyz = evalin('base','bead_tnxyz');
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
jet_div=jet(7);
cmap_ref=[transpose(length) jet(7)];
cmap=[];
headWidth = 8;
headLength = 8;
LineLength = 0.08;

figure,
hold on
for i=1:number
    x=[x x1(i)];
    y=[y y1(i)];
    z=[z z1(i)];
    u=[u u1(i)];
    v=[v v1(i)];
    w=[w w1(i)];
    
    ll=floor(sqrt(u1(i)^2+v1(i)^2+w1(i)^2))+1; 
   
        if ll>7
            count(7)=count(7)+1;
            cmap(i,:)=cmap_ref(7,[2 3 4]);
        else
             count(ll)=count(ll)+1;
             cmap(i,:)=cmap_ref(ll,[2 3 4]);      
        end
        
    small_move=quiver3(x(i),y(i),z(i),u(i),v(i),w(i),'color',cmap(i,:),'MaxHeadSize',70,'AutoScaleFactor',0.89,'AutoScale','off');
    alpha = 0.33; 
    beta = 0.23; 
    
end

set(gca,'color',[0 0 0]);
xlabel('x axis','fontsize',14);
ylabel('y axis','fontsize',14);
zlabel('z axis','fontsize',14);
axis([0 168 0 168]);


c=colorbar;
caxis([0 7]);
colorbar('Ticks',[0 1 2 3 4 5 6 7],'TickLabels',{'0','1','2','3','4','5','6','7'},'fontsize',14)
ylabel(colorbar,'Displacement(um)','fontsize',14);

end
