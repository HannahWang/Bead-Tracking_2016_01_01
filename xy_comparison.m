%xy_comparison(min_z,max_z)
%comparison bead position of different z

bead_tnxyz=evalin('base', 'bead_tnxyz');
min_z=3;
max_z=5;
z_layer_num=max_z-min_z+1;

for z=min_z:max_z

timestep=1;
count=z-min_z+1; 
n_row_idx = (bead_tnxyz(:,5) == z);
n1_filtered = bead_tnxyz(n_row_idx,:);
t_row_idx = (n1_filtered(:,1) == 1);
nt1_filtered = n1_filtered(t_row_idx,:); 
num_idx=bead_tnxyz(n_row_idx,2);

timestep=2;
count2=z-min_z+1+z_layer_num;
A=bead_tnxyz(t_row_idx ,:);
t2_row_idx = (A(:,1) == 2);
nt2_filtered = A(t2_row_idx,:); 
%take out x y z position
eval(['x',int2str(count),'= transpose(nt1_filtered(:,3));'])
eval(['y',int2str(count),'= transpose(nt1_filtered(:,4));'])
eval(['z',int2str(count),'= transpose(nt1_filtered(:,5));'])
%plot
eval(['subplot(2,',int2str(z_layer_num),',',int2str(count),'), scatter(x',int2str(count),' ,y',int2str(count),');'])
title(sprintf('z%d t1',z));

% for timestep =2
eval(['x',int2str(count2),'= transpose(nt2_filtered(:,3));'])
eval(['y',int2str(count2),'= transpose(nt2_filtered(:,4));'])
eval(['z',int2str(count2),'= transpose(nt2_filtered(:,5));'])

% plot
eval(['subplot(2,',int2str(z_layer_num),',',int2str(count2),'), scatter(x',int2str(count2),' ,y',int2str(count2),');'])
title(sprintf('z%d t2',z));


end

