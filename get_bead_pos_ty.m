function A = get_bead_pos_ty(timestep, bead_num)

%source = 'bead_tnxyz';

%existStr=sprintf('exist(''%s'', ''var'')==1', source);
%isVarExist=evalin('base', existStr);
%if ~isVarExist
%   get_beadtnxyz_han; 
%end

bead_tnxyz = evalin('base','bead_tnxyz');
%load('Copy_of_StrainEnergy3D_SD_2016-01-01/bead_tnxyz.mat','bead_tnxyz');

n_row_idx = (bead_tnxyz(:,2) == bead_num);
n_filtered = bead_tnxyz(n_row_idx,:);
t_row_idx = (n_filtered(:,1) == timestep);
nt_filtered = n_filtered(t_row_idx,:);

x_cor = nt_filtered(3);
y_cor = nt_filtered(4);
z_cor = nt_filtered(5);

x_cor = x_cor/6;
y_cor = y_cor/6;
A=[x_cor y_cor z_cor];

end
