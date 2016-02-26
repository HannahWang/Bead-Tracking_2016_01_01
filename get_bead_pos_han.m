function [x, y, z] = get_bead_pos(timestep, bead_num)

source = 'bead_tnxyz'

existStr=sprintf('exist(''%s'', ''var'')==1', source);
isVarExist=evalin('base', existStr);
if ~isVarExist
   get_beadnxyz_han; 
end
bead_tnxyz = evalin('base','bead_tnxyz');

n_row_idx = (bead_tnxyz(:,2) == bead_num);
n_filtered = bead_tnxyz(n_row_idx,:);

t_row_idx = (n_filtered(:,1) == timestep);
nt_filtered = n_filtered(t_row_idx,:);

x = nt_filtered(3);
y = nt_filtered(4);
z = nt_filtered(5);

end