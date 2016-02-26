fin = fopen('bead_pos.txt','r');
tline = fgets(fin);
bead_tnxyz = zeros(0,5);
while ischar(tline)
    element = str2num(tline);
    bead_tnxyz = [bead_tnxyz; element(1:5)];
    tline = fgets(fin);
end
save('bead_tnxyz.mat','bead_tnxyz');


