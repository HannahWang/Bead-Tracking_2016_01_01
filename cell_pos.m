
    for ztmp = 7:20
   I = imread(sprintf('Cell_1_%d.tif', ztmp));
     for i = 1:1002
          for j = 1:1004
           if I(i,j) == 1
             cell_XYZ(i,j,ztmp) = 1;
             %X(end+1) = i/6;
             %Y(end+1) = j/6;
             %Z(end+1) = ztmp;
           end
          end
     end
    end
    