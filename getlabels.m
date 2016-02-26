function getlabels(t,z)
signal = imread(sprintf('CELL/Cell_%d_%d.tif',t,z));
raw = imread(getBFfileName(z,t));
square_width = 17;
padding = (square_width-1)/2;
sampling = 3;
nosig_sampling = 0;
for i = 1+padding:sampling:1002-padding
   for j = 1+padding:sampling:1004-padding
       now = raw(i-padding:i+padding,j-padding:j+padding);
       if signal(i,j)
           %save(sprintf('Bead-Tracking/LABELS/with-signal/t%d_z%d_%04i_%04i.mat',t,z,i,j), 'now');
       else
           nosig_sampling = nosig_sampling + 1;
           if nosig_sampling == 15
                %save(sprintf('Bead-Tracking/LABELS/no-signal/t%d_z%d_%04i_%04i.mat',t,z,i,j), 'now');
                nosig_sampling = 0;
           end
       end
       
   end
end

end