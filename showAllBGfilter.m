t = 1;
threshold = 320;
eraseSs = [1, 2, 5];
compare = 2;
[x_size, y_size] = size(imread(getBFfileName(0,t)));
last = zeros(x_size, y_size);
now = zeros(x_size, y_size);
next = BGfilter(t, 0, threshold, eraseSs);
for z = 0:399
    last = now;
    now  = next;
    if z < 399
        next = BGfilter(t, z+1, threshold, eraseSs);
    else
        next = zeros(x_size, y_size);
    end
    
    for x=1:x_size
        for y=1:y_size
            if now(x,y)
                minx = max(1, x-compare);
                maxx = min(x_size, x+compare);
                miny = max(1, y-compare);
                maxy = min(y_size, y+compare);
                if (sum(last(minx:maxx,miny:maxy))+sum(next(minx:maxx,miny:maxy))) == 0
                    now(x,y)=0;
                end
            end
        end
    end
    out_matrix(:,:,z+1) = now; 
    imshowpair(now,imread(getBFfileName(z,t)),'montage');
    title(['fig',num2str(z)],'color','b');
end
save(sprintf('cell_%d.mat',t), 'out_matrix');
