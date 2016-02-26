function browseStack(targetdir)
% browseStack(targetdir, bp)
% Navigate images of the stacks in targetdir, defaulting to pwd.
% 
% Pressing 'n' or 'down arrow' displays next z image 
%          'p' or 'up arrow' displays the previous z image
%          'N' or 'right arrow' displays next t image
%          'P' or 'left arrow' displays previous t image
%          'B' toggles between brightfield and fluorecent images
%          'e' or 'home/pos 1' jump to beginning
%          'E' or 'end' jump to end
%          'C' jump to center
%          'f' choose a bead to follow by mouse input, press again to stop
%          'F' choose a bead to follow by number, press again to stop
%          'g' get bead numberB
%          'o' for outlining the cell and writing cell_shape_XXXXXX.mat where XXXXXX corresponds to current time step
%          '1..0' chooses increment in the range of 1 to 10
%          'O' opens bead position file and plots bead positions. 
%          '+' increases gamma for gamma correction
%          '-' decreases gammaO
%          '=' resets gamma
%          'h' add top (highest) point to height profile of the cell
%          'H' delete top point from height profile
%          'm' add mid/center point (cell body) to height profile
%          'M' delete mid point from height profile
%          'l' add bottom (lowest) point to height profile
%          'L' delete bottom point
%          'w' write height profile to height_profile_XXXXXX.mat where XXXXXX corresponds to current time step
%          't' choose cell position for current time and append to list
%          'T' write list of cell positions for all times to cellpos.mat
%          'q' quits.
% -------------------------------------------------------------------------
% This file is part of the method published in
%
% Koch TM, Münster S, Bonakdar N, Butler JP, Fabry B (2012) 3D traction 
% forces in cancer cell invasion. PLoS ONE
%
% If you use any part of it, please cite this paper.
% -------------------------------------------------------------------------
cwd=pwd;
if(nargin>=1)
    cd(targetdir);
end
%
% handles for bead plotting
%
bp=[];
point_handles_good=[];
point_handles_bad=[];
point_handles_above1=[];
point_handles_above2=[];
point_handles_above3=[];
point_handles_below1=[];
point_handles_below2=[];
point_handles_below3=[];
follow_bead_handle=[];
follow_bead=[];
%
% handles for height profile
%
low=[];
mid=[];
high=[];
hprofile=[];
hprofile_handles=[];
hprofile_handles_above=[];
hprofile_handles_below=[];
%
% determine stack size and duration assuming that 
% the stack size is constant for all timesteps and that both are identical
% for BF and FL files.
%
z=0;
while(exist(getFLfileName(z,1),'file'))
    z=z+1;
end
stacksize=z;
t=1;
while(exist(getFLfileName(0,t),'file'))
    t=t+1;
end
duration=t-1;
cellposVsTime=zeros(duration,3);
%
% frame buffer to avoid reloading
%
BFbufferIdxMap=zeros(stacksize,duration,'uint16');
BFbufferIdx=1;
BFbuffer=[];
FLbufferIdxMap=zeros(stacksize,duration,'uint16');
FLbufferIdx=1;
FLbuffer=[];
%
% remaining init
%
current_z=0;
current_t=1;
mode='FL';
inc=1;
FLgamma=1.0;
BFgamma=1.0;
gammaval=1.0;
%
% init figure
%
fh=figure;
clf reset;
colormap(gray(256));
%
fname=getFLfileName(0,1); 
frame=double(imread(fname));
frame=frame./max(frame(:));
frame=uint8(255*frame);
[imgh imgw]=size(frame);
FLbuffer(FLbufferIdx).frame=frame;
FLbuffer(FLbufferIdx).fname=fname;
FLbufferIdxMap(1,1)=FLbufferIdx;
FLbufferIdx=FLbufferIdx+1;
imh=imagesc(frame,[0 255]); 
axis off;
set(gca, 'Position', [0 0 1 1]);
set(gca, 'YDir', 'normal');
th=text(10, 20,strrep(fname,'\','\\'), 'Color', [0 1 0]);
set(fh,'Name',strrep(fname,'\','\\'));
%
set(fh, 'KeyPressFcn', @handlekeypress);
uiwait;

%
% Callback functions
%
function handlekeypress(src, event)    
    if(isempty(event.Character))
        event.Character='ß';
    end
    if(event.Character=='p' || strcmp(event.Key, 'uparrow'))
        current_z=current_z-inc;
        updateimg();
    elseif(event.Character=='n' || strcmp(event.Key, 'downarrow'))
        current_z=current_z+inc;
        updateimg();
    elseif(event.Character=='P' || strcmp(event.Key, 'leftarrow'))
        current_t=current_t-inc;
        updateimg();
    elseif(event.Character=='N' || strcmp(event.Key, 'rightarrow'))
        current_t=current_t+inc;
        updateimg();    
    elseif(event.Character=='e' || strcmp(event.Key, 'home'))
        current_z=0;
        updateimg();    
    elseif(event.Character=='E' || strcmp(event.Key, 'end'))
        current_z=Inf;
        updateimg();  
    elseif(event.Character=='C')
        current_z=round(stacksize/2);
        updateimg();     
    elseif(event.Character=='B')
        if(strcmp(mode,'BF'))
            mode='FL';            
        else
            mode='BF';            
        end
        updateimg();        
    elseif(event.Character=='g')
        [cx cy]=ginput(1);
        bn=find( (bp(:,1)==current_t) & (bp(:,5)>(current_z-0.5)) & (bp(:,5)<=(current_z+0.5)) & (bp(:,3)>(cx-4)) & (bp(:,3)<(cx+4)) & (bp(:,4)>(cy-4)) & (bp(:,4)<(cy+4)));
        if(~isempty(bn))         
            fprintf('Chosen bead # %d.\n', bn(1));
        end
    elseif(event.Character=='f') % choose bead with mouse
        if(isempty(follow_bead))
            [cx cy]=ginput(1);
            bn=find( (bp(:,1)==current_t) & (bp(:,5)>(current_z-0.5)) & (bp(:,5)<=(current_z+0.5)) & (bp(:,3)>(cx-4)) & (bp(:,3)<(cx+4)) & (bp(:,4)>(cy-4)) & (bp(:,4)<(cy+4)));
            if(~isempty(bn)) 
                follow_bead=bp(bn(1),2);
                fprintf('Chosen bead # %d.\n', bn(1));
            end    
        else
            follow_bead=[];
        end
    elseif(event.Character=='F') % choose bead by number
        if(isempty(follow_bead))
            bni=inputdlg('Enter bead number:');
            if(~isempty(bni))
                try
                    bn=find( (bp(:,1)==current_t) & (bp(:,2)==str2num(bni{1})) );
                    if(~isempty(bn))
                        follow_bead=bp(bn(1),2);
                    end
                end
            end
        else
            follow_bead=[];  
        end
    elseif(event.Character=='o')
        outlineCell();     
    elseif(event.Character=='1')
        inc=1;
    elseif(event.Character=='2')
        inc=2;
    elseif(event.Character=='3')
        inc=3;
    elseif(event.Character=='4')
        inc=4;
    elseif(event.Character=='5')
        inc=5;
    elseif(event.Character=='6')
        inc=6;
    elseif(event.Character=='7')
        inc=7;
    elseif(event.Character=='8')
        inc=8;
    elseif(event.Character=='9')
        inc=9;
    elseif(event.Character=='0')
        inc=10;
    elseif(event.Character=='O') % load bead positions
        [fname pname]=uigetfile('*.txt', 'Open bead position file');
        if(fname)
            bp=dlmread(sprintf('%s/%s',pname,fname),'\t');
            updateimg();
        end
    elseif(event.Character=='+') % gamma correction
        changeGamma(0.05);
        updateimg();
    elseif(event.Character=='-') % gamma correction
        changeGamma(-0.05);
        updateimg();
    elseif(event.Character=='=') % gamma correction
        resetGamma();
        updateimg();
    elseif(event.Character=='h') % add highest point to height profile
        [x y, button]=ginput(1); % choose one point with the mouse, right button cancels       
        if(button~=3)              
            high=[x y current_z];
        end
        plotHeightProfile();
    elseif(event.Character=='H') % delete highest point to height profile
        high=[];
        plotHeightProfile();
    elseif(event.Character=='m') % add mid point to height profile
        [x y, button]=ginput(1); % choose one point with the mouse, right button cancels       
        if(button~=3)              
            mid=[x y current_z];
        end
        plotHeightProfile();
    elseif(event.Character=='M') % delete mid point to height profile
        mid=[];
        plotHeightProfile();
    elseif(event.Character=='l') % add lowest point to height profile
        [x y, button]=ginput(1); % choose one point with the mouse, right button cancels       
        if(button~=3)              
            low=[x y current_z];
        end
        plotHeightProfile();
    elseif(event.Character=='L') % delete lowest point to height profile
        low=[];
        plotHeightProfile();    
    elseif(event.Character=='w') % writes height profile to height_profile_XXXXXX.mat
        if(isempty(high) || isempty(mid) || isempty(low))
            errordlg('You need to set three points (high/mid/low) to write height profile.', 'Not all points set', 'modal');
        else
            oname=sprintf('height_profile_%06d.mat',current_t);
            ok=false;
            if(exist(oname,'file')==2)
                cont=questdlg(sprintf('File %s already exists. Overwrite?', oname), 'Overwrite file?', 'Overwrite', 'Cancel', 'Cancel');
                if(isequal(cont,'Overwrite'))
                    ok=true;
                else
                    ok=false;     
                end
            else
                ok=true;
            end
            if(ok)
                save(oname, 'hprofile');
            end
        end
    elseif(event.Character=='t') % add cell pos to list
        [x y, button]=ginput(1); % choose one point with the mouse, right button cancels       
        if(button~=3)              
            cellposVsTime(current_t,:)=[x y current_z-1];
        end
    elseif(event.Character=='T') % writes cellposVsTime to cellpos.mat        
        oname='cellpos.mat';
        ok=false;
        if(exist(oname,'file')==2)
            cont=questdlg(sprintf('File %s already exists. Overwrite?', oname), 'Overwrite file?', 'Overwrite', 'Cancel', 'Cancel');
            if(isequal(cont,'Overwrite'))
                ok=true;
            else
                ok=false;     
            end
        else
            ok=true;
        end
        if(ok)
            save(oname, 'cellposVsTime');
        end        
    elseif(event.Character=='q')
        cd(cwd);
        uiresume;
        close(gcf);
    end    
end
%
function updateimg()
    if(current_t<1)
        current_t=1;
    elseif(current_t>duration)
        current_t=duration;
    end
    if(current_z<1)
        current_z=1;
    elseif(current_z>stacksize)
        current_z=stacksize;
    end
    if(strcmp(mode,'FL'))
        if(FLbufferIdxMap(current_z,current_t)==0)
            fname=getFLfileName(current_z-1, current_t);
            if(exist(fname,'file'))
                frame=double(imread(fname));
            else
                frame=ones(imgh, imgw);
                warning('File %s does not exist.', fname);
            end
            frame=frame./max(frame(:));  
            sl=stretchlim(frame);
            sl=[1 2]'.*sl;
            if(sl(2)>1) 
                sl(2)=1;
            end
            frame=imadjust(frame, sl);
            frame=uint8(255.*frame);
            FLbufferIdxMap(current_z, current_t)=FLbufferIdx;
            FLbuffer(FLbufferIdx).frame=frame;
            FLbuffer(FLbufferIdx).fname=fname;
            FLbufferIdx=FLbufferIdx+1;
        else
            frame=FLbuffer(FLbufferIdxMap(current_z,current_t)).frame;
            fname=FLbuffer(FLbufferIdxMap(current_z,current_t)).fname;
        end           
        gammaval=FLgamma;
    else                        
        if(BFbufferIdxMap(current_z,current_t)==0)
            fname=getBFfileName(current_z-1, current_t);
            if(exist(fname,'file'))
                frame=double(imread(fname));
            else
                frame=ones(imgh, imgw);
                warning('File %s does not exist.', fname);
            end
            frame=frame./max(frame(:));                                       
            frame=uint8(255.*frame);             
            BFbufferIdxMap(current_z, current_t)=BFbufferIdx;
            BFbuffer(BFbufferIdx).frame=frame;
            BFbuffer(BFbufferIdx).fname=fname;
            BFbufferIdx=BFbufferIdx+1;
        else
            frame=BFbuffer(BFbufferIdxMap(current_z,current_t)).frame;
            fname=BFbuffer(BFbufferIdxMap(current_z,current_t)).fname;
        end                  
        gammaval=BFgamma;
    end
    set(imh, 'CData', imadjust(frame,[],[],gammaval));     
    if(~isempty(follow_bead))
        fname=sprintf('%s - bead %d', fname, follow_bead);
    end
    set(th, 'String', strrep(fname,'\','\\'));
    set(fh,'Name',strrep(fname,'\','\\'));
    if(~isempty(bp))
        plotbeads();
    end   
    if(~isempty(hprofile))
        plotHeightProfile();
    end
end  
%
function plotbeads()
    if(~isempty(point_handles_good))
        delete(point_handles_good);
        point_handles_good=[];
    end
    if(~isempty(point_handles_bad))
        delete(point_handles_bad);
        point_handles_bad=[];
    end
    if(~isempty(follow_bead_handle))
        delete(follow_bead_handle);
        follow_bead_handle=[];
    end
    if(~isempty(point_handles_above1))
        delete(point_handles_above1);
        point_handles_above1=[];
    end
    if(~isempty(point_handles_above2))
        delete(point_handles_above2);
        point_handles_above2=[];
    end
    if(~isempty(point_handles_above3))
        delete(point_handles_above3);
        point_handles_above3=[];
    end
    if(~isempty(point_handles_below1))
        delete(point_handles_below1);
        point_handles_below1=[];
    end
    if(~isempty(point_handles_below2))
        delete(point_handles_below2);
        point_handles_below2=[];
    end
    if(~isempty(point_handles_below3))
        delete(point_handles_below3);
        point_handles_below3=[];
    end
    if(isempty(follow_bead))
        bn_good=find( (bp(:,1)==current_t) & (bp(:,5)>(current_z-0.5)) & (bp(:,5)<=(current_z+0.5)) & bp(:,6)==1);
        hold on;
        point_handles_good=plot(bp(bn_good,3), bp(bn_good,4), 'go');
        hold off;
        bn_bad=find( (bp(:,1)==current_t) & (bp(:,5)>(current_z-0.5)) & (bp(:,5)<=(current_z+0.5)) & bp(:,6)==0);
        hold on;
        point_handles_bad=plot(bp(bn_bad,3), bp(bn_bad,4), 'kx');
        hold off;
        bn_above1=find( (bp(:,1)==current_t) & (bp(:,5)>((current_z-1)-0.5)) & (bp(:,5)<=((current_z-1)+0.5)));
        bn_above2=find( (bp(:,1)==current_t) & (bp(:,5)>((current_z-2)-0.5)) & (bp(:,5)<=((current_z-2)+0.5)));
        bn_above3=find( (bp(:,1)==current_t) & (bp(:,5)>((current_z-3)-0.5)) & (bp(:,5)<=((current_z-3)+0.5)));
        hold on;
        point_handles_above3=plot(bp(bn_above3,3), bp(bn_above3,4), 'wo');
        point_handles_above2=plot(bp(bn_above2,3), bp(bn_above2,4), 'mo');       
        point_handles_above1=plot(bp(bn_above1,3), bp(bn_above1,4), 'ro');
        hold off;
        bn_below1=find( (bp(:,1)==current_t) & (bp(:,5)>((current_z+1)-0.5)) & (bp(:,5)<=((current_z+1)+0.5)));
        bn_below2=find( (bp(:,1)==current_t) & (bp(:,5)>((current_z+2)-0.5)) & (bp(:,5)<=((current_z+2)+0.5)));
        bn_below3=find( (bp(:,1)==current_t) & (bp(:,5)>((current_z+3)-0.5)) & (bp(:,5)<=((current_z+3)+0.5)));
        hold on;
        point_handles_below1=plot(bp(bn_below1,3), bp(bn_below1,4), 'bo');
        point_handles_below2=plot(bp(bn_below2,3), bp(bn_below2,4), 'co');
        point_handles_below3=plot(bp(bn_below3,3), bp(bn_below3,4), 'yo');
        hold off;
    else
        bn=find( bp(:,1)==current_t & bp(:,2)==follow_bead );
        if( (bp(bn,5)>(current_z-0.5)) & (bp(bn,5)<=(current_z+0.5)) )
            col='go';
        else
            if( bp(bn,5)<current_z )
                col='ro';
            else
                col='bo';
            end
        end
        hold on;
        follow_bead_handle=plot(bp(bn,3), bp(bn,4), col);
        hold off;
    end
end
%
    function plotHeightProfile()
        hprofile=[];
        if(~isempty(high))
            hprofile=[hprofile;high];
        end
        if(~isempty(mid))
            hprofile=[hprofile;mid];
        end
        if(~isempty(low))
            hprofile=[hprofile;low];
        end
        if(~isempty(hprofile_handles))
            delete(hprofile_handles);
            hprofile_handles=[];
        end
        if(~isempty(hprofile_handles_below))
            delete(hprofile_handles_below);
            hprofile_handles_below=[];
        end
        if(~isempty(hprofile_handles))
            delete(hprofile_handles_above);
            hprofile_handles_above=[];
        end
        if~(isempty(hprofile))
            hold on;
            idx=(hprofile(:,3)<current_z);          
            hprofile_handles_above=plot(hprofile(idx,1), hprofile(idx,2), 'rx');
            idx=(hprofile(:,3)==current_z);          
            hprofile_handles=plot(hprofile(idx,1), hprofile(idx,2), 'gx');
            idx=(hprofile(:,3)>current_z);          
            hprofile_handles_below=plot(hprofile(idx,1), hprofile(idx,2), 'bx');
            hold off;
        end
    end
%
function changeGamma(val)
    if(strcmp(mode,'FL'))
        FLgamma=FLgamma+val;
        if(FLgamma<0)
            FLgamma=0;
        end
    else
        BFgamma=BFgamma+val;
        if(BFgamma<0)
            BFgamma=0;
        end
    end
end
%
function resetGamma()
    if(strcmp(mode,'FL'))
        FLgamma=1.0;            
    else
        BFgamma=1.0;           
    end
end
%
    function outlineCell()
        [cell_bw cell_xi cell_yi]=roipoly;
        if(~isempty(cell_bw))
            cell_area=polyarea(cell_xi, cell_yi);
            cell_perimeter=sum(sqrt((cell_xi-circshift(cell_xi,1)).^2+(cell_yi-circshift(cell_yi,1)).^2));
            cell_circularity=cell_perimeter^2/(4*pi*cell_area);        
            %
            [xx yy]=meshgrid(1:imgw, 1:imgh);
            cell_x=sum(sum(xx.*cell_bw))/sum(cell_bw(:));
            cell_y=sum(sum(yy.*cell_bw))/sum(cell_bw(:));
            cell_z=current_z;
            %
            % orientation
            %
            xx=xx-cell_x;
            yy=yy-cell_y;
            cxx=sum(sum(cell_bw.*(xx.^2)));
            cyy=sum(sum(cell_bw.*(yy.^2)));
            cxy=sum(sum(cell_bw.*(xx.*yy)));

            cell_mom=[ cxx cxy; ...
                       cxy cyy];
            [cell_eigvecs cell_eigvals]=eig(cell_mom);

            cell_isoidx=max(cell_eigvals([1 4]))/min(cell_eigvals([1 4]));
            cax=cell_eigvecs(2,:);
            cell_theta=atan2(cax(2), cax(1));
            %
            figure();
            imagesc(imadjust(frame,[],[],gammaval));
            colormap('gray');
            axis off;
            set(gca, 'Position', [0 0 1 1]);
            set(gca, 'YDir', 'normal');
            hold on;
            plot(cell_xi, cell_yi, 'g');
            plot(cell_x, cell_y, 'rp');
            plot([cell_x cell_x+200*cos(cell_theta)], [cell_y cell_y+200*sin(cell_theta)], 'r');
            %
            oname=sprintf('cell_shape_%06d.mat',current_t);
            ok=false;
            if(exist(oname,'file')==2)
                cont=questdlg(sprintf('File %s already exists. Overwrite?', oname), 'Overwrite file?', 'Overwrite', 'Cancel', 'Cancel');
                if(isequal(cont,'Overwrite'))
                    ok=true;
                else
                    ok=false;     
                end
            else
                ok=true;
            end
            if(ok)
                save(oname, 'cell_x', 'cell_y', 'cell_z', 'cell_xi', 'cell_yi', 'cell_bw', 'cell_area', 'cell_perimeter', 'cell_circularity', 'cell_theta', 'cell_eig*', 'cell_isoidx');
            end
        end
    end   
%
end
