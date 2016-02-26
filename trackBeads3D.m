function trackBeads3D(cfg, trackDir)
% 3D bead tracking algorithm.
%
% -------------------------------------------------------------------------
% This file is part of the method published in
%
% Koch TM, Münster S, Bonakdar N, Butler JP, Fabry B (2012) 3D traction 
% forces in cancer cell invasion. PLoS ONE
%
% If you use any part of it, please cite this paper.
% -------------------------------------------------------------------------

% ------------------------------------------------------------------------
% Don't edit this file
% ------------------------------------------------------------------------

% goto dir
saveDir=pwd;
if(~isempty(trackDir))
    cd(trackDir);
end

% ------------------------------------------------------------------------
% Init
% ------------------------------------------------------------------------

% determine number of z and time steps
t=cfg.tstart;
z=0;
while(exist(cfg.FLfileName(z,t),'file'))
    z=z+1;
end
stacksize=z;
buffersize=stacksize;

if(cfg.tend==0)
    t=cfg.tstart;
    while(exist(cfg.FLfileName(0,t),'file'))
        t=t+1;
    end
    tend=t-1;
    tsteps=cfg.tstart:tend;
else
    tsteps=[cfg.tstart cfg.tend];
end

% frame size
ii=imfinfo(cfg.FLfileName(0,cfg.tstart));
imagesize_x=ii.Width;
imagesize_y=ii.Height;

% frame buffers
stack1=[];
stack2=[];
stack_next=1;
stack1_t=-1;
stack2_t=-1;
stack_keep=-1;
% ------------------------------------------------------------------------
% Identify beads
% ------------------------------------------------------------------------
start_coord=[];

%
% bordercutter, i.e. boundary lines of search volume at z
% position of bead
%
bordercutter=zeros(cfg.xsize, cfg.ysize, cfg.zsize);
bordercutter(1,1:end,1+(cfg.zsize-1)/2)=1;
bordercutter(end,1:end,1+(cfg.zsize-1)/2)=1;
bordercutter(1:end,1,1+(cfg.zsize-1)/2)=1;
bordercutter(1:end,end,1+(cfg.zsize-1)/2)=1;

tstep=cfg.tstart;

if(exist(cfg.startposfilename,'file'))    
    bp=dlmread(cfg.startposfilename,'\t');
    start_coord=bp(:,2:5);
    clear bp;
else
    findbar=waitbar(0,'Finding beads...initializing candidates');
    hsteps=length(1:cfg.xsize:(imagesize_x-cfg.xsize))*length(1:cfg.ysize:(imagesize_y-cfg.ysize))*length(1:cfg.zsize:(stacksize-cfg.zsize));
    start_coord=zeros(1E6, 3);
    scc=0;
    %
    threshold=cfg.beadCandidateThreshold; % for decon version
    lcounter=1;
    for z=1:cfg.zsize:(stacksize-cfg.zsize)
        for x=1:cfg.xsize:(imagesize_x-cfg.xsize)
            for y=1:cfg.ysize:(imagesize_y-cfg.ysize)                
                small=getSmallDeconFLStack( [y (y+cfg.ysize-1)], [x (x+cfg.xsize-1)], [z (z+cfg.zsize-1)], tstep, true);
                maxval=max(small(:));
                if(maxval>threshold)
                    idx=find(small==maxval);
                    [sy sx sz] = ind2sub(size(small), idx);
                    sx=sx+x-1;
                    sy=sy+y-1;
                    sz=sz+z-1;                
                    scc=scc+1;
                    start_coord(scc, 1:3)=[sx(1) sy(1) sz(1)];                   
                end
                waitbar(lcounter/hsteps*0.4, findbar);
                lcounter=lcounter+1;
            end
        end
    end    
    resetSmallStack; %
    %
    % sorting z
    %
    [tmp, idx]=sort(start_coord(1:scc,3));
    start_coord2=start_coord;
    start_coord=[];
    start_coord=start_coord2(idx,:);    
    %
    % correcting start coordinates by setting them on the maximum intensity 
    %
    start_coord2=[];
    waitbar(0.4, findbar, 'Finding beads...identifying candidates');
    for m = 1:scc;
        x=(start_coord(m,1));
        y=(start_coord(m,2));
        z=(start_coord(m,3));
        xold = -100;
        yold = -100;
        zold = -100;
        trash = 0;
        lcounter = 0;
        while((xold ~= x || yold ~= y || zold ~= z) && ~trash && lcounter < 100)
            lcounter=lcounter+1;
            xold = x;
            yold = y;
            zold = z;

            if(x <= (cfg.xsize-1)/2 || x >= imagesize_x-(cfg.xsize-1)/2)
                trash = 1;
                continue;
            end
            if(y <= (cfg.ysize-1)/2 || y >= imagesize_y-(cfg.ysize-1)/2)
                trash = 1;
                continue;
            end
            if(z <= (cfg.zsize-1)/2 || z >= stacksize-(cfg.zsize-1)/2)
                trash = 1;
                continue;
            end

            small = getSmallFLStack( [y-(cfg.ysize-1)/2 y+(cfg.ysize-1)/2], [x-(cfg.xsize-1)/2 x+(cfg.xsize-1)/2], [z-(cfg.zsize-1)/2 z+(cfg.zsize-1)/2], tstep, true);
            
            maxidx = find(small==max(small(:)));
            if(length(maxidx)>1)                
                maxidx=maxidx(1);
            end
            [dy, dx, dz] = ind2sub(size(small), maxidx);            
            x=x-(cfg.xsize-1)/2+dx-1;
            y=y-(cfg.ysize-1)/2+dy-1;
            z=z-(cfg.zsize-1)/2+dz-1;            
        end       
        if(~trash)            
            start_coord2=[start_coord2; x y z];
        else
            %
        end
        waitbar(m/scc*0.4+0.4, findbar);
    end  
    [b, m, n] = unique(start_coord2, 'rows');
    start_coord3=start_coord2(m,:);

    % 
    % sort in descending z order
    % 
    tmp=sortrows(start_coord3(:,[3 1 2]));
    start_coord3=zeros(size(start_coord3,1),6);
    start_coord3=tmp(end:-1:1,[2 3 1]);

    %
    % good beads?
    %
    waitbar(0.8, findbar, 'Finding beads...eliminating bad candidates');
    for n=1:size(start_coord3,1)
        xr = round(start_coord3(n,1));
        yr = round(start_coord3(n,2));
        zr = round(start_coord3(n,3));

        small = getSmallFLStack( [(yr-(cfg.ysize-1)/2) (yr+(cfg.ysize-1)/2)], [(xr-(cfg.xsize-1)/2) (xr+(cfg.xsize-1)/2)], [(zr-(cfg.zsize-1)/2) (zr+(cfg.zsize-1)/2)], tstep, false);

        %
        % Good Beads: Intesity has to drop in all directions from central pixel
        %
        % Criterium: 
        % 1. Integrate small along z
        % 2. Max ratio of borderintensity / central intensity not allowed to
        % exceed threshold (0.7)
        %
        projected_small=sum(small,3);
        I_center=projected_small(1+(cfg.ysize-1)/2, 1+(cfg.xsize-1)/2);
        projected_border=[projected_small(1,:) projected_small(:,end)' projected_small(end,:) projected_small(:,1)'];    
        Imax_borderratio=max(projected_border./I_center);
        %
        start_coord3(n,5)=I_center;
        start_coord3(n,6)=Imax_borderratio;
        if(Imax_borderratio>cfg.goodBeadThreshold) 
            start_coord3(n,4)=0; % bad bead
        else
            start_coord3(n,4)=1; % good bead
        end
        waitbar(n/size(start_coord3,1)*0.2+0.8, findbar);
    end
    close(findbar);
    %
    % Eliminating bad start positions
    %
    idx=find(start_coord3(:,4)==0);
    %bad_start_coord=start_coord3(idx, :);
    start_coord3(idx,:)=[];

    %
    % assign bead numbers and write out
    %
    start_coord=[(1:size(start_coord3,1))' start_coord3(:,1:3)];          
    ofile = fopen(cfg.startposfilename, 'w'); % Warning! File will be overwritten without question at this point.
    for n=1:size(start_coord,1)
        fprintf(ofile, '%f\t%f\t%f\t%f\t%f\t%f\t%f\n', 0, start_coord(n,1), start_coord(n,2), start_coord(n,3), start_coord(n,4), start_coord3(n,5), start_coord3(n,6));
    end
    fclose(ofile);
end

% ------------------------------------------------------------------------
% Track beads
% ------------------------------------------------------------------------

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interolation setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r_sign=[1 1 -1 -1 1 1 -1 -1];
    s_sign=[1 1 1 1 -1 -1 -1 -1];
    t_sign=[1 -1 1 -1 1 -1 1 -1];
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Global XCORR Setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gxc_sizex=floor(imagesize_y/4);
    gxc_sizey=floor(imagesize_y/4);
    gxc_startx=floor(imagesize_x/2-gxc_sizex/2);
    gxc_starty=floor(imagesize_y/2-gxc_sizey/2);
    gxc_sizez=floor(buffersize/4); %floor(0.4*stacksize); %40;
    gxc_startz=floor(stacksize/2-gxc_sizez/2); %30;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initializing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bead_coord = [];
    start_coord(:,5)=1;
    if(~exist('sequentialTracking', 'var'))
        sequentialTracking=false;
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % tracking all times
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    for tstep=tsteps
        trackbar=waitbar(0,sprintf('Tracking %d beads at timestep %d.\n', size(start_coord,1),tstep));
        if(tstep>cfg.tstart)
            if(sequentialTracking)
                tprev=tstep-1;
            else
                tprev=cfg.tstart;
                stack_keep=tprev;
            end
            prevstack=getSmallFLStack([gxc_starty (gxc_starty+gxc_sizey)], [gxc_startx (gxc_startx+gxc_sizex)], [gxc_startz (gxc_startz+gxc_sizez)], tprev, true);
        end      
        bead_coord = zeros(length(start_coord),10);       

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % global cross correlation, i.e. drift correction
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        doGlobalCorrection=true;
        x_correct=0;
        y_correct=0;
        z_correct=0;
        if( tstep>cfg.tstart && doGlobalCorrection)           
            %
            % cross correlation
            %                       
            xc_d=[0 0 0; -1 0 0; 1 0 0; 0 -1 0; 0 1 0; 0 0 -1; 0 0 1];
            xcval=ones(7,1).*1E99;
            busy=true;
            while(busy)
                for n=1:7
                    xcsmall=getSmallFLStack( [(gxc_starty+xc_d(n,2)+y_correct) (gxc_starty+y_correct+xc_d(n,2)+gxc_sizey)], [(gxc_startx+x_correct+xc_d(n,1)) (gxc_startx+x_correct+xc_d(n,1)+gxc_sizex)], [(gxc_startz+z_correct+xc_d(n,3)) (gxc_startz+z_correct+xc_d(n,3)+gxc_sizez)], tstep, true);
                    xcval(n)=sum(sum(sum( (prevstack-xcsmall).^2 )));                 
                end            
                idx=find(xcval==min(xcval));
                if(idx==1) 
                    busy=false;
                end
                x_correct=x_correct+xc_d(idx,1);
                y_correct=y_correct+xc_d(idx,2);
                z_correct=z_correct+xc_d(idx,3);
            end            
        else
            if(tstep>cfg.tstart)
                fprintf('GXC switched off.\n');
            end
        end    
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        

        x_lastcorrect=0;
        y_lastcorrect=0;
        z_lastcorrect=0;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % track bead
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for n = 1:size(start_coord,1);
            bead_num = start_coord(n,1);
            x = start_coord(n,2);
            y = start_coord(n,3);
            z = start_coord(n,4);

            if(start_coord(n,5)) % good bead?
                trash = 0;
            else 
                trash = 1;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % individual cross correlation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            doIndividualCorrection=true;
            if(tstep>cfg.tstart && ~trash && doIndividualCorrection)
                %
                % xcorr
                %
                xpx=round(x);
                ypx=round(y);
                zpx=round(z);               
                %
                xc_xsize=25;
                xc_ysize=25;
                xc_zsize=7;
                %
                % cross correlation
                %
                busy=true;                
                xbmin=(xpx)-(xc_xsize-1)/2;
                xbmax=(xpx)+(xc_xsize-1)/2;
                ybmin=(ypx)-(xc_ysize-1)/2;
                ybmax=(ypx)+(xc_ysize-1)/2;
                zbmin=(zpx)-(xc_zsize-1)/2;
                zbmax=(zpx)+(xc_zsize-1)/2;
                if(xbmin < 1 || xbmax > imagesize_x)
                    trash = 1;
                    busy=false;                
                end
                if(ybmin < 1 || ybmax > imagesize_y)
                    trash = 1;
                    busy=false;                
                end
                if(zbmin < 1 || zbmax > stacksize)
                    trash = 1;
                    busy=false;                
                end
                if(~trash)
                    prevsmall=getSmallFLStack( [ybmin ybmax], [xbmin xbmax], [zbmin zbmax], tprev, true);            
                end                
                x_tmpcorr=x_correct; % initialize individual correction vector with global dc
                y_tmpcorr=y_correct;
                z_tmpcorr=z_correct;
                                                                       
                xc_d=[0 0 0; -1 0 0; 1 0 0; 0 -1 0; 0 1 0; 0 0 -1; 0 0 1];
                xcval=ones(7,1).*1E99;
                while(busy)
                    for xcn=1:7
                        xbmin=(xpx+x_tmpcorr+xc_d(xcn,1))-(xc_xsize-1)/2;
                        xbmax=(xpx+x_tmpcorr+xc_d(xcn,1))+(xc_xsize-1)/2;
                        ybmin=(ypx+y_tmpcorr+xc_d(xcn,2))-(xc_ysize-1)/2;
                        ybmax=(ypx+y_tmpcorr+xc_d(xcn,2))+(xc_ysize-1)/2;
                        zbmin=(zpx+z_tmpcorr+xc_d(xcn,3))-(xc_zsize-1)/2;
                        zbmax=(zpx+z_tmpcorr+xc_d(xcn,3))+(xc_zsize-1)/2;
                        if(xbmin < 1 || xbmax > imagesize_x)
                            trash = 1;
                            busy=false;
                            continue;
                        end
                        if(ybmin < 1 || ybmax > imagesize_y)
                            trash = 1;
                            busy=false;
                            continue;
                        end
                        if(zbmin < 1 || zbmax > stacksize)
                            trash = 1;
                            busy=false;
                            continue;
                        end
                        xcsmall=getSmallFLStack( [ybmin ybmax], [xbmin xbmax], [zbmin zbmax], tstep, true); 
                        xcval(xcn)=sum(sum(sum( (prevsmall-xcsmall).^2 )));                          
                    end
                    if(busy)
                        idx=find(xcval==min(xcval));
                        if(length(idx)>1)
                            idx=idx(1);
                            trash=1;
                            busy=false;
                        end
                        if(idx==1) 
                            busy=false;
                            x_lastcorrect=x_tmpcorr;
                            y_lastcorrect=y_tmpcorr;
                            z_lastcorrect=z_tmpcorr;    
                        else
                            x_tmpcorr=x_tmpcorr+xc_d(idx,1);
                            y_tmpcorr=y_tmpcorr+xc_d(idx,2);
                            z_tmpcorr=z_tmpcorr+xc_d(idx,3);
                        end
                    end
                end
                %
                % fail safe: fall back to finding global minimum if correction
                % is too excessive
                %
                if( abs(x_lastcorrect-x_correct)>5 || abs(y_lastcorrect-y_correct)>5 || abs(z_lastcorrect-z_correct)>3) 
                    xbmin=(xpx+x_correct-20)-(xc_xsize-1)/2;
                    xbmax=(xpx+x_correct+20)+(xc_xsize-1)/2;
                    ybmin=(ypx+y_correct-20)-(xc_ysize-1)/2;
                    ybmax=(ypx+y_correct+20)+(xc_ysize-1)/2;
                    zbmin=(zpx+z_correct-6)-(xc_zsize-1)/2;
                    zbmax=(zpx+z_correct+6)+(xc_zsize-1)/2;
                    if(xbmin < 1 || xbmax > imagesize_x)
                        trash = 1;                                        
                    end
                    if(ybmin < 1 || ybmax > imagesize_y)
                        trash = 1;                    
                    end
                    if(zbmin < 1 || zbmax > stacksize)
                        trash = 1;                    
                    end
                    if(~trash)
                        ssd=zeros(41,41,13);                
                        for fsdc_i=-20:20
                            for fsdc_j=-20:20
                                for fsdc_k=-6:6
                                    xbmin=(xpx+x_correct+fsdc_i)-(xc_xsize-1)/2;
                                    xbmax=(xpx+x_correct+fsdc_i)+(xc_xsize-1)/2;
                                    ybmin=(ypx+y_correct+fsdc_j)-(xc_ysize-1)/2;
                                    ybmax=(ypx+y_correct+fsdc_j)+(xc_ysize-1)/2;
                                    zbmin=(zpx+z_correct+fsdc_k)-(xc_zsize-1)/2;
                                    zbmax=(zpx+z_correct+fsdc_k)+(xc_zsize-1)/2;
                                    xcsmall=getSmallFLStack( [ybmin ybmax], [xbmin xbmax], [zbmin zbmax], tstep, true);
                                    ssd(fsdc_i+21,fsdc_j+21,fsdc_k+7)=sum(sum(sum( (prevsmall-xcsmall).^2 )));                              
                                end
                            end
                        end               
                        %
                        % minimum and update
                        %
                        [fsdc_i fsdc_j fsdc_k]=ind2sub(size(ssd), find(ssd==min(ssd(:))));
                        x=x+x_correct+fsdc_i-21;
                        y=y+y_correct+fsdc_j-21;
                        z=z+z_correct+fsdc_k-7;             
                        %
                    end
                else
                    %
                    % px update
                    %
                    x=x+x_lastcorrect;
                    y=y+y_lastcorrect;
                    z=z+z_lastcorrect; 
                end
            else               
                 if(tstep>cfg.tstart && ~trash)
                    if(doGlobalCorrection)
                        x=x+x_correct;
                        y=y+y_correct;
                        z=z+z_correct;
                    end
                end                
            end     
            %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Interpolation in xyz-direction
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            xsave=x;
            ysave=y;
            zsave=z;
            x_prev = round(start_coord(n,2));
            y_prev = round(start_coord(n,3));
            z_prev = round(start_coord(n,4));
            x_curr=round(x);
            y_curr=round(y);
            z=round(z);
            %
            %
            if(tstep>cfg.tstart && ~trash) % Good bead
                %
                % track
                %
                if(x_prev <= ((cfg.xsize-1)/2+1) || x_prev >= imagesize_x-(cfg.xsize-1)/2-1)
                    trash = 1;               
                elseif(y_prev <= ((cfg.ysize-1)/2+1) || y_prev >= imagesize_y-(cfg.ysize-1)/2-1)
                    trash = 1;
                elseif(z_prev <= ((cfg.zsize-1)/2+1) || z_prev >= stacksize-(cfg.zsize-1)/2-1)
                    trash = 1;
                elseif(x_curr <= (cfg.xsize-1)/2 || x_curr >= imagesize_x-(cfg.xsize-1)/2)
                    trash = 1;
                elseif(y_curr <= (cfg.ysize-1)/2 || y_curr >= imagesize_y-(cfg.ysize-1)/2)
                    trash = 1;
                elseif(z <= ((cfg.zsize-1)/2) || z >= stacksize-(cfg.zsize-1)/2)
                    trash = 1;                
                else
                    %
                    % bead pos now is in px accuracy
                    % subpixel by interpolation
                    %        
                    %
                    % interpolation
                    %
                    prevsmall=getSmallFLStack( round([(y_prev-(cfg.ysize-1)/2)-1 (y_prev+(cfg.ysize-1)/2)+1]), round([(x_prev-(cfg.xsize-1)/2-1) (x_prev+(cfg.xsize-1)/2)+1]), ([(z_prev-(cfg.zsize-1)/2)-1 (z_prev+(cfg.zsize-1)/2)+1]), tprev, true); 
                    small = getSmallFLStack( round([(y_curr-(cfg.ysize-1)/2) (y_curr+(cfg.ysize-1)/2)]), round([(x_curr-(cfg.xsize-1)/2) (x_curr+(cfg.xsize-1)/2)]), ([(z-(cfg.zsize-1)/2) (z+(cfg.zsize-1)/2)]), tstep, true);     
                    [r(1),s(1),t(1),f(1)]=minimizeInterpolatedSmallSD(prevsmall(1:cfg.xsize+1,1:cfg.ysize+1,1:cfg.zsize+1),small,0,0,0);  
                    [r(2),s(2),t(2),f(2)]=minimizeInterpolatedSmallSD(prevsmall(1:cfg.xsize+1,1:cfg.ysize+1,2:cfg.zsize+2),small,0,0,0);
                    [r(3),s(3),t(3),f(3)]=minimizeInterpolatedSmallSD(prevsmall(1:cfg.xsize+1,2:cfg.ysize+2,1:cfg.zsize+1),small,0,0,0);
                    [r(4),s(4),t(4),f(4)]=minimizeInterpolatedSmallSD(prevsmall(1:cfg.xsize+1,2:cfg.ysize+2,2:cfg.zsize+2),small,0,0,0);
                    [r(5),s(5),t(5),f(5)]=minimizeInterpolatedSmallSD(prevsmall(2:cfg.xsize+2,1:cfg.ysize+1,1:cfg.zsize+1),small,0,0,0);
                    [r(6),s(6),t(6),f(6)]=minimizeInterpolatedSmallSD(prevsmall(2:cfg.xsize+2,1:cfg.ysize+1,2:cfg.zsize+2),small,0,0,0);
                    [r(7),s(7),t(7),f(7)]=minimizeInterpolatedSmallSD(prevsmall(2:cfg.xsize+2,2:cfg.ysize+2,1:cfg.zsize+1),small,0,0,0);
                    [r(8),s(8),t(8),f(8)]=minimizeInterpolatedSmallSD(prevsmall(2:cfg.xsize+2,2:cfg.ysize+2,2:cfg.zsize+2),small,0,0,0);                    
                    %
                    % min
                    %
                    idx=(abs(r)<=1.001) & (abs(s)<=1.001) & (abs(t)<=1.001);
                    sub_dx=-0.5*r(idx)+r_sign(idx)*0.5;
                    sub_dy=-0.5*s(idx)+s_sign(idx)*0.5;
                    sub_dz=-0.5*t(idx)+t_sign(idx)*0.5;
                    if(isempty(sub_dx))                        
                        idx=(f==min(f));
                        sub_dx=-0.5*r(idx)+r_sign(idx)*0.5;
                        sub_dy=-0.5*s(idx)+s_sign(idx)*0.5;
                        sub_dz=-0.5*t(idx)+t_sign(idx)*0.5;                        
                        if(isempty(sub_dx)) % happens if f is all NaN, i.e. Newton-Raphson diverged
                            sub_dx=1E6;
                            sub_dy=1E6;
                            sub_dz=1E6;
                            trash=1;                            
                        elseif( abs(sub_dx(1))>=1 || abs(sub_dy(1))>=1 || abs(sub_dz(1))>=1) % More than subpix shift                             
                            trash=1;                            
                        end
                    else
                        if(length(sub_dx)>1)
                            idx2=(f(idx)==min(f(idx)));
                            sub_dx=sub_dx(idx2);
                            sub_dy=sub_dy(idx2);
                            sub_dz=sub_dz(idx2);                       
                        end                         
                    end
                    x=xsave+sub_dx(1);
                    y=ysave+sub_dy(1);
                    z=zsave+sub_dz(1);                          
                end
            end
            %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % saving bead pos in memory
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
            if(~trash)
                xr = round(x);
                yr = round(y);
                zr = round(z);
                small = getSmallFLStack( [(yr-(cfg.ysize-1)/2) (yr+(cfg.ysize-1)/2)], [(xr-(cfg.xsize-1)/2) (xr+(cfg.xsize-1)/2)], [(zr-(cfg.zsize-1)/2) (zr+(cfg.zsize-1)/2)], tstep, true); 

                projected_small=sum(small,3);
                I_center=projected_small(1+(cfg.ysize-1)/2, 1+(cfg.xsize-1)/2);
                projected_border=[projected_small(1,:) projected_small(:,end)' projected_small(end,:) projected_small(:,1)'];    
                Imax_borderratio=max(projected_border./I_center);

                bead_coord(n,:) = [bead_num x y z 1 I_center Imax_borderratio x_lastcorrect y_lastcorrect z_lastcorrect];
                %
            else % bad bead
                bead_coord(n,:) = [bead_num x y z 0 0 0 x_lastcorrect y_lastcorrect z_lastcorrect];
            end
            waitbar(n/size(start_coord,1),trackbar);
            %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        end 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % writing bead coord to file
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(tstep==cfg.tstart && exist(cfg.beadposfilename))
            delete(cfg.beadposfilename);
        end
        ofile = fopen(cfg.beadposfilename, 'a');
        for n=1:size(bead_coord,1)
            % File format
            % 1:timestep  2:bead number 3:x 4:y 5:z 6:valid bead 7:I_center
            % 8:Imax_borderratio 9:x_lastcorrect 10:y_lastcorrect 11:z_lastcorrect
            fprintf(ofile, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', tstep, bead_coord(n,1), bead_coord(n,2), bead_coord(n,3), bead_coord(n,4), bead_coord(n,5), bead_coord(n,6), bead_coord(n,7), bead_coord(n,8), bead_coord(n,9), bead_coord(n,10));
        end
        fclose(ofile);       
        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % preparing variables for next step
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        if(sequentialTracking)
            start_coord = (bead_coord(:, 1:5));
        else
            if(tstep==cfg.tstart)
                start_coord = (bead_coord(:, 1:5));
            end
        end            
        close(trackbar);
    end
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


% ------------------------------------------------------------------------
% Subfunctions
% ------------------------------------------------------------------------

function small=getSmallFLStack(xrange, yrange, zrange, t, normalize)
%
% access to small subvolume in stack at time t.
% loads images, if necessary
%

if( (t~=stack1_t) && (t~=stack2_t) ) 
    %
    % init buffer
    %
    if(stack_next==1)
        h=waitbar(0,sprintf('Filling image buffer 1 with timestep %d', t));        
        for i=1:(buffersize);             
            fname = cfg.FLfileName(i-1,t);
            im=uint16(imread(fname)-cfg.CCDoffset);            
            stack1(i).im = im;
            stack1(i).mean=mean(im(:));
            waitbar(i/buffersize,h);
        end
        stack1_t=t;
        if(stack2_t==stack_keep)
            stack_next=1;
        else
            stack_next=2;
        end
        close(h);
    else
        h=waitbar(0,sprintf('Filling image buffer 2 with timestep %d', t));        
        for i=1:(buffersize); 
            fname = cfg.FLfileName(i-1,t);
            im=uint16(imread(fname)-cfg.CCDoffset);
            stack2(i).im = im;
            stack2(i).mean=mean(im(:));
            waitbar(i/buffersize,h);
        end        
        stack2_t=t;
        if(stack1_t==stack_keep)
            stack_next=2;
        else
            stack_next=1;
        end
        close(h);
    end
end
   
%
% extract small
%
smx=length(xrange(1):xrange(2));
smy=length(yrange(1):yrange(2));
smz=length(zrange(1):zrange(2));
small=zeros(smx, smy, smz, 'double');
k=1;
if(t==stack1_t)    
    for i=zrange(1):zrange(2)
        if(normalize)
            small(:,:,k)=double(stack1(i).im(xrange(1):xrange(2), yrange(1):yrange(2)))./stack1(i).mean;
        else
            small(:,:,k)=double(stack1(i).im(xrange(1):xrange(2), yrange(1):yrange(2)));
        end
        k=k+1;
    end
elseif(t==stack2_t)    
    for i=zrange(1):zrange(2)
        if(normalize)
            small(:,:,k)=double(stack2(i).im(xrange(1):xrange(2), yrange(1):yrange(2)))./stack2(i).mean;
        else
            small(:,:,k)=double(stack2(i).im(xrange(1):xrange(2), yrange(1):yrange(2)));
        end
        k=k+1;
    end
else
    error('Fatal error!');
end

end

function small=getSmallDeconFLStack(xrange, yrange, zrange, t, normalize)
%
% access to small subvolume in stack at time t.
% loads images, if necessary
% processes images with a simple deconvolution scheme first.
%

if( (t~=stack1_t) && (t~=stack2_t) ) 
    %
    % init buffer
    %
    if(stack_next==1)              
        h=waitbar(0,sprintf('Filling image buffer 1 with timestep %d applying deconvolution', t));
        for i=1:(buffersize); 
            fname = cfg.FLfileName(i-1,t);
            im=uint16(imread(fname)-cfg.CCDoffset); 
            im=im-imfilter(im,fspecial('gauss',[5 5],1));
            stack1(i).im = im;
            stack1(i).mean=mean(im(:));
            waitbar(i/buffersize,h);
        end
        stack1_t=t;
        stack_next=2;
        close(h);
    else
        h=waitbar(0,sprintf('Filling image buffer 2 with timestep %d applying deconvolution', t));
        for i=1:(buffersize); 
            fname = cfg.FLfileName(i-1,t);
            im=uint16(imread(fname)-cfg.CCDoffset); 
            im=im-imfilter(im,fspecial('gauss',[5 5],1));
            stack2(i).im = im;
            stack2(i).mean=mean(im(:));
            waitbar(i/buffersize,h);
        end
        stack2_t=t;
        stack_next=1;
        close(h);
    end
end
   
%
% extract small
%
smx=length(xrange(1):xrange(2));
smy=length(yrange(1):yrange(2));
smz=length(zrange(1):zrange(2));
small=zeros(smx, smy, smz, 'double');
k=1;
if(t==stack1_t)
    %fprintf('Using stack 1\n');
    for i=zrange(1):zrange(2)
        if(normalize)
            small(:,:,k)=double(stack1(i).im(xrange(1):xrange(2), yrange(1):yrange(2)))./stack1(i).mean;
        else
            small(:,:,k)=double(stack1(i).im(xrange(1):xrange(2), yrange(1):yrange(2)));
        end
        k=k+1;
    end
elseif(t==stack2_t)
    %fprintf('Using stack 2\n');
    for i=zrange(1):zrange(2)
        if(normalize)
            small(:,:,k)=double(stack2(i).im(xrange(1):xrange(2), yrange(1):yrange(2)))./stack2(i).mean;
        else
            small(:,:,k)=double(stack2(i).im(xrange(1):xrange(2), yrange(1):yrange(2)));
        end
        k=k+1;
    end
else
    error('Fatal error!');
end
end

function resetSmallStack()
    stack1=[];
    stack2=[];
    stack_next=1;
    stack1_t=-1;
    stack2_t=-1;
    stack_keep=-1;
end

function [r, s, t, f]=minimizeInterpolatedSmallSD(prevsmall, small, start_r, start_s, start_t)
%
% [r, s, t]=minimizeInterpolatedSmallSD(prevsmall, small, start_r, start_s, start_t)
%
% minimizes squared difference between interpolated prevsmall shifted by
% r, s, t and small using Newton-Raphson method.
% 
% prevsmall has size(small)+1 
%
    persistent PI PJ PK;
    
    if(isempty(PI))
        PI(:,1)=2:size(prevsmall,1);
        PI(:,2)=2:size(prevsmall,1);
        PI(:,3)=1:size(prevsmall,1)-1;
        PI(:,4)=1:size(prevsmall,1)-1;
        PI(:,5)=2:size(prevsmall,1);
        PI(:,6)=2:size(prevsmall,1);
        PI(:,7)=1:size(prevsmall,1)-1;
        PI(:,8)=1:size(prevsmall,1)-1;

        PJ(:,1)=2:size(prevsmall,2);
        PJ(:,2)=1:size(prevsmall,2)-1;
        PJ(:,3)=1:size(prevsmall,2)-1;
        PJ(:,4)=2:size(prevsmall,2);
        PJ(:,5)=2:size(prevsmall,2);
        PJ(:,6)=1:size(prevsmall,2)-1;
        PJ(:,7)=1:size(prevsmall,2)-1;
        PJ(:,8)=2:size(prevsmall,2);

        PK(:,1)=2:size(prevsmall,3);
        PK(:,2)=2:size(prevsmall,3);
        PK(:,3)=2:size(prevsmall,3);
        PK(:,4)=2:size(prevsmall,3);
        PK(:,5)=1:size(prevsmall,3)-1;
        PK(:,6)=1:size(prevsmall,3)-1;
        PK(:,7)=1:size(prevsmall,3)-1;
        PK(:,8)=1:size(prevsmall,3)-1;
    end

    a=zeros(1,8);
    A=zeros(8,8);
    for i=1:8
        a(i)=sum(sum(sum( small.*prevsmall(PI(:,i), PJ(:,i), PK(:,i)) )));
        for j=1:8
            A(i,j)=sum(sum(sum( prevsmall(PI(:,i), PJ(:,i), PK(:,i)) .* prevsmall(PI(:,j), PJ(:,j), PK(:,j)) )));
        end
    end

    r_new=start_r;
    s_new=start_s;
    t_new=start_t;
    %
    % Newton Raphson scheme
    %
    iter=0;
    crit=1;
    while(crit>0.001)
        r_old=r_new;
        s_old=s_new;
        t_old=t_new;
        %
        rp=(1+r_old);
        rm=(1-r_old);
        sp=(1+s_old);
        sm=(1-s_old);
        tp=(1+t_old);
        tm=(1-t_old);
        h(1)=1/8*(rp)*(sp)*(tp);
        h(2)=1/8*(rm)*(sp)*(tp);
        h(3)=1/8*(rm)*(sm)*(tp);
        h(4)=1/8*(rp)*(sm)*(tp);
        h(5)=1/8*(rp)*(sp)*(tm);
        h(6)=1/8*(rm)*(sp)*(tm);
        h(7)=1/8*(rm)*(sm)*(tm);
        h(8)=1/8*(rp)*(sm)*(tm);
        %
        % first derivatives
        %
        hr(1)= 1/8*(sp)*(tp);
        hr(2)=-1/8*(sp)*(tp);
        hr(3)=-1/8*(sm)*(tp);
        hr(4)= 1/8*(sm)*(tp);
        hr(5)= 1/8*(sp)*(tm);
        hr(6)=-1/8*(sp)*(tm);
        hr(7)=-1/8*(sm)*(tm);
        hr(8)= 1/8*(sm)*(tm);
        hs(1)= 1/8*(rp)*(tp);
        hs(2)= 1/8*(rm)*(tp);
        hs(3)=-1/8*(rm)*(tp);
        hs(4)=-1/8*(rp)*(tp);
        hs(5)= 1/8*(rp)*(tm);
        hs(6)= 1/8*(rm)*(tm);
        hs(7)=-1/8*(rm)*(tm);
        hs(8)=-1/8*(rp)*(tm);
        ht(1)= 1/8*(rp)*(sp);
        ht(2)= 1/8*(rm)*(sp);
        ht(3)= 1/8*(rm)*(sm);
        ht(4)= 1/8*(rp)*(sm);
        ht(5)=-1/8*(rp)*(sp);
        ht(6)=-1/8*(rm)*(sp);
        ht(7)=-1/8*(rm)*(sm);
        ht(8)=-1/8*(rp)*(sm);
        %
        % 2nd derivatives
        %
        hrs(1)= 1/8*(tp);
        hrs(2)=-1/8*(tp);
        hrs(3)= 1/8*(tp);
        hrs(4)=-1/8*(tp);
        hrs(5)= 1/8*(tm);
        hrs(6)=-1/8*(tm);
        hrs(7)= 1/8*(tm);
        hrs(8)=-1/8*(tm);
        hrt(1)= 1/8*(sp);
        hrt(2)=-1/8*(sp);
        hrt(3)=-1/8*(sm);
        hrt(4)= 1/8*(sm);
        hrt(5)=-1/8*(sp);
        hrt(6)= 1/8*(sp);
        hrt(7)= 1/8*(sm);
        hrt(8)=-1/8*(sm);
        hst(1)= 1/8*(rp);
        hst(2)= 1/8*(rm);
        hst(3)=-1/8*(rm);
        hst(4)=-1/8*(rp);
        hst(5)=-1/8*(rp);
        hst(6)=-1/8*(rm);
        hst(7)= 1/8*(rm);
        hst(8)= 1/8*(rp);
        %
        % F and the Jacobian matrix
        %
        F(1)=0;
        F(2)=0;
        F(3)=0;
        J=zeros(3,3);
        for i=1:8
            F(1)=F(1)-2*(a(i)*hr(i));
            F(2)=F(2)-2*(a(i)*hs(i));
            F(3)=F(3)-2*(a(i)*ht(i));
            %       
            J(1,2)=J(1,2)-2*a(i)*hrs(i);
            J(1,3)=J(1,3)-2*a(i)*hrt(i);        
            J(2,3)=J(2,3)-2*a(i)*hst(i);        
            for j=1:8
                F(1)=F(1)+2*A(i,j)*h(i)*hr(j);        
                F(2)=F(2)+2*A(i,j)*h(i)*hs(j);      
                F(3)=F(3)+2*A(i,j)*h(i)*ht(j);
                %
                J(1,1)=J(1,1)+2*A(i,j)*(hr(i)*hr(j));
                J(1,2)=J(1,2)+2*A(i,j)*(hs(i)*hr(j)+h(i)*hrs(j));
                J(1,3)=J(1,3)+2*A(i,j)*(ht(i)*hr(j)+h(i)*hrt(j));           
                J(2,2)=J(2,2)+2*A(i,j)*(hs(i)*hs(j));
                J(2,3)=J(2,3)+2*A(i,j)*(ht(i)*hs(j)+h(i)*hst(j));            
                J(3,3)=J(3,3)+2*A(i,j)*(ht(i)*ht(j));
            end
        end
        J(2,1)=J(1,2);
        J(3,1)=J(1,3);
        J(3,2)=J(2,3);
        %
        % Newton-Raphson update
        %
        dx=J\(-F)';
        r_new=r_old+dx(1);
        s_new=s_old+dx(2);
        t_new=t_old+dx(3);
        crit=max(abs([r_new-r_old s_new-s_old t_new-t_old]));
        iter=iter+1;
        if(iter>100)
            crit=0;
            r_new=100;
            s_new=100;
            t_new=100;
        end
    end
    r=r_new;
    s=s_new;
    t=t_new;
    %
    % Calculate last F
    %
    r_old=r_new;
    s_old=s_new;
    t_old=t_new;
    rp=(1+r_old);
    rm=(1-r_old);
    sp=(1+s_old);
    sm=(1-s_old);
    tp=(1+t_old);
    tm=(1-t_old);
    %
    h(1)=1/8*(rp)*(sp)*(tp);
    h(2)=1/8*(rm)*(sp)*(tp);
    h(3)=1/8*(rm)*(sm)*(tp);
    h(4)=1/8*(rp)*(sm)*(tp);
    h(5)=1/8*(rp)*(sp)*(tm);
    h(6)=1/8*(rm)*(sp)*(tm);
    h(7)=1/8*(rm)*(sm)*(tm);
    h(8)=1/8*(rp)*(sm)*(tm);             
    %
    % Calculate f
    %
    f=sum(sum(sum(...
      (small-h(1)*prevsmall(PI(:,1), PJ(:,1), PK(:,1))...
            -h(2)*prevsmall(PI(:,2), PJ(:,2), PK(:,2))...
            -h(3)*prevsmall(PI(:,3), PJ(:,3), PK(:,3))...
            -h(4)*prevsmall(PI(:,4), PJ(:,4), PK(:,4))...
            -h(5)*prevsmall(PI(:,5), PJ(:,5), PK(:,5))...
            -h(6)*prevsmall(PI(:,6), PJ(:,6), PK(:,6))...
            -h(7)*prevsmall(PI(:,7), PJ(:,7), PK(:,7))...
            -h(8)*prevsmall(PI(:,8), PJ(:,8), PK(:,8)) ).^2 )));
end

% ------------------------------------------------------------------------
% End
% ------------------------------------------------------------------------
% return to prev dir
cd(saveDir);
end    