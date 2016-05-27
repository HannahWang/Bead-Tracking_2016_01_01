%
% Strain Energy Results
%
% After the bead tracking and strain energy computation is finished, this
% script computes the strain energy in a subvolume around the cell and
% visualizes results. 
%
% If the location of the cell is known and stored in 'cell_shape_*.mat',
% this location will be the center of the coordinate systems used here.
% Otherwise the center of the subvolume will be used.
%
% -------------------------------------------------------------------------
% This file is part of the method published in
%
% Koch TM, Münster S, Bonakdar N, Butler JP, Fabry B (2012) 3D traction 
% forces in cancer cell invasion. PLoS ONE
%
% If you use any part of it, please cite this paper.
% -------------------------------------------------------------------------

% x/y/z axis length [µm] for SE density plot
sedensplot_length_x_mum=170; 
sedensplot_length_y_mum=170; 
sedensplot_length_z_mum=400;

% SE density iso value for isosurface, relative to maximum.
sedensplot_relativeIsoval=0.36;

% filename for exported SE density figure, leave '' to switch off
% use '.png' file extension for high quality export, '*.fig' to save as
% Matlab figure.
sedensplot_filename='sedens.fig'; 

% filename for plot of strain energy vs. integration volume, leave '' to
% switch off. currently only Matlab '*.fig' files
sevol_filename='sevol.fig';

% box length of integration volume for total strain energy
se_length_x_mum=170;
se_length_y_mum=170;
se_length_z_mum=15;

% mat file name to save cell strain energy, leave '' to switch off
se_filename='cell_strain_energy.mat';

% px dimensions of displacement image around cell
celldisplot_range_x_px=801;
celldisplot_range_y_px=801;    % was 401

% scale factor for projected cell displacments, set to 0 to autoscale
celldisplot_scalef=0;

% filename for plot of projected displacements around cell, leave '' to
% switch off. currently only Matlab '*.fig' files
celldisplot_filename='cell_displacements.fig';

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% Don't edit below this line
% ------------------------------------------------------------------------
cfg=getStrainEnergyAssayConfig();
set(0, 'Units', 'centimeters')
scrsz=get(0,'ScreenSize');
%
seres=load(cfg.resultsfilename);
%
% SE density projection figure
%
Ugs=smooth3(seres.Ucdg_near,'box',[5 5 5]);
%
sedensfig=figure();
set(gcf, 'Units', 'centimeters', 'Position', [scrsz(3)-12 scrsz(4)-12 8 8]);
set(gca, 'Position', [0 0 1 1]);
set(gcf,'PaperPositionMode', 'auto', 'Color', 'w');
sedensIsoProjectionFigure(cfg,Ugs,seres.xg,seres.yg,seres.zg,sedensplot_length_x_mum,sedensplot_length_y_mum,sedensplot_length_z_mum,sedensplot_relativeIsoval,sedensplot_filename);
%
% integrated SE vs integration volume
%
volirange=cfg.gridsize:cfg.gridsize:(0.2*(max(seres.xg(:))-min(seres.xg(:))));
Uvolsym=zeros(length(volirange)+1,1);
Vvolsym=zeros(length(volirange)+1,1);
k=2; % first entry remains 0
for i=volirange
    [nx, ny, nz, Udens_sub]=subvolume(seres.xg,seres.yg,seres.zg,seres.Ucdg_near, [-2*i,2*i,-2*i,2*i,-i,i]);
    Uvolsym(k)=sum(Udens_sub(:).*cfg.gridsize^3);
    Vvolsym(k)=length(Udens_sub(:))*cfg.gridsize^3; 
    k=k+1;
end
sevolfig=figure();
set(gcf, 'Units', 'centimeters', 'Position', [2 scrsz(4)-8 7 5.25]);
set(gcf,'PaperPositionMode', 'auto', 'Color', 'w');
plot(Vvolsym,Uvolsym,'g.', 'MarkerSize', 16);
hold on;
plot(Vvolsym,Uvolsym,'k--');
xlabel('Volume [µm^3]');
ylabel('Strain Energy [pJ]');
if(~isempty(sevol_filename))
    if(strcmpi(sevol_filename((end-3):end), '.fig'))
        saveas(gcf,sevol_filename,'fig');
    end      
end
%
% Integrated strain energy around cell 
%
utmp=subvolume(seres.xg,seres.yg,seres.zg,seres.Ucdg_near,[-se_length_x_mum/2, se_length_x_mum/2, -se_length_y_mum/2, se_length_y_mum/2, -se_length_y_mum/2, se_length_y_mum/2]);
Ucell=sum(utmp(:))*cfg.gridsize^3;
if(~isempty(se_filename))
%     warning off MATLAB:xlswrite:AddSheet;
%     xlsout={'Directory', 'Integration box length x [µm]', 'y [µm]', 'z [µm]', 'Strain Energy [pJ]'; pwd sedensplot_length_x_mum sedensplot_length_y_mum sedensplot_length_z_mum Ucell};
%     xlswrite(se_filename, xlsout, 'CellSE');
    save(se_filename, 'se_length_x_mum', 'se_length_y_mum', 'se_length_z_mum', 'Ucell');
end
%
% Projected displacement figure
%
ddc=seres.d-repmat(mean(seres.d),size(seres.d,1),1);
[sbss_xg,sbss_yg,sbss_zg]=meshgrid( 20:20:1344, 20:20:1024, round(seres.cell_z_mum/(cfg.refractiveIdx*cfg.scale_Px2muM(3))-10):round(seres.cell_z_mum/(cfg.refractiveIdx*cfg.scale_Px2muM(3))+10) );
%
% scaling before interpolation  
%
sbss_xg_mum=cfg.scale_Px2muM(1).*sbss_xg;
sbss_yg_mum=cfg.scale_Px2muM(2).*sbss_yg;
sbss_zg_mum=cfg.refractiveIdx*cfg.scale_Px2muM(3).*sbss_zg;    
sbss_dxg_mum=griddata(cfg.scale_Px2muM(1).*seres.x(:,1), cfg.scale_Px2muM(2).*seres.x(:,2), cfg.refractiveIdx*cfg.scale_Px2muM(3).*seres.x(:,3), cfg.scale_Px2muM(1).*ddc(:,1), sbss_xg_mum, sbss_yg_mum, sbss_zg_mum);
sbss_dyg_mum=griddata(cfg.scale_Px2muM(1).*seres.x(:,1), cfg.scale_Px2muM(2).*seres.x(:,2), cfg.refractiveIdx*cfg.scale_Px2muM(3).*seres.x(:,3), cfg.scale_Px2muM(2).*ddc(:,2), sbss_xg_mum, sbss_yg_mum, sbss_zg_mum);
sbss_dzg_mum=griddata(cfg.scale_Px2muM(1).*seres.x(:,1), cfg.scale_Px2muM(2).*seres.x(:,2), cfg.refractiveIdx*cfg.scale_Px2muM(3).*seres.x(:,3), cfg.refractiveIdx*cfg.scale_Px2muM(3).*ddc(:,3), sbss_xg_mum, sbss_yg_mum, sbss_zg_mum);
%
sbss_dxg=(1/cfg.scale_Px2muM(1)).*sbss_dxg_mum;
sbss_dyg=(1/cfg.scale_Px2muM(2)).*sbss_dyg_mum;
sbss_dzg=(1/(cfg.refractiveIdx*cfg.scale_Px2muM(3))).*sbss_dzg_mum;
%
sbss_dxg(isnan(sbss_dxg(:)))=0;
sbss_dyg(isnan(sbss_dyg(:)))=0;
sbss_dzg(isnan(sbss_dzg(:)))=0;
%
dc=mean([sbss_dxg(:) sbss_dyg(:) sbss_dzg(:)]);
dxg=smooth3(sbss_dxg-dc(1));
dyg=smooth3(sbss_dyg-dc(2));
dzg=smooth3(sbss_dzg-dc(3));
xg=sbss_xg-seres.cell_x_mum/cfg.scale_Px2muM(1);
yg=sbss_yg-seres.cell_y_mum/cfg.scale_Px2muM(2);
zg=sbss_zg-seres.cell_z_mum/(cfg.refractiveIdx*cfg.scale_Px2muM(3));
%
pdx=round((celldisplot_range_x_px-1)/2);
pdy=round((celldisplot_range_y_px-1)/2);
%
[cx,cy,cz,dx]=subvolume(xg,yg,zg,dxg,[-pdx,pdx,-pdy,pdy,-10,10]);
[cx,cy,cz,dy]=subvolume(xg,yg,zg,dyg,[-pdx,pdx,-pdy,pdy,-10,10]);
[cx,cy,cz,dz]=subvolume(xg,yg,zg,dzg,[-pdx,pdx,-pdy,pdy,-10,10]);
%
if(exist('CELLframe_adjusted.jpg', 'file')==2)
    bf=double(imread('CELLframe_adjusted.jpg'));
    bf=bf(:,:,1);
else
    bf=double(imread(getBFfileName(round(seres.cell_z_mum/(cfg.refractiveIdx*cfg.scale_Px2muM(3))),cfg.t_deformed))); 
end
figure();
set(gcf,'Position', [100 100 401 401]);
ih=imagesc(bf( round(seres.cell_y_mum/cfg.scale_Px2muM(2)-pdy):round(seres.cell_y_mum/cfg.scale_Px2muM(2)+pdy), round(seres.cell_x_mum/cfg.scale_Px2muM(1)-pdx):round(seres.cell_x_mum/cfg.scale_Px2muM(1)+pdx) ));
axis square;
axis off;
colormap(gray);
set(gca,'Position', [0 0 1 1])
set(gca, 'YDir', 'normal');
hold on;
set(gcf,'PaperPositionMode', 'auto', 'Color', 'w');
if(celldisplot_scalef==0)
    quiver(cx(:,:,11)+pdx,cy(:,:,11)+pdy,dx(:,:,11),dy(:,:,11),'w');    
else
    quiver(cx(:,:,11)+pdx,cy(:,:,11)+pdy,celldisplot_scalef.*dx(:,:,11),celldisplot_scalef.*dy(:,:,11),0,'w');    
end
if(~isempty(celldisplot_filename))
    if(strcmpi(celldisplot_filename((end-3):end), '.fig'))
        saveas(gcf,celldisplot_filename,'fig');
    end      
end