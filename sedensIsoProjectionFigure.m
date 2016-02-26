function sedensIsoProjectionFigure(cfg, sd, xg, yg, zg, region_length_x_mum, region_length_y_mum, region_length_z_mum, relativeIsoval, exportFilename)
%
% This function creates plots of strain energy density.
%
% -------------------------------------------------------------------------
% This file is part of the method published in
%
% Koch TM, Münster S, Bonakdar N, Butler JP, Fabry B (2012) 3D traction 
% forces in cancer cell invasion. PLoS ONE
%
% If you use any part of it, please cite this paper.
% -------------------------------------------------------------------------

rx=region_length_x_mum/2;
ry=region_length_y_mum/2;
rz=region_length_z_mum/2;
[nx,ny,nz,nsd]=subvolume(xg,yg,zg,sd,[-rx,rx,-ry,ry,-rz, rz]);
[nnx,nny,nnz,nnsd]=subvolume(nx,ny,nz,nsd,[min(nx(:))+20, max(nx(:))-20, min(ny(:))+20, max(ny(:))-20, min(nz(:))+0, max(nz(:))-0]);
planes=min(nsd(:)).*ones(size(nnsd));
planes(1)=1.001*min(nsd(:));
%
contour_x=0;
contour_y=0;
contour_z=0;
%
hz=slice(nx,ny,nz+abs(contour_z-min(nz(:))),nsd,[],[],max(nz(:)));
set(hz,'FaceColor','interp','EdgeColor','none');
xlim([min(nx(:)) max(nx(:))]);
ylim([min(ny(:)) max(ny(:))]);
zlim([min(nz(:)) max(nz(:))]);
hold on;
%
view(3);
%
hx=slice(nx-abs(contour_x-min(nx(:))),ny,nz,nsd,min(nx(:)),[],[]);
set(hx,'FaceColor','interp','EdgeColor','none');
hy=slice(nx,ny-abs(contour_y-min(ny(:))),nz,nsd,[],min(ny(:)),[]);
set(hy,'FaceColor','interp','EdgeColor','none');
%
set(gca,'XTickLabel',{});
set(gca,'YTickLabel',{});
set(gca,'ZTickLabel',{});
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'ZTick',[]);
set(gca,'ZDir', 'reverse');
daspect([1 1 1]);
%
az=148.5;
el=24;
view(az,el);
%
cmin=min(nsd(:)) ;
cmax=max(nsd(:));
caxis([cmin cmax]);
%
idx=find( xg(:)<-rx | xg(:)>rx | yg(:)<-ry | yg(:)>ry | zg(:)<-rz | zg(:)>rz );
sd2=sd;
sd2(idx)=0;
%
hold on;
piso=relativeIsoval*cmax; 
p=patch(isosurface(xg,yg,zg,sd2,piso));
isonormals(xg,yg,zg,sd2,p);
set(p,'FaceColor','cyan','EdgeColor','none', 'FaceAlpha', 0.7);
%
for az=-180:30:0
    for el=-100:10:100
        camlight(az,el);
    end
end
lighting none;
set(p,'FaceLighting','phong');
%
% adding cell outline
%
project_la=false;
if(exist('cell_shape_000001.mat', 'file')==2)
    load cell_shape_000001.mat;
    cell_x_mum=cfg.scale_Px2muM(1)*cell_x;
    cell_y_mum=cfg.scale_Px2muM(2)*cell_y;
    cell_z_mum=cfg.scale_Px2muM(3)*cell_z;
    cell_xi_mum=cfg.scale_Px2muM(1).*cell_xi-cell_x_mum;
    cell_yi_mum=cfg.scale_Px2muM(2).*cell_yi-cell_y_mum;
    R=[ cos(-cell_theta) -sin(-cell_theta) 0; ...
        sin(-cell_theta)  cos(-cell_theta) 0; ...
        0           0          1];    
    for i=1:size(cell_xi,1)
        cxtmp=R*[cell_xi_mum(i) cell_yi_mum(i) 1]';
        cell_xi_mum(i)=cxtmp(1);
        cell_yi_mum(i)=cxtmp(2);
    end
    plot3(cell_xi_mum, cell_yi_mum, (max(nz(:))-1).*ones(size(cell_xi_mum)), 'w', 'LineWidth', 1); % cell outline on xy plane
    if(exist('height_profile_000001.mat', 'file')==2)
        load height_profile_000001.mat;
        hprofile_mum=hprofile.*repmat(cfg.scale_Px2muM, size(hprofile,1), 1)-repmat([cell_x_mum cell_y_mum cell_z_mum], size(hprofile,1), 1);
        for i=1:size(hprofile,1)
            hprofile_mum(i,:)=(R*hprofile_mum(i,:)')';    
        end
        project_la=true;
    end 
end               
% long axis
if(project_la)
    m=(hprofile_mum(3,3)-hprofile_mum(1,3))/(hprofile_mum(3,1)-hprofile_mum(1,1));
    c=hprofile_mum(1,3)-m*hprofile_mum(1,1);
    z0=m*hprofile_mum(2,1)+c;
    c0=z0+1/m*hprofile_mum(2,1);
    l=2.5; % mum
    s=sqrt( l^2/(1/m^2-1) );
    ox1=hprofile_mum(2,1)-s;
    ox2=hprofile_mum(2,1)+s;
    oz1=-1/m*ox1+c0;
    oz2=-1/m*ox2+c0;
    xx=linspace(hprofile_mum(1,1), hprofile_mum(3,1), 21);
    % 
    zz1=spline(real([hprofile_mum(1,1) ox1 hprofile_mum(3,1)]),real([hprofile_mum(1,3) oz1 hprofile_mum(3,3)]), xx);
    zz2=spline(real([hprofile_mum(1,1) ox2 hprofile_mum(3,1)]),real([hprofile_mum(1,3) oz2 hprofile_mum(3,3)]), xx);
    plot3(xx,(min(ny(:))+1).*ones(size(xx)), zz1, 'w', 'LineWidth', 1);
    plot3(xx,(min(ny(:))+1).*ones(size(xx)), zz2, 'w', 'LineWidth', 1);
else
    xi=linspace(min(cell_xi_mum), max(cell_xi_mum), 3);
    xx=linspace(min(cell_xi_mum), max(cell_xi_mum), 21);
    zz1=spline(xi,[0 2.5 0], xx);
    zz2=spline(xi,[0 -2.5 0], xx);
    plot3(xx,(min(ny(:))+1).*ones(size(xx)), zz1, 'w', 'LineWidth', 1);
    plot3(xx,(min(ny(:))+1).*ones(size(xx)), zz2, 'w', 'LineWidth', 1);
end
% short axis
yi=linspace(min(cell_yi_mum), max(cell_yi_mum), 3);
yy=linspace(min(cell_yi_mum), max(cell_yi_mum), 21);
zz3=spline(yi,[0 2.5 0], yy);
zz4=spline(yi,[0 -2.5 0], yy);
plot3((min(nx(:))+1).*ones(size(yy)), yy, zz3, 'w', 'LineWidth', 1);
plot3((min(nx(:))+1).*ones(size(yy)), yy, zz4, 'w', 'LineWidth', 1);
%
% scale bar [µm]
%
sb_length_mum=50;
sb_offset_x=region_length_x_mum*0.05;
sb_offset_y=region_length_y_mum*0.05;
plot3([max(nx(:))-(sb_length_mum+sb_offset_x) max(nx(:))-sb_offset_x], [max(ny(:))-sb_offset_y max(ny(:))-sb_offset_y], [max(nz(:))-1 max(nz(:))-1], 'w', 'LineWidth', 3);
%
% exporting
%
if(~isempty(exportFilename))
    if(strcmpi(exportFilename((end-3):end), '.fig'))
        saveas(gcf,exportFilename,'fig');
    end
    if(strcmpi(exportFilename((end-3):end), '.png'))
        print('-dpng', '-opengl', '-r1200', exportFilename);
    end    
end