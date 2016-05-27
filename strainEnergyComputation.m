function strainEnergyComputation(cfg, targetDir)
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
if(~isempty(targetDir))
    cd(targetDir);
end

% ------------------------------------------------------------------------
% Init
% ------------------------------------------------------------------------
shapefile=sprintf('cell_shape_%06d.mat', cfg.t_deformed);
if(exist(shapefile, 'file')==2)
    cs=load(shapefile);
else
    cs.cell_x=0;
    cs.cell_y=0;
    cs.cell_z=0;
    cs.cell_theta=0;
end
   
% loading bead positions and calculate displacements
statusbar=waitbar(0,'Calculating bead displacements...');
bp=dlmread(cfg.beadposfilename,'\t');
gb=extractGoodBeads(bp, cfg.goodBeadThreshold);

% refraction correction
gb(:,5)=cfg.refractiveIdx.*gb(:,5);

[x d bn]=extractBeadDisplacements(gb, cfg.t_undeformed, cfg.t_deformed);
x_mum=convertPx2muM(x, cfg.scale_Px2muM);
d_mum=convertPx2muM(d, cfg.scale_Px2muM);

% number of beads
nb=size(x_mum,1);

% meshing
waitbar(0.025, statusbar, 'Generating element mesh...');
elements=tessellateBeadsWithAuxPoints(x_mum);

% rotate and recenter around cell
cell_x_mum=cfg.scale_Px2muM(1)*cs.cell_x;
cell_y_mum=cfg.scale_Px2muM(2)*cs.cell_y;
cell_z_mum=cfg.refractiveIdx*cfg.scale_Px2muM(3)*cs.cell_z;
x_mum_orig=x_mum;
d_mum_orig=d_mum;    
x_mum=x_mum_orig-repmat([cell_x_mum cell_y_mum cell_z_mum], nb,1);
if(cs.cell_theta~=0)
    R=rotateAboutZ(-cs.cell_theta);
    for i=1:nb
        x_mum(i,:)=( R*x_mum(i,:)' )';
        d_mum(i,:)=( R*d_mum(i,:)' )';
    end
end

% mesh quality 
mq=meshQuality(x_mum, elements);
elements(find(mq<cfg.shapeCutoff),:)=[];

% compute strain energy
waitbar(0.2, statusbar, 'Computing stiffness matrix...');
Emod=cfg.Gmod*2*(1+cfg.poisson);
[sm sdm Vrawel]=stiffnessMatrixPerElement(elements, x_mum, Emod, cfg.poisson);
waitbar(0.29, statusbar, 'Computing element strain energies...');
Urawel=strainEnergyPerElement(sm, elements, d_mum);
[StressRaw, StrainRaw]=stressStrainPerElement(sdm, elements, d_mum, Emod, cfg.poisson);
%
% filter unrealistic strain
% 
[principalStrainDirections principalStrainValues]=diagonalizeStrainTensor(StrainRaw);
idx_min=find(principalStrainValues(:,1)<=-cfg.strainCutoff | principalStrainValues(:,1)>=cfg.strainCutoff);
idx_mid=find(principalStrainValues(:,2)<=-cfg.strainCutoff | principalStrainValues(:,2)>=cfg.strainCutoff);
idx_max=find(principalStrainValues(:,3)<=-cfg.strainCutoff | principalStrainValues(:,3)>=cfg.strainCutoff);
idx_strainfilter=unique([idx_min; idx_mid; idx_max]);
%
% interpolation
%
waitbar(0.37, statusbar, 'Interpolating strain energy densities...');
elsf=elements; % strain filtered elements
elsf(idx_strainfilter,:)=[];
Urawel(idx_strainfilter)=[];
Vrawel(idx_strainfilter)=[];
Urawdens=Urawel./Vrawel;
tetcentroids=[ 1/4*( x_mum(elsf(:,1),1) + x_mum(elsf(:,2),1) + x_mum(elsf(:,3),1) + x_mum(elsf(:,4),1) ) ...
               1/4*( x_mum(elsf(:,1),2) + x_mum(elsf(:,2),2) + x_mum(elsf(:,3),2) + x_mum(elsf(:,4),2) ) ...
               1/4*( x_mum(elsf(:,1),3) + x_mum(elsf(:,2),3) + x_mum(elsf(:,3),3) + x_mum(elsf(:,4),3) ) ];
%
xtmp=x_mum(unique(elsf(:)), :);
[xg yg zg]=meshgrid( (min(xtmp(:,1))):cfg.gridsize:(max(xtmp(:,1))), (min(xtmp(:,2))):cfg.gridsize:(max(xtmp(:,2))), (min(xtmp(:,3))):cfg.gridsize:(max(xtmp(:,3))) );
clear xtmp;
nearest_element=griddata(tetcentroids(:,1), tetcentroids(:,2), tetcentroids(:,3), 1:size(elsf,1), xg, yg, zg, 'nearest');
Urawdensg_near=reshape(Urawdens(nearest_element(:)), size(xg));
%
% noise correction
%
waitbar(0.92, statusbar, 'Noise correction...');
dU=sePertubationExpectedValue(sm,cfg.deltax,cfg.deltay,cfg.deltaz);
dU(idx_strainfilter)=[];
Ucel=Urawel-dU;
Ucel(find(Ucel<0))=0;
Ucd=Ucel./Vrawel;
Ucdg_near=reshape(Ucd(nearest_element(:)), size(xg));
%
% saving
%
clear sm sdm;
waitbar(0.94, statusbar, 'Saving results...');
save(cfg.resultsfilename);
close(statusbar);
% ------------------------------------------------------------------------
% Subfunctions
% ------------------------------------------------------------------------

function gb=extractGoodBeads(bp, threshold)
% 
% Throw out invalid beads (e.g. beads drifting out of field of view, beads
% above the cfg.goodBeadThreshold)
%   
    %
    % eliminate beads above goodBeadThreshold
    %
    crit=(bp(:,8)>threshold);        
    bp(crit,6)=0; 
    %
    % sort bead numbers
    %
    tsteps=sort(unique(bp(:,1)));
    for t=tsteps'
        idx=(bp(:,1)==t);
        bp_tmp=bp(idx,:);
        [dummy idx_tmp]=sort(bp_tmp(:,2));
        bp(idx,:)=bp_tmp(idx_tmp,:);
    end
    %
    % check if we have the same beads throughout each timestep    
    %
    ok=true;
    bn1=bp((bp(:,1)==tsteps(1)), 2);
    valid=bp((bp(:,1)==tsteps(1)), 6);
    for i=2:length(tsteps)
        if(~isequal(bn1, bp((bp(:,1)==tsteps(i)), 2)))
            ok=false;
        else
            valid=valid & bp((bp(:,1)==tsteps(i)), 6);
        end
    end
    if(ok)
        bp(:,6)=repmat(valid,length(tsteps),1);
    else
        warning('Bead numbers not consistent at all times. Finding common beads of all timesteps');
        %
        bpv=bp(bp(:,6)==1,:);
        bni=intersect(bpv(bpv(:,1)==tsteps(1),2), bpv(bpv(:,1)==tsteps(2),2));
        for i=3:length(tsteps)
            bni=intersect(bni, bpv(bpv(:,1)==tsteps(i),2));
        end
        warning('Found %d of %d initial beads', length(bni), length(bpv(bpv(:,1)==tsteps(1),2)));
        bp=bp(ismember(bp(:,2), bni),:);        
    end
    %
    % good beads are those still having valid bit set
    %
    gb=bp( (bp(:,6)==1),:);
    %
    % remove double entries
    %
    gb=removeDoubleBeads(gb);
end

function rb=removeDoubleBeads(bp)
%
% removes beads that occupy same position.
%
    t1=bp(1,1);
    t2=bp(end,1);
    for t=t1:t2    
        idx=find(bp(:,1)==t);
        x=bp(idx,3:5);
        m=[];
        [y m]=unique(x,'rows');
        bn=bp(idx,2);
        bn(m)=[];    
        for i=1:length(bn)
             idx=find(bp(:,2)==bn(i));
             bp(idx,6)=0;
        end
    end
    rb=bp(find(bp(:,6)==1),:);
end

function [xdata displacements bead_numbers]=extractBeadDisplacements(bp, t1, t2)
%
% Calculates displacements of beads at positions xdata between time t1 and
% t2.
%
    idx1=find(bp(:,1)==t1);
    idx2=find(bp(:,1)==t2);
    if(length(idx1)~=length(idx2))
        error('Not the same number of beads at times t1 and t2.');
    end

    bn1=bp(idx1,2);
    bn2=bp(idx2,2);
    if(~isequal(bn1, bn2))
        error('Not the same beads at times t1 and t2 or not in the same order.');
    end

    xdata=bp(idx1, 3:5);
    displacements=bp(idx2,3:5)-xdata;
    bead_numbers=bn1;
    disp(bn1);
end

function converted=convertPx2muM(xdata, scale)
%
% converts position vectors [x y z; ...] from pixel to µm given the scale,
% i.e. conversion factors from pixel to µm in [fx fy fz].
%
    converted(:,1)=scale(1).*xdata(:,1);
    converted(:,2)=scale(2).*xdata(:,2);
    converted(:,3)=scale(3).*xdata(:,3);
end

function [elements elements2 x2]=tessellateBeadsWithAuxPoints(x_mum)
%
% mesh generation with auxiliary points around the boundary surfaces of
% bead beads positions. those aux points are removed after mesh generation.
%
    gridsize=25;
    griddisp=5;

    x_min=min(x_mum(:,1))-griddisp;
    x_max=max(x_mum(:,1))+4.1*griddisp;
    y_min=min(x_mum(:,2))-griddisp;
    y_max=max(x_mum(:,2))+4.1*griddisp;
    z_min=min(x_mum(:,3))-griddisp;
    z_max=max(x_mum(:,3))+4.1*griddisp;

    [xg yg zg]=meshgrid( x_min:gridsize:x_max, y_min:gridsize:y_max, z_min:gridsize:z_max );

    x_tra=unique([xg(:) yg(:) zg(:)], 'rows');
    x_max=max(xg(:));
    y_max=max(yg(:));
    z_max=max(zg(:));
    x_tra=x_tra(find(x_tra(:,1)==x_min | x_tra(:,1)==x_max | x_tra(:,2)==y_min | x_tra(:,2)==y_max | x_tra(:,3)==z_min | x_tra(:,3)==z_max),:);

    x2=[x_mum; x_tra];
    elements2=delaunay(x2(:,1), x2(:,2), x2(:,3));

    N=length(x_mum);
    idx=find(elements2(:,1)>N | elements2(:,2)>N | elements2(:,3)>N | elements2(:,4)>N);
    elements=elements2;
    elements(idx,:)=[];
    %
    % Checking correct order of nodes
    %
    wrong_order=0;
    for j=1:size(elements,1)
        if(det([1,1,1,1;x_mum(elements(j,:),:)'])<0)
            elements(j,:)=elements(j,[1 3 2 4]);
            wrong_order=wrong_order+1;
        end
    end
end

function shape_quality=meshQuality(coordinates, elements)
%
% mesh quality of tetrahedral meshes
%
%
% 1 -> regular tetrahedron, best
% <1 -> decreasing quality
% <0 -> negative Jacobian (already taken care of with correct node numbering)
    shape_quality=zeros(size(elements,1),1);
    for j=1:size(elements,1)
        shape_quality(j) = ...
                            (6*sqrt(2)*det([1,1,1,1;coordinates(elements(j,:),:)'])) / ...
                            (sum((coordinates(elements(j,1),:)-coordinates(elements(j,2),:)).^2)^(3/2) + ...
                             sum((coordinates(elements(j,2),:)-coordinates(elements(j,3),:)).^2)^(3/2) + ...
                             sum((coordinates(elements(j,3),:)-coordinates(elements(j,1),:)).^2)^(3/2) + ...
                             sum((coordinates(elements(j,1),:)-coordinates(elements(j,4),:)).^2)^(3/2) + ...
                             sum((coordinates(elements(j,2),:)-coordinates(elements(j,4),:)).^2)^(3/2) + ...
                             sum((coordinates(elements(j,3),:)-coordinates(elements(j,4),:)).^2)^(3/2)); 
    end
end

function U=rotateAboutZ(alpha)
% gives the matrix corresponding to a rotation of alpha around Z axis 
%
    U=[ cos(alpha) -sin(alpha) 0; ...
        sin(alpha)  cos(alpha) 0; ...
        0           0          1];
end

function [sm sdm Vel]=stiffnessMatrixPerElement(elements, x_mum, Emod, poisson)
% Calculates stiffness matrix sm, strain-displacement-matrix sdm and the
% element volume for each element, given elastic contants Young's modulus
% Emod and Poisson's ratio poisson.
    Nel=size(elements,1);
    sm=zeros(Nel, 144);
    sdm=zeros(Nel,72);
    Vel=zeros(Nel,1);

    dc=Emod/((1+poisson)*(1-2*poisson));
    D=[(1-poisson)*dc poisson*dc poisson*dc 0 0 0; ...
        poisson*dc (1-poisson)*dc poisson*dc 0 0 0; ...
        poisson*dc poisson*dc (1-poisson)*dc 0 0 0; ...
        0 0 0 (1-2*poisson)*(dc/2) 0 0; ...
        0 0 0 0 (1-2*poisson)*(dc/2) 0;
        0 0 0 0 0 (1-2*poisson)*(dc/2)];

    for n=1:Nel
        J=[x_mum(elements(n,2),:)-x_mum(elements(n,1),:); ...
           x_mum(elements(n,3),:)-x_mum(elements(n,1),:); ...
           x_mum(elements(n,4),:)-x_mum(elements(n,1),:)]';
        dJ=det(J);
        R=inv(J);
        B=[-(R(1,1)+R(2,1)+R(3,1)) 0 0 R(1,1) 0 0 R(2,1) 0 0 R(3,1) 0 0; ....
           0 -(R(1,2)+R(2,2)+R(3,2)) 0 0 R(1,2) 0 0 R(2,2) 0 0 R(3,2) 0; ...
           0 0 -(R(1,3)+R(2,3)+R(3,3)) 0 0 R(1,3) 0 0 R(2,3) 0 0 R(3,3); ...
           0 -(R(1,3)+R(2,3)+R(3,3)) -(R(1,2)+R(2,2)+R(3,2)) 0 R(1,3) R(1,2) 0 R(2,3) R(2,2) 0 R(3,3) R(3,2); ...
           -(R(1,3)+R(2,3)+R(3,3)) 0 -(R(1,1)+R(2,1)+R(3,1)) R(1,3) 0 R(1,1) R(2,3) 0 R(2,1) R(3,3) 0 R(3,1); ...
           -(R(1,2)+R(2,2)+R(3,2)) -(R(1,1)+R(2,1)+R(3,1)) 0 R(1,2) R(1,1) 0 R(2,2) R(2,1) 0 R(3,2) R(3,1) 0];
       K=dJ/6* B'*D*B;
       %
       Vel(n)=dJ/6;
       sdm(n,:)=B(:)';
       sm(n,:)=K(:)';
    end
end

function U=strainEnergyPerElement(stiffnessMatrix, elements, d)
% calculates per element strain energy from stiffnessMatrix and
% displacements d.
    Nel=size(elements,1);
    U=zeros(Nel,1);
    for n=1:Nel
        K=reshape(stiffnessMatrix(n,:), 12,12);
        dd=reshape(d([elements(n,1) elements(n,2) elements(n,3) elements(n,4)],:)',12,1);
        U(n)=0.5.*dd'*K*dd;
    end
end

function [stress, strain]=stressStrainPerElement(strainDisplacementMatrix, elements, d, Emod, poisson)
% calculates per element stress and strain tensors from the
% strainDisplacementMatrix of tetrahedral elements with 
% displacements d. stress is computed by Hooke's law for an isotropic
% linear elastic material with Young's modulus Emod and Poisson's ratio
% poisson.
    dc=Emod/((1+poisson)*(1-2*poisson));
    D=[(1-poisson)*dc poisson*dc poisson*dc 0 0 0; ...
            poisson*dc (1-poisson)*dc poisson*dc 0 0 0; ...
            poisson*dc poisson*dc (1-poisson)*dc 0 0 0; ...
            0 0 0 (1-2*poisson)*(dc/2) 0 0; ...
            0 0 0 0 (1-2*poisson)*(dc/2) 0;
            0 0 0 0 0 (1-2*poisson)*(dc/2)];

    Nel=size(elements,1);
    strain=zeros(Nel,6);
    stress=zeros(Nel,6);
    for n=1:Nel
        B=reshape(strainDisplacementMatrix(n,:), 6,12);    
        dd=reshape(d([elements(n,1) elements(n,2) elements(n,3) elements(n,4)],:)',12,1);
        tmp=(B*dd)';
        strain(n,:)=tmp;
        stress(n,:)=(D*tmp')';
    end
end

function [eigvec eigval]=diagonalizeStrainTensor(T)
% returns eigenvectors and eigenvalues of strain tensor T with size T N x
% 6, i.e. T is symmetric with factor 2 in off-diagonal values.
    N=size(T,1);
    %
    % Diagonalization
    %
    eigvec=zeros(N, 3,3);
    eigval=zeros(N, 3);
    for i=1:N
        [eigvec(i,:,:) dummy]=eig( ...
                          [T(i,1) T(i,6)/2 T(i,5)/2; ...
                           T(i,6)/2 T(i,2) T(i,4)/2; ...
                           T(i,5)/2 T(i,4)/2 T(i,3)] );                                     
        eigval(i,:)=dummy([1 5 9]);
        if(~issorted(eigval(i,:))) error('Eigenvalues not sorted'); end
    end  
end

function dU=sePertubationExpectedValue(stiffnessMatrix, sig_x, sig_y, sig_z)
% calculates <U(u+du)>-<U(u)>.
    Nel=size(stiffnessMatrix,1);
    dU=zeros(Nel,1);
    sig=repmat([sig_x sig_y sig_z]', 4,1).^2;
    for n=1:Nel
        K=stiffnessMatrix(n, 1:13:end);    
        dU(n)=0.5*K*sig;
    end
end

% ------------------------------------------------------------------------
% End
% ------------------------------------------------------------------------
% return to prev dir
cd(saveDir);
end    