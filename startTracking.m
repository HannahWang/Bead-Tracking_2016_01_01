function startTracking
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
% Configuraton for Tracking algorithm resides in
% getStrainEnergyAssayConfigDefaults()
% Make adjustments in that file.
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% Don't edit below this line
% ------------------------------------------------------------------------

trackDir = uigetdir('', 'Choose directory for bead tracking.');

if(trackDir)

    cwd=pwd;
    cd(trackDir);
    cfg=getStrainEnergyAssayConfig();
    cd(cwd);
    
    if(exist(sprintf('%s',trackDir,filesep,cfg.startposfilename),'file'))
        cont=questdlg(sprintf('File %s already exists. Overwrite or use instead of finding beads?', cfg.startposfilename), 'Overwrite file?', 'Overwrite', 'Use', 'Use');
        if(isequal(cont,'Overwrite'))
            delete(sprintf('%s',trackDir,filesep,cfg.startposfilename));
        end
    end

    if(exist(sprintf('%s',trackDir,filesep,cfg.beadposfilename),'file'))
        cont=questdlg(sprintf('File %s already exists. Overwrite?', cfg.beadposfilename), 'Overwrite file?', 'Overwrite', 'Cancel', 'Cancel');
        if(isequal(cont,'Overwrite'))
            delete(sprintf('%s',trackDir,filesep,cfg.beadposfilename));
        else
            return;
        end
    end    
    
    trackBeads3D(cfg, trackDir);

end    
