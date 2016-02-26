function startStrainEnergyComputation
% Starts the computation of strain energy.
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
% Configuraton for strain energy algorithm resides in
% getStrainEnergyAssayConfigDefaults()
% Make adjustments in that file.
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
% Don't edit below this line
% ------------------------------------------------------------------------

seDir = uigetdir('', 'Choose directory for strain energy calculation.');

if(seDir)

    cwd=pwd;
    cd(seDir);
    cfg=getStrainEnergyAssayConfig();
    cd(cwd);
    
    if(~exist(sprintf('%s',seDir,filesep,cfg.beadposfilename),'file'))
        errordlg(sprintf('Tracking results in %s are missing. Run bead tracking first.',cfg.beadposfilename), 'Tracking results missing.');
        return;
    end

    if(exist(sprintf('%s',seDir,filesep,cfg.resultsfilename),'file'))
        cont=questdlg(sprintf('File %s already exists. Overwrite?', cfg.resultsfilename), 'Overwrite file?', 'Overwrite', 'Cancel', 'Cancel');
        if(isequal(cont,'Overwrite'))
            delete(sprintf('%s',seDir,filesep,cfg.resultsfilename));
        else
            return;
        end
    end

    strainEnergyComputation(cfg, seDir);

end
