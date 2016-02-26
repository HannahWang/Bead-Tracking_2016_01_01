function cfg=getStrainEnergyAssayConfig
%
% This function returns the configuration for all algorithms. 
% Adjust accordingly.
% -------------------------------------------------------------------------
% This file is part of the method published in
%
% Koch TM, Münster S, Bonakdar N, Butler JP, Fabry B (2012) 3D traction 
% forces in cancer cell invasion. PLoS ONE
%
% If you use any part of it, please cite this paper.
% -------------------------------------------------------------------------

% ------------------------------------------------------------------------
% Configuraton for Tracking algorithm
% =================================== 
%
% Please edit below.
% ------------------------------------------------------------------------

% Filename for storing initial bead positions (txt file)
cfg.startposfilename='start_pos.txt';

% Filename for storing tracked bead positions (txt file)
cfg.beadposfilename='bead_pos.txt';

% Size of subvolume containing a bead
cfg.xsize=7;
cfg.ysize=7;
cfg.zsize=7;

% Timestep to start tracking, usually 1
cfg.tstart=1;

% Timestep to end tracking, 0 to track all available
cfg.tend=0; 

% Handle to function returning the filename of a fluorescent slice at depth
% z (where 0 is the top surface) and time step t
cfg.FLfileName=@getFLfileName;

% CCD camera offset (darkness)
cfg.CCDoffset=20;

% Threshold for identification of bead candidates
cfg.beadCandidateThreshold=10;

% Threshold for identification of bad beads used to eliminate tracking
% errors. Good beads have a "maximum intensity border ratio" below this
% threshold, where the "maximum intensity border ratio" is calculated as
% followed:
% First the intensity of the bead containing subvolume is summed along the
% z direction giving a 2D image of size cfg.xsize/cfg.ysize.
% The ratio is then the maximum inensity of the boundary divided by the
% central intensity. Set to Inf to disable
cfg.goodBeadThreshold=0.7;

% ------------------------------------------------------------------------
% Configuraton for Strain Energy Computations
% =========================================== 
%
% Please edit below.
% ------------------------------------------------------------------------

% Timestep of deformed configuration
cfg.t_deformed=1;

% Timestep of undeformed (reference) configuration
cfg.t_undeformed=2;

% Conversion from Pixel to µm in x/y/z
cfg.scale_Px2muM=[645E-3 645E-3 527*3.8E-3];

% Refractive index of the medium, used to correct the z components
cfg.refractiveIdx=1.333;

% Elastic shear modulus of ECM in µN/µm²
cfg.Gmod=118E-6;

% Poisson's ratio of ECM 
cfg.poisson=0.35;

% Shape correction cutoff, between 0 (all elements) and 1 (no elements),
% 0.5 recommended
cfg.shapeCutoff=0.5;

% Cutoff for unrealistic strains. Set to Inf to switch off
cfg.strainCutoff=0.5;

% gridsize for nearest neighbor interpolation im µm
cfg.gridsize=5;

% bead displacement uncertainties in µm, used for noise correction
cfg.deltax=0.022;
cfg.deltay=0.022;
cfg.deltaz=0.130;

% file for storing results (mat file)
cfg.resultsfilename='strainEnergyResults.mat';


