function fname=getFLfileName(z,t)
%
% function returning the filename of a fluorecent slice at depth z below
% the top surface (counting from 0) at timestep t.
% -------------------------------------------------------------------------
% This file is part of the method published in
%
% Koch TM, Münster S, Bonakdar N, Butler JP, Fabry B (2012) 3D traction 
% forces in cancer cell invasion. PLoS ONE
%
% If you use any part of it, please cite this paper.
% -------------------------------------------------------------------------

% ------------------------------------------------------------------------
% adjust to your naming convention
% ------------------------------------------------------------------------
fname = sprintf('FLFRAME/FLframe_t%06i_%04i.tif', t, z);

end


