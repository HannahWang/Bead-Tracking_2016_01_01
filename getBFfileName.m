function fname=getBFfileName(z,t)
%
% function returning the filename of a brightfield slice at depth z below
% the top surface (counting from 0) at timestep t.
% -------------------------------------------------------------------------
% This file is part of the method published in
%
% Koch TM, M�nster S, Bonakdar N, Butler JP, Fabry B (2012) 3D�traction�
% forces�in�cancer�cell�invasion. PLoS ONE
%
% If you use any part of it, please cite this paper.
% -------------------------------------------------------------------------

% ------------------------------------------------------------------------
% adjust to your naming convention
% ------------------------------------------------------------------------
if t==1
fname = sprintf('BFFRAME/BFframe_t%06i_%04i.tif', t, z);
elseif t==2
fname = sprintf('BFFRAME/BFframe_t%06i_%08i.tif', t, z);
end
end