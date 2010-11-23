%% corr_init()
%
% 
% Version Martin Schorb 101123
% Initializes pathnames and parameters for LM/EM Correlation script
% martin_correlate.m



% Location to search for fiducial coordinate file - *.pickspots.mat

loc_pickspots = '/struct/briggs/wanda/DataLightMicroscopy';


% Location to search for shift correction coordinate file - *.shiftcoos.mat

loc_shiftcoos = '/struct/briggs/wanda/DataLightMicroscopy';


% Location to search for Highmag Fiducial coordinate file - *.shiftcoos.mat

loc_hmcoos = '/struct/briggs/wanda/DataLightMicroscopy';




% other parameters that can be adjusted

% emboxsize=57; % must be odd number   -  size of box for EM-image subpixel fitting
fmboxsize=19; % must be odd number   -  size of box for bead-image subpixel fitting
imboxsize=19; % must be odd number   -  size of box for fluo-image subpixel fitting

accuracy=36;  % Size of prediction circle in pixel



