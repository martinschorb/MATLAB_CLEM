%% corr_init()
%
% 
% Version Martin Schorb 120222
% Initializes pathnames and parameters for LM/EM Correlation script
% martin_correlate.m


% Location to search for fiducial coordinate file - *.pickspots.mat

loc_pickspots = '/struct/briggs/wanda/DataLightMicroscopy';

% Location to search for shift correction coordinate file - *.shiftcoos.mat

loc_shiftcoos = '/struct/briggs/wanda/DataLightMicroscopy';

% Location to search for Highmag Fiducial coordinate file - *.lmhmcoos.mat

loc_hmcoos = '/struct/briggs/wanda/DataLightMicroscopy';

% Flip fluorescence images (different grid orientation in LM and EM)
% flip = 1 if images should be flipped // flip = 0 if images are in same orientation

flip = 0;

% Adjust the contrast of display of the fluorescence images (blue, green, red) (0 - auto;1 - open adjustment window)

contr_b = 1;
contr_g = 0;
contr_r = 0;

% Skip the shift adjustment between channels (active if 1)

shift_skip = 0;

% Size of prediction circle in pixel

accuracy=36;  

% write overlay images for high-mag correlation (files will be written if 1)

hm_overlays = 0;

% multiple spots of interest?

multispot=0;

%  ------------------------------------------------------------------------
%  

% other parameters that can be adjusted, if unsure leave them as they are

% emboxsize=57; % must be odd number   -  size of box for EM-image subpixel fitting
fmboxsize=19; % must be odd number   -  size of box for bead-image subpixel fitting
imboxsize=19; % must be odd number   -  size of box for fluo-image subpixel fitting

pixelsize_lm = 5.068; % pixel size lowmag tomogram in nm
pixelsize_hm = 1.18 ; % pixel size highmag tomogram in nm

hmaccuracy=accuracy*pixelsize_lm/pixelsize_hm;



