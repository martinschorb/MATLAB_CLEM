%% corr_init()
%
% Version Martin Schorb 120725
% Initializes pathnames and parameters for LM/EM Correlation script
% martin_correlate.m

% Location to search for fiducial coordinate file - *.pickspots.mat

loc_pickspots = pwd;%'/struct/briggs/schorb/beadcoord';

% Location to search for shift correction coordinate file - *.shiftcoos.mat

loc_shiftcoos = pwd;%'/struct/briggs/schorb/';

% Location to search for Highmag Fiducial coordinate file - *.lmhmcoos.mat

loc_hmcoos = pwd;%'/struct/briggs/wanda/DataLightMicroscopy/';

% Flip fluorescence images (different grid orientation in LM and EM)
% flip = 1 if images should be flipped // flip = 0 if images are in same orientation

flip = 0;

% Adjust the contrast of display of the fluorescence images (blue, green, red) (0 - auto;1 - open adjustment window)

contr_fid = 1;
contr_poi = 1;
contr_other = 0;

% Skip the shift adjustment between channels (active if 1)

shift_skip = 0;

% Size of prediction circle in nm

accuracy=50;  

% write overlay images for high-mag correlation (files will be written if 1)

hm_overlays = 0;

% multiple spots of interest?

multispot = 1;

% subpixel localization of fluorophores enable/disable (Optimization Toolbox required!!)
% options:
% 0 - no fitting
% 1 - fit signal of interest only
% 2 - fit fiducials only
% 3 - fit both fiducials and signal of interest

gaussloc = 3;

% interactive mode for fitting (0 - inactive, 1 - active if bad fit detected, 2 - always active)

fit_interactive = 2;


%  ------------------------------------------------------------------------

% other parameters that can be adjusted, if unsure leave them as they are

% trafo = 'linear conformal';
trafo = 'affine';
% trafo = 'projective';

% emboxsize=57; % must be odd number   -  size of box for EM-image subpixel fitting
fmboxsize=7; % must be odd number  - size of box for bead-image subpixel fitting
imboxsize=7; % must be odd number  -  size of box for fluo-image subpixel fitting

fmfilter = 45;

pixelsize_lm = 3.8936; % pixel size lowmag tomogram in nm
pixelsize_hm = 0.3780 ; % pixel size highmag tomogram in nm
hmaccuracy=accuracy/pixelsize_hm;
minbeads = 3;

% check Toolbox

a=ver('optim');
if isempty(a)
    warning('Optimization Toolbox not found! Not using subpixel localization.');
    gaussloc = 0;
end
clear a

% have all init variables ready in a structure to pass to succeeding
% functions 

vars=who;
for i=1:length(vars)
    init.(vars{i})=evalin('caller', vars{i});
end
