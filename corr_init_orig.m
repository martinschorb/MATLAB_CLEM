% corr_init()
% ----------------------------------------
% version MartinSchorb 130312
% Copyright EMBL 2013, All rights reserved
% 
% Initializes pathnames and parameters for LM/EM Correlation script
% martin_correlate.m


% ---------  replace > pwd < by the directory of choice in parentheses ----

% Location to search for fiducial coordinate file - *.pickspots.mat 

loc_pickspots = pwd; %    '/directory/to/correlation';

% Location to search for shift correction coordinate file - *.shiftcoos.mat

loc_shiftcoos = pwd; %    '/directory/to/correlation';

% Location to search for Highmag Fiducial coordinate file - *.lmhmcoos.mat

loc_hmcoos = pwd; %    '/directory/to/HMcorrelation/';

% Flip fluorescence images (different grid orientation in LM and EM)
% flip = 1 if images should be flipped // flip = 0 if images are in same orientation

flip = 0;

% automatic hm correlation (hmauto=2 means interactive)

hmauto = 2; hmcrop = 1; magx=0;

% Adjust the contrast of display of the fluorescence images (blue, green, red) (0 - auto;1 - open adjustment window)

contr_fid = 1;
contr_poi = 0;
% obsolete contr_other = 0;

% Skip the shift adjustment between channels (active if 1)

shift_skip = 0;

% Size of prediction circle in nanometers

accuracy=50;  

% write overlay images for high-mag correlation (files will be written if 1)

hm_overlays = 0;

% multiple spots of interest?

multispot = 0;

% subpixel localization of fluorophores enable/disable (Optimization Toolbox required!!)
% options:
% 0 - no fitting
% 1 - fit signal of interest only
% 2 - fit fiducials only
% 3 - fit both fiducials and signal of interest

gaussloc = 3;

% interactive mode for fitting (0 - inactive, 1 - active, 2 - always active)

fit_interactive = 2;

%  ------------------------------------------------------------------------

% other parameters that can be adjusted, if unsure leave them as they are

trafo = 'linear conformal';
% trafo = 'affine';
% trafo = 'projective';

fmboxsize=11; % must be odd number   -  size of box for bead-image subpixel fitting
imboxsize=11; % must be odd number   -  size of box for fluo-image subpixel fitting

fmfilter = 45;

pixelsize_lm = 5.068; % pixel size lowmag tomogram in nm
pixelsize_hm = 1.18 ; % pixel size highmag tomogram in nm
hmaccuracy=accuracy*pixelsize_lm/pixelsize_hm;
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