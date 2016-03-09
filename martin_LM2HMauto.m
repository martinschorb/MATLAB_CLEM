function [tfmpts,hm_insert]  = martin_LM2HMauto(pts,inlm,inhm,crop,magx,interactive,slicenums)

% % version MartinSchorb 160107
% % Copyright EMBL 2016, All rights reserved
%
%usage is martin_LM2HMauto('coordinate list','lowmag_image','highmag_image', mag_factor)
%
% designed for correlating em images of different magnifications.
%
%
%outputs transformed coordinates
%

if nargin<6
    slices =[0 0];
end


if nargin<5
    interactive=1;
end

if nargin<5
    magx=0;
end

if nargin<4
    crop = 0;
    % region is in the center
end



[T,hm_insert,scale] = martin_emwarp(inlm,inhm,crop,magx,interactive,slicenums);



tfmpts=transformPointsInverse(T,pts);

