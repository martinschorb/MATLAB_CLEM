function [tfmpts,hm_insert]  = martin_LM2HMauto(pts,inlm,inhm,crop,magx,interactive,slicenums)

% % version MartinSchorb 160627
% % Copyright EMBL 2011-2016, 
% %
%      This file is part of EMBL script for high-accuracy CLEM.
%  
%      EMBL script for high-accuracy CLEM is free software: you can redistribute it and/or modify
%      it under the terms of the GNU General Public License as published by
%      the Free Software Foundation, either version 3 of the License, or
%      (at your option) any later version.
%  
%      EMBL script for high-accuracy CLEM is distributed in the hope that it will be useful,
%      but WITHOUT ANY WARRANTY; without even the implied warranty of
%      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%      GNU General Public License for more details.
%  
%      You should have received a copy of the GNU General Public License
%      along with EMBL script for high-accuracy CLEM.  If not, see <http://www.gnu.org/licenses/>.
%  

%
%usage is martin_LM2HMauto('coordinate list','lowmag_image','highmag_image', mag_factor)
%
% designed for correlating em images of different magnifications.
%
%
%outputs transformed coordinates
%

if nargin<7
    slicenums =[0 0];
end


if nargin<6
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

