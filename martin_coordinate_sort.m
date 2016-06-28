function out=martin_coordinate_sort(coos,rotid,imsize)

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

% manages coordinates after #of 90degree rotations (rotid)

rotid=rem(rotid,4);
l=size(coos,1);

if ~rem(rotid,2)
    imsize=fliplr(imsize);
end

switch rotid
    case 0
        out=coos;
    case 1
        out(:,1)=repmat(imsize(2),[l,1])-coos(:,2)+1;
        out(:,2)=coos(:,1);
    case 2
        out(:,2)=repmat(imsize(2),[l,1])-coos(:,2)+1;
        out(:,1)=repmat(imsize(1),[l,1])-coos(:,1)+1;
    case 3
        out(:,1)=coos(:,2);
        out(:,2)=repmat(imsize(1),[l,1])-coos(:,1)+1;
end
      
        