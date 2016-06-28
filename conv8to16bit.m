function im=conv8to16bit(im)

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


%converts an image from 8bit to 16 bit


s=size(im,3);

im1=uint16(im);
im2=uint16(zeros(size(im)));

if isa(im,'uint8')
    for i=1:s
        im2(:,:,i)=256*im1(:,:,i);
    end
    im=im2;
else
<<<<<<< HEAD
%     disp(' image format')
=======
%     disp('ERROR: wrong image format')
>>>>>>> master
end
