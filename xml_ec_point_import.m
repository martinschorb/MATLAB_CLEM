function pos=xml_ec_point_import(file)

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

if exist(file,'file')

st=xml2struct(file);

for i=1:length(st.root.rois.roi)    
    pos(i,1)=str2num(st.root.rois.roi{i}.position.pos_x.Text);
	pos(i,2)=str2num(st.root.rois.roi{i}.position.pos_y.Text);
end

else
    pos=[];
    
end