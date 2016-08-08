function [im,head]=martin_loadim(varargin)

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


% opens EM image files

if nargin==0
        [in,indir] = uigetfile({'*.tif;*.mrc;*.st;*.rec','Image files';'*.*',  'All Files (*.*)'},'select image file',pwd);
        in=[indir,in];
else
    in=varargin{1};
end
    

if nargin==2   
    v=varargin{2};
else
    v=0;
end


if ~ischar(in) && ismatrix(in)
    im=in;
    head=[];
    
elseif exist(in) == 2
    if ~isempty([strfind(in,'.mrc'),strfind(in,'.st'),strfind(in,'.rec')])
        imf=tom_mrcread(in);
        head=imf.Header;
        if length(head.Size)>2
        if head.Size(3)>1
            if v==0
                kk=1:head.Size(3);
                [v,s]=listdlg('PromptString', 'Select slice','SelectionMode', 'single', 'ListString', num2str(kk') , 'Name' , imf);
            end
            im=rot90((imf.Value(:,:,v)));
        else    
            im=rot90(imf.Value);
        end
        end
    elseif ~isempty(strfind(in,'tif'))
        
        head=imfinfo(in);
        
        if strcmp(head(1).PhotometricInterpretation,'RGB')
            im=imread(in);
        elseif v==0;
            kk=size(head);v=1;
            if kk(1)>1
            [v,s]=listdlg('PromptString', 'Select slice','SelectionMode', 'single', 'ListString', num2str((1:kk(1))') , 'Name' , in);
            end
            im=imread(in,v);
        else
            im=imread(in,v);
        end
        head=[];
    else
        a=msgbox('Please provide either an image file in tif or mrc format!','Error','modal');
        error('Please provide either an image file in tif or mrc format!')   ;
    end
else
    a=msgbox('Please provide either an image matrix or a valid file name!','Error','modal');
    error('Please provide either an image matrix or a valid file name!')    ;
end

