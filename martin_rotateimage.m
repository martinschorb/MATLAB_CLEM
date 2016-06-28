function [outim,rotid] = martin_rotateimage(em,im)

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
% usage is [output_image2, rotation parameter] = martin_rotateimage(input_image1,input_image2-to be turned); 
% 
%
% graphical interface to turn images with respect to each other


% generate title
imname = inputname(2);

switch imname
    case 'em'
        imtitle='EM';
    case 'fm'
        imtitle='Fluorescent fiducial';
    case 'gm' 
        imtitle='GFP';
    case 'rm'
        imtitle='RFP';
    case 'im'
        imtitle='Fluorescence spot of interest'; 
    case 'om'
        imtitle='Fluorescence spot of interest - second channel';     
    otherwise
        imtitle='Fluorescence';
end

% im=uint16(im);
im_orig = im;
outim=im;
rotid=0;

hFig = figure('Toolbar','none',...
              'Menubar','none',...
              'IntegerHandle','off',...
              'NumberTitle','off',...
              'Name','EM image');
            set(hFig,  'Position',[0 160 500 500]);
          em_sh=imshow(em);

szim=size(im);
          
hFig2 = figure('Toolbar','none',...
              'Menubar','none',...    
              'IntegerHandle','off',...
              'NumberTitle','off',...
              'Name',[imtitle,' image - rotation']);
set(hFig2,'Position',[521 160  800 800]);
          
im_sh=imshow(im);
hold on
hSp = imscrollpanel(hFig2,im_sh);
api = iptgetapi(hSp);
api.setMagnification(2);

% set(hSp,'Units','normalized','Position',[0.1 0.1 0.9 0.9])


h_rotbutton = uicontrol('Parent',hFig2,'Style','PushButton','Units','Normalized','Position',[0.5 0.02 0.2 0.03],'Callback',@rotbutton,'String','Rotate counterclockwise');

h_closebutton = uicontrol('Parent',hFig2,'Style','PushButton','Units','Normalized','Position',[0.72 0.02 .1 0.03],'Callback',@closebutton,'String','Done');

uiwait
uiresume
close all

rotid=rem(rotid,4);

% -----------------------------------

function rotbutton(h_rotbutton,event)
     im=rot90(im);

     delete im_sh
     delete(hSp)
     im_sh=imshow(im);
     hSp = imscrollpanel(hFig2,im_sh);
    api = iptgetapi(hSp);
    api.setMagnification(2);
     szim=size(im);
     set(hFig2,'Position',[521 160 800 800]);
%      if rem(rotid,2)
%          set(hSp,'Position',[0 0 800 800]); %      set(hFig2,'Position',[2000 1801 szim szim+50]);
%      end
     
    
     drawnow     
     outim=im;
     rotid=rotid + 1;
end

% -----------------------------------

function closebutton(h_closebutton,event)
     outim=im;close all;
end


end