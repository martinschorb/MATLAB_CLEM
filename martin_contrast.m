function outim = martin_contrast(im)

% version MartinSchorb 120312
% 
% usage is output_image = martin_conrast(input_image); 
% 
%
% graphical interface to adjust the contrast of images


% generate title
imname = inputname(1);

switch imname
    case 'em'
        imtitle='EM';
    case 'fm'
        imtitle='Fluorescent fiducial';
    case 'gm' 
        imtitle='GFP';
    case 'rm'
        imtitle='RFP';
    otherwise
        imtitle='Fluorescence';
end

outim=im;
im=double(im);
im_orig = im;
outim=uint16(imadjust(outim));


ev=0;

hFig = figure('MenuBar','none','Toolbar','none','NumberTitle','off','Name',[imtitle,' image - contrast adjustment']);
im_area=imshow(outim);hold on

hSp = imscrollpanel(hFig,im_area);
api = iptgetapi(hSp);
api.setMagnification(2);

hMagBox = immagbox(hFig,im_area);
set(hMagBox,'Position',[10 10 50 25]);
% hover = imoverview(im_area);
% set(hover,'MenuBar','none')

set(hSp,'Units','normalized','Position',[0 .1 1 0.9])
hOvPanel = imoverview(im_area);
% set(hOvPanel,'Units','Normalized','Position',[0 0.1 1 .3])

h_text =  uicontrol('Parent',hFig,'Style','text','Units','Normalized','Position',[0.05 0.03 0.1 .02],'String','Adjust Contrast:');
h_text2 =  uicontrol('Parent',hFig,'Style','text','Units','Normalized','Position',[0.05 0.01 0.1 .02],'String','Adjust Offset:');

h_slider = uicontrol('Parent',hFig,'Style','slider','Units','Normalized','Position',[0.2 0.03 0.2 0.015],'Callback',@contrast_slider,'Value', 0.5);
h_slider2 = uicontrol('Parent',hFig,'Style','slider','Units','Normalized','Position',[0.2 0.01 0.2 0.015],'Callback',@offset_slider,'Value', 0.5);

h_autobutton = uicontrol('Parent',hFig,'Style','PushButton','Units','Normalized','Position',[0.5 0.02 0.2 0.03],'Callback',@autobutton,'String','Auto Contrast/Reset');

h_closebutton = uicontrol('Parent',hFig,'Style','PushButton','Units','Normalized','Position',[0.72 0.02 .1 0.03],'Callback',@closebutton,'String','Done');

or_max=max(max(im_orig));
contr = 1-get(h_slider,'Value');
offset = get(h_slider2,'Value');

uiwait
uiresume
close gcf


% uiwait(hFig)
% -----------------------------------

function contrast_slider(h_slider,event)
     contr = 1-get(h_slider,'Value');
     im=im_orig.*(65535./(min([contr,1])*or_max))+(offset-0.5)*(65535);
%      im_contr=im;
     outim=uint16(im);
     im_area=imshow(outim);drawnow;
     ev=1;
end

% -----------------------------------

function offset_slider(h_slider2,event)
     offset = get(h_slider2,'Value');
     im=im_orig.*(65535./(min([contr,1])*or_max))+(offset-0.5)*(65535);
%      im_offs=im;
     outim=uint16(im);
     im_area=imshow(outim);drawnow;
     ev=1;
     
end


% -----------------------------------

function autobutton(h_autobutton,event)
     outim=imadjust(uint16(im_orig));
     im_area=imshow(outim);drawnow
     set(h_slider,'Value',0.5);
     set(h_slider2,'Value',0.5);
     ev=0;
end

% -----------------------------------

function closebutton(h_closebutton,event)
     close gcf;
     if ~ev
      outim=imadjust(outim);
     end
     
end


end


