function outim = martin_contrast(im)

im=uint16(im);
im_orig = im;
im_contr=im;
im_offs=im;
outim=im;

hFig = figure('MenuBar','none','Name','Fluorescence image contrast adjustment');
im_area=imshow(im);hold on
h_text =  uicontrol('Parent',hFig,'Style','text','Position',[10 30 200 25],'String','Adjust Contrast:');
h_text2 =  uicontrol('Parent',hFig,'Style','text','Position',[10 10 200 20],'String','Adjust Offset:');

h_slider = uicontrol('Parent',hFig,'Style','slider','Position',[350 35 350 20],'Callback',@contrast_slider,'Value', 0.5);
h_slider2 = uicontrol('Parent',hFig,'Style','slider','Position',[350 10 350 20],'Callback',@offset_slider,'Value', 0.5);

h_autobutton = uicontrol('Parent',hFig,'Style','PushButton','Position',[750 15 150 35],'Callback',@autobutton,'String','Auto Contrast/Reset');

h_closebutton = uicontrol('Parent',hFig,'Style','PushButton','Position',[950 15 120 35],'Callback',@closebutton,'String','Done');

uiwait(hFig)
% -----------------------------------

function contrast_slider(h_slider,event)
     contr = 1-get(h_slider,'Value');
     im=im_offs.*(65535./(min([contr,1])*max(max(im_offs))));
     im_contr=im;
     im_area=imshow(im);drawnow;
     outim=im;
end

% -----------------------------------

function offset_slider(h_slider2,event)
     offset = get(h_slider2,'Value');
     im=uint16(im_contr+(offset-0.5)*(65535));
     im_offs=im;
     im_area=imshow(im);drawnow;
     
end


% -----------------------------------

function autobutton(h_autobutton,event)
     contr = get(h_slider,'Value');
     im=imadjust(im_orig);
     im_area=imshow(im);drawnow
     set(h_slider,'Value',0.5);
     set(h_slider2,'Value',0.5);
     outim=im;
end

% -----------------------------------

function closebutton(h_closebutton,event)
     outim=im;close gcf;
end


end


