function [init,emf,fmf,imf,omf,outfile,fluorsel,omfluor,slices]=martin_correlate_init(varargin)

% % version MartinSchorb 130312
% % Copyright EMBL 2013, All rights reserved
% %
% =========================================================================
%         DO NOT MODIFY !!!!!  use init script to set up parameters!!!  
% =========================================================================    
% %
%     intitialisation of parameters for correlation using GUI 


close all
% need to run this outside because of nested functions
% if exist('corr_init')==2
%     corr_init();
% elseif exist('corr_init_orig')==2
%     corr_init_orig(); 
% else 
%     a=msgbox('No initialization script found!','Error','modal');uiwait(a);
%     a=msgbox('Please update algorithms!','Error','modal');uiwait(a);
%     return 
% end 



% init values

init = varargin{1};
outfile = varargin{2};

pxs = init.pixelsize_lm;
slices = struct;
% minbeads = 3;
slices.em = 1;
slices.fm = 1;
slices.im = 1;
slices.om = 1;
backcol1 = [0.85 0.85 0.85];
backcol2 = [0.95 0.95 0.95];


if nargin<3
    emf='';
    fmf='';
    imf='';
    omf='';
    fluorsel = 'GFP';   
elseif nargin<4
    error('wrong number of input arguments')
else
    fluorsel = varargin{3};
    emf = varargin{4};
    fmf = varargin{5};
    imf = varargin{6};
    omf='';
end

if nargin<8
    omfluor = 'RFP';
else
    omf = varargin{7};
    omfluor = varargin{8};
end

if nargin>8
    slices = varargin{9};
end

% GUI

f = figure('Visible','off','Position',[0,120,1000,1000],'NumberTitle','off','Menubar','none','Name','EM/FM Correlation - Initial Configuration');


t_outfile = uipanel('Title','output file name prefix','FontSize',12,'BackgroundColor',backcol1,...
          'Position',[0.05,0.9,0.4,0.08]);
c_outfile = uicontrol('Parent',t_outfile,'Style','edit','BackgroundColor',backcol2,...
           'Position',[10,10,380,40],'String',outfile); 
       
t_load = uipanel('Title','','FontSize',1,'BackgroundColor',backcol1,'BorderType','none',...
          'Position',[0.55,0.91,0.01,0.03]);       
go = uicontrol('Style','pushbutton','Parent',t_load,'BackgroundColor',backcol2,...
           'Position',[-1 -1 120 40],'String','load settings','Callback',{@load_Callback}); 
       
       
%    --- EM

t_emdir = uipanel('Title','EM image file','FontSize',12,'BackgroundColor',backcol1,...
           'Position',[0.05,0.8,0.4,0.08]);
ct_emdir = uicontrol('Parent',t_emdir,'Style','edit','BackgroundColor',backcol2,...
           'Position',[10,10,300,40],'String',emf);       
but_emdir = uicontrol('Parent',t_emdir,'Style','pushbutton','BackgroundColor',backcol2,...
           'Position',[315,10,65,40],'String','Browse','Callback',{@but_emdir_Callback});         
     
t_emslice = uipanel('Title','Slice number','FontSize',10,'BackgroundColor',backcol1,...
           'Position',[0.61,0.8,0.1,0.08]);
ct_emslice = uicontrol('Parent',t_emslice,'Style','edit','BackgroundColor',backcol2,...
           'Position',[10,10,50,40],'String',slices.em);       
       
       
       
%    --- Fiducials

t_fmdir = uipanel('Title','Fiducial image file','FontSize',12,'BackgroundColor',backcol1,...
           'Position',[0.05,0.7,0.4,0.08]);
ct_fmdir = uicontrol('Parent',t_fmdir,'Style','edit','BackgroundColor',backcol2,...
           'Position',[10,10,300,40],'String',fmf);       
but_fmdir = uicontrol('Parent',t_fmdir,'Style','pushbutton','BackgroundColor',backcol2,...
           'Position',[315,10,65,40],'String','Browse','Callback',{@but_fmdir_Callback});   
       
       
g_flip = uibuttongroup('Title','FM images flipped','FontSize',12,'BackgroundColor',backcol1,...
          'Position',[0.45,0.75,0.15,0.08]);
flip0 = uicontrol('Style','Radio','String','no','BackgroundColor',backcol1,...
    'pos',[10 30 110 25],'parent',g_flip,'HandleVisibility','off');
flip1 = uicontrol('Style','Radio','String','yes','BackgroundColor',backcol1,...
    'pos',[10 5 110 25],'parent',g_flip);
 % Initialize button group properties. 
set(g_flip,'SelectionChangeFcn',@flipselcbk);      
    switch init.flip
        case 0
            set(g_flip,'SelectedObject',flip0);
        case 1
            set(g_flip,'SelectedObject',flip1);
    end    
       

t_fmslice = uipanel('Title','Slice number','FontSize',10,'BackgroundColor',backcol1,...
           'Position',[0.61,0.7,0.1,0.08]);
ct_fmslice = uicontrol('Parent',t_fmslice,'Style','edit','BackgroundColor',backcol2,...
           'Position',[10,10,50,40],'String',slices.fm);
    

       
              
%    ---  IM

t_imdir = uipanel('Title','Flourescence image of interest','FontSize',12,'BackgroundColor',backcol1,...
           'Position',[0.05,0.6,0.4,0.08]);
ct_imdir = uicontrol('Parent',t_imdir,'Style','edit','BackgroundColor',backcol2,...
           'Position',[10,10,300,40],'String',imf);       
but_imdir = uicontrol('Parent',t_imdir,'Style','pushbutton','BackgroundColor',backcol2,...
           'Position',[315,10,65,40],'String','Browse','Callback',{@but_imdir_Callback});          
       
   
g_im1 = uibuttongroup('Title','Fluorophore','FontSize',12,'BackgroundColor',backcol1,...
          'Position',[0.45,0.6,0.15,0.08],'SelectionChangeFcn',@im1selcbk);
imu0 = uicontrol('Style','Radio','String','GFP','BackgroundColor',backcol1,...
    'pos',[10 38 110 18],'parent',g_im1,'HandleVisibility','off');
imu1 = uicontrol('Style','Radio','String','RFP','BackgroundColor',backcol1,...
    'pos',[10 20 110 18],'parent',g_im1);
imu2 = uicontrol('Style','Radio','String','other','BackgroundColor',backcol1,...
    'pos',[10 2 110 18],'parent',g_im1);

t_imslice = uipanel('Title','Slice number','FontSize',10,'BackgroundColor',backcol1,...
           'Position',[0.61,0.6,0.1,0.08]);
ct_imslice = uicontrol('Parent',t_imslice,'Style','edit','BackgroundColor',backcol2,...
           'Position',[10,10,50,40],'String',slices.im);

       
%    ---

t_omdir = uipanel('Title','additional Flourescence image (optional)','FontSize',12,'BackgroundColor',backcol1,...
           'Position',[0.05,0.5,0.4,0.08]);
ct_omdir = uicontrol('Parent',t_omdir,'Style','edit','BackgroundColor',backcol2,...
           'Position',[10,10,300,40],'String',omf);       
but_omdir = uicontrol('Parent',t_omdir,'Style','pushbutton','BackgroundColor',backcol2,...
           'Position',[315,10,65,40],'String','Browse','Callback',{@but_omdir_Callback});          
       
   
g_im2 = uibuttongroup('Title','Fluorophore','FontSize',12,'BackgroundColor',backcol1,...
          'Position',[0.45,0.5,0.15,0.08],'SelectionChangeFcn',@im2selcbk);
omu0 = uicontrol('Style','Radio','String','GFP','BackgroundColor',backcol1,...
    'pos',[10 38 110 18],'parent',g_im2,'HandleVisibility','off');
omu1 = uicontrol('Style','Radio','String','RFP','BackgroundColor',backcol1,...
    'pos',[10 20 110 18],'parent',g_im2);
omu2 = uicontrol('Style','Radio','String','other','BackgroundColor',backcol1,...
    'pos',[10 2 110 18],'parent',g_im2);
     
switch omfluor
    case 'GFP'
        set(g_im2,'SelectedObject',omu0);
    case 'RFP'
        set(g_im2,'SelectedObject',omu1);
    otherwise
        set(g_im2,'SelectedObject',omu2);
end
    
switch fluorsel
        case 'GFP'
            set(g_im1,'SelectedObject',imu0);
            set(omu0,'Visible','off');
            set(omu1,'Visible','on');
            if strcmp(omfluor,'GFP')
                set(g_im2,'SelectedObject',omu2);
                omfluor='other';
            end
        case 'RFP'
            set(g_im1,'SelectedObject',imu1);
            set(omu0,'Visible','on');
            set(omu1,'Visible','off');
            if strcmp(omfluor,'RFP')
                set(g_im2,'SelectedObject',omu2);
                omfluor='other';
            end
    otherwise
            set(g_im1,'SelectedObject',omu2);
            set(omu0,'Visible','on');
            set(omu1,'Visible','on');
end
    


t_omslice = uipanel('Title','Slice number','FontSize',10,'BackgroundColor',backcol1,...
           'Position',[0.61,0.5,0.1,0.08]);
ct_omslice = uicontrol('Parent',t_omslice,'Style','edit','BackgroundColor',backcol2,...
           'Position',[10,10,50,40],'String',slices.om);

       
       
%   ---      
      
t_px = uipanel('Title','EM pixel size [nm]','FontSize',12,'BackgroundColor',backcol1,...
          'Position',[0.05,0.4,0.4,0.08]);
ct_px = uicontrol('Parent',t_px,'Style','edit','BackgroundColor',backcol2,...
          'Position',[10,10,300,40],'String',pxs); 

%   ---
    
g_trafo = uibuttongroup('Title','Transformation type to use','FontSize',12,'BackgroundColor',backcol1,...
          'Position',[0.05,0.25,0.4,0.13]);
u0 = uicontrol('Style','Radio','String','linear conformal (default)','BackgroundColor',backcol1,...
    'pos',[10 70 300 30],'parent',g_trafo,'HandleVisibility','off');
u1 = uicontrol('Style','Radio','String','affine','BackgroundColor',backcol1,...
    'pos',[10 40 300 30],'parent',g_trafo);
u2 = uicontrol('Style','Radio','String','projective','BackgroundColor',backcol1,...
    'pos',[10 10 300 30],'parent',g_trafo);
% Initialize button group properties. 
set(g_trafo,'SelectionChangeFcn',@selcbk);      
    switch init.trafo
        case 'linear conformal'
            set(g_trafo,'SelectedObject',u0);
        case 'affine'
            set(g_trafo,'SelectedObject',u1);
        case 'projective'
            set(g_trafo,'SelectedObject',u2);
    end




%   ---
t_bds = uipanel('Title','Minimum number of beads','FontSize',12,'BackgroundColor',backcol1,...
          'Position',[0.05,0.15,0.4,0.08]);
ct_bds = uicontrol('Parent',t_bds,'Style','edit','BackgroundColor',backcol2,...
           'Position',[160,10,40,40],'String',init.minbeads); 



%   ---
go = uicontrol('Style','pushbutton','BackgroundColor',backcol2,...
           'Position',[50 50 300 50],'String','Go','Callback',{@go_Callback}); 


% movegui(f,'center')      
set(f,'Visible','on');
%
 uiwait(f);

% -------------------------- Callback functions  ----------------------

% load

function load_Callback(source,eventdata)
    [matfile,matdir] = uigetfile('*.pickspots1.mat','select existing correlation');
    if ~isnumeric(matfile)
        if exist([matdir,matfile],'file')
            a = load([matdir,matfile]);
        else return
        end
    else
        return
    end
    
    evalfields(source,eventdata);
    close(f);
    
    if ~isfield(a,'fluorsel')
        a.fluorsel = 'GFP';
        a.omfluor = 'RFP';
    end
    
    if ~isfield(a,'imf')
        a.imf = a.gmf;
        a.omf = a.rmf;
        a.slices = slices;       
        
        if isfield(a,'medshift_GFP')
            a.fluorsel = 'GFP';
            a.omfluor = 'RFP';
        else
            a.fluorsel = 'RFP';
            a.omfluor = 'GFP';
        end
    end
          
    
    [init,emf,fmf,imf,omf,outfile,fluorsel,omfluor,slices]=martin_correlate_init(init,outfile,a.fluorsel,a.emf,a.fmf,a.imf,a.omf,a.omfluor,a.slices);
end


% em

function but_emdir_Callback(source,eventdata)
    emf = get(ct_emdir,'String');
    if ~isempty(emf)
        dirpos = strfind(emf,'/');
        emdir = emf(1:dirpos(end));
    else
        emdir=pwd;
    end
    
    [emf,emdir] = uigetfile('*.tif','select EM image file to correlate',emdir);
    emf = [emdir,emf];
    set(ct_emdir,'String',emf)
end


% Fiducials

function but_fmdir_Callback(source,eventdata)
    fmf = get(ct_fmdir,'String');
    if ~isempty(fmf)
        dirpos = strfind(fmf,'/');
        fmdir = emf(1:dirpos(end));
    else
        fmdir=pwd;
    end
   
    [fmf,fmdir] = uigetfile('*.tif','select Fiducial image file to correlate',fmdir);
    fmf = [fmdir,fmf];
    set(ct_fmdir,'String',fmf)
end



function  flipselcbk(source,eventdata)
    flip = get(eventdata.NewValue,'String');
    switch flip
        case 'no'
            init.flip=0;
        case 'yes'
            init.flip = 1;
    end
end



% im

function but_imdir_Callback(source,eventdata)
    imf = get(ct_imdir,'String');
    if ~isempty(imf)
        dirpos = strfind(imf,'/');
        imdir = imf(1:dirpos(end));
    else
        imdir=pwd;
    end
   
    [imf,imdir] = uigetfile('*.tif','select Fluorescence image of interest to correlate',imdir);
    imf = [imdir,imf];
    set(ct_imdir,'String',imf)
end



function  im1selcbk(source,eventdata)
    
    fluorsel = get(eventdata.NewValue,'String');

    switch fluorsel
        case 'GFP'
            set(omu0,'Visible','off');
            set(omu1,'Visible','on');
            if strcmp(omfluor,'GFP')
                set(g_im2,'SelectedObject',omu2);
                omfluor='other';
            end
        case 'RFP'
            set(omu0,'Visible','on');
            set(omu1,'Visible','off');
            if strcmp(omfluor,'RFP')
                set(g_im2,'SelectedObject',omu2);
                omfluor='other';
            end
        otherwise
            set(omu0,'Visible','on');
            set(omu1,'Visible','on');
    end
    
end

% om

function but_omdir_Callback(source,eventdata)
    omf = get(ct_omdir,'String');
    omf = get(ct_imdir,'String');
    if ~isempty(omf)
        dirpos = strfind(omf,'/');
        omdir = omf(1:dirpos(end));
    else
        omdir=pwd;
    end
    
    [omf,omdir] = uigetfile('*.tif','select Fluorescence image of interest to correlate',omdir);
    omf = [omdir,omf];
    set(ct_omdir,'String',omf)
end



function  im2selcbk(source,eventdata)
    omfluor = get(eventdata.NewValue,'String');
end


% ---



function  selcbk(source,eventdata)
    init.trafo = get(eventdata.NewValue,'String');
    switch init.trafo
        case 'linear conformal (default)'
            init.trafo = 'linear conformal';
            init.minbeads = 3;
            set(ct_bds,'String','3');
        case 'affine'
            init.minbeads = 3;
            set(ct_bds,'String','3');
        case 'projective'
            init.minbeads = 4;
            set(ct_bds,'String','4')
    end
end


function go_Callback(source,eventdata)
    evalfields(source,eventdata);
    close(f);
end



function evalfields(source,eventdata)
    outfile = get(c_outfile,'String');
    emf = get(ct_emdir,'String');
    fmf = get(ct_fmdir,'String');
    imf = get(ct_imdir,'String');
    omf = get(ct_omdir,'String');
    
    slices.em = str2num(get(ct_emslice,'String'));
    if isempty(slices.em)
        slices.em = 1;
    end
    slices.fm = str2num(get(ct_fmslice,'String'));
    if isempty(slices.fm)
        slices.fm = 1;
    end
    slices.im = str2num(get(ct_imslice,'String'));
    if isempty(slices.im)
        slices.im = 1;
    end
    slices.om = str2num(get(ct_omslice,'String'));
    if isempty(slices.om)
        slices.om = 1;
    end
    
    init.pixelsize_lm = str2double(get(ct_px,'String'));
    init.minbeads = str2double(get(ct_bds,'String'));   
    
end


end