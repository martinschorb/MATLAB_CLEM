function [outfile,in_dir,pxs,trafo,minbeads]=martin_corr_accuracy_init(loc_pickspots,pixelsize_lm,trafo)

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
% 

% init values

in_dir = loc_pickspots;
pxs = pixelsize_lm;
outfile = 'corr_accuracy';
minbeads = 3;


f = figure('Visible','off','Position',[0,120,1000,1000],'NumberTitle','off','Menubar','none','Name','Correlation Accuracy Determination - Initial Configuration');


t_outfile = uipanel('Title','output file name','FontSize',12,'BackgroundColor',[0.85 0.85 0.85],...
          'Position',[0.05,0.8,0.4,0.08]);
c_outfile = uicontrol('Parent',t_outfile,'Style','edit','BackgroundColor',[0.95 0.95 0.95],...
           'Position',[10,10,380,40],'String',outfile);  
%    ---

t_dir = uipanel('Title','starting input directory','FontSize',12,'BackgroundColor',[0.85 0.85 0.85],...
          'Position',[0.05,0.7,0.4,0.08]);
ct_dir = uicontrol('Parent',t_dir,'Style','edit','BackgroundColor',[0.95 0.95 0.95],...
           'Position',[10,10,300,40],'String',in_dir);       
but_dir = uicontrol('Parent',t_dir,'Style','pushbutton','BackgroundColor',[0.95 0.95 0.95],...
           'Position',[315,10,65,40],'String','Browse','Callback',{@but_dir_Callback});      
%   ---      
      
t_px = uipanel('Title','EM pixel size [nm]','FontSize',12,'BackgroundColor',[0.85 0.85 0.85],...
          'Position',[0.05,0.6,0.4,0.08]);
ct_px = uicontrol('Parent',t_px,'Style','edit','BackgroundColor',[0.95 0.95 0.95],...
           'Position',[10,10,300,40],'String',pxs); 

%   ---
    
g_trafo = uibuttongroup('Title','Transformation type to use','FontSize',12,'BackgroundColor',[0.85 0.85 0.85],...
          'Position',[0.05,0.45,0.4,0.13]);
u0 = uicontrol('Style','Radio','String','linear conformal (default)','BackgroundColor',[0.85 0.85 0.85],...
    'pos',[10 70 300 30],'parent',g_trafo,'HandleVisibility','off');
u1 = uicontrol('Style','Radio','String','affine','BackgroundColor',[0.85 0.85 0.85],...
    'pos',[10 40 300 30],'parent',g_trafo,'Callback',{@aff_Callback});
u2 = uicontrol('Style','Radio','String','projective','BackgroundColor',[0.85 0.85 0.85],...
    'pos',[10 10 300 30],'parent',g_trafo,'Callback',{@proj_Callback});
% Initialize button group properties. 
set(g_trafo,'SelectionChangeFcn',@selcbk);      
      
%   ---
t_bds = uipanel('Title','Minimum number of beads','FontSize',12,'BackgroundColor',[0.85 0.85 0.85],...
          'Position',[0.05,0.35,0.4,0.08]);
ct_bds = uicontrol('Parent',t_bds,'Style','edit','BackgroundColor',[0.95 0.95 0.95],...
           'Position',[160,10,40,40],'String',minbeads); 



%   ---
go = uicontrol('Style','pushbutton','BackgroundColor',[0.95 0.95 0.95],...
           'Position',[50 50 30 50],'String','Go','Callback',{@go_Callback}); 


% movegui(f,'center')      
set(f,'Visible','on');

uiwait(f);

% -------------------------- Callback functions  ----------------------

function but_dir_Callback(source,eventdata)   
    in_dir = uigetdir('select previously picked beads',in_dir);
    set(ct_dir,'String',in_dir)
end



function  selcbk(source,eventdata)
    trafo = get(eventdata.NewValue,'String');
    switch trafo
        case 'linear conformal (default)'
            trafo = 'linear conformal';
            minbeads = 3;
            set(ct_bds,'String','3');
        case 'affine'
            minbeads = 3;
            set(ct_bds,'String','3');
        case 'projective'
            minbeads = 4;
            set(ct_bds,'String','4')
    end
end

function  aff_Callback(source,eventdata)
    minbeads = 3;
    set(ct_bds,'String','3')
end


function  proj_Callback(source,eventdata)
    minbeads = 4;
    set(ct_bds,'String','4')
end

function go_Callback(source,eventdata)
    outfile = get(c_outfile,'String');
    in_dir = get(ct_dir,'String');
    pxs = str2num(get(ct_px,'String'));
    minbeads = str2num(get(ct_bds,'String'));
    close(f)
end

end