function [im,head]=martin_loadim(varargin)

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
        if head.Size(3)>1
            if v==0
                kk=1:head.Size(3);
                [s,v]=listdlg('PromptString', 'Select slice','SelectionMode', 'single', 'ListString', num2str(kk') , 'Name' , imf);
            end
            im=rot90((imf.Value(:,:,v)));
        else    
            im=rot90(imf.Value);
        end
    elseif ~isempty(strfind(in,'tif'))
        
        head=imfinfo(in);
        
        if head.PhotometricInterpretation=='RGB'
            im=imread(in);
        elseif v==0;
            im=imread(in);
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