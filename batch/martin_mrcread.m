function a = martin_mrcread(varargin)

% Version Martin Schorb 101128
% 
%  based on readMRCfile (MATLABcentral)
% 
% usage is outvol = martin_mrcread('filename',[xmin xmax ymin ymax zmin zmax])
% 
%  second parameter describes region of the subvolume to be extracted from
%  MRC-file.




if nargin>1
    fname=varargin{1};
    res = varargin{2};
elseif nargin < 2
    if nargin==1
        if isstr(varargin{1})
            fname=varargin{1};
            res=[];
         elseif isnumeric(varargin{1})
             res=varargin{1};
        end
    else 
        res=[];
    end

    [filename, pathname] = uigetfile({'*.rec';'*.st'},'select MRC file');
    if filename==0
        error('no file specified');
        a= -1;
        return;
    end
    fname=[pathname,filename];
end    


%% read MRCfile 
[fid,message]=fopen(fname,'r');

if fid == -1
    error('can''t open file');
    a= -1;
    return;
end

%%

nx=fread(fid,1,'long');
ny=fread(fid,1,'long');
nz=fread(fid,1,'long');
type= fread(fid,1,'long');
% fprintf(1,'nx= %d ny= %d nz= %d type= %d', nx, ny,nz,type);


%%

% Seek to start of data
status=fseek(fid,1024,-1);

%%


if isempty(res)

    % Shorts
    if type== 1
        a=fread(fid,nx*ny*nz,'int16');

    %floats
    elseif type == 2
        a=fread(fid,nx*ny*nz,'float32');


    else
        a=fread(fid,nx*ny*nz,'uint8');
    end

    fclose( fid);
      a= reshape(a, [nx ny nz]);
    if nz == 1
        a= reshape(a, [nx ny]);
    end

else
    % Shorts
    
    xmin=res(1)-1;
    xmax=res(2);
    ymin=res(3)-1;
    ymax=res(4);
       
    x=xmin+1:xmax;
    y=ymin+1:ymax;

    lx=length(x);
    ly=length(y);
    
    
    if length(res)>5;
        zmin=res(5)-1;
        zmax=res(6);
        z=zmin+1:zmax;
    else
        zmin=0;
        zmax=1;
        z=1;
    end
    
    lz=length(z);
    
 

    
   if type == 1
        ftype='int16';
        mult=2;
        a = int16(zeros (lx,ly,lz));
    %floats
    elseif type == 2
       ftype='float32';
       mult=4;
       a = double(zeros (lx,ly,lz));
    else
      ftype='uint8';
      mult=1;
      a = uint8(zeros (lx,ly,lz));
    end 
    
    
%%    

% discard values outside the specified region

    junk=fseek(fid,(zmin*nx*ny)*mult,'cof');
    
    for iz=1:lz
        junk=fseek(fid,(ymin*nx)*mult,'cof');
        for iy=1:ly
            junk=fseek(fid,(xmin)*mult,'cof');
            a(1:lx,iy,iz)=(fread(fid,lx,ftype))';
            junk=fseek(fid,(nx-xmax)*mult,'cof');
        end
        junk=fseek(fid,(ny-ymax)*nx*mult,'cof'); 
    end   
        fclose( fid);         
end             

   



