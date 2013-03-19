function martin_LMtoHM(hmf,smf,outfileroot,fit)

% % version MartinSchorb 130312
% % Copyright EMBL 2013, All rights reserved
%
%usage is martin_LMtoHM('highmaggold_image','highmag slice of interest', 'outputfileroot')
%
%designed for correlating em images using electron
%dense fiducials such as the gold beads used for tomogram reconstruction.
%
%calls cpselect for control point registration and uses cp2tform to
%calculate the transform
%
%outputs overlayer registered image to screen
%outputs files in tif format containing (transformed) images
%outputs files in tif format representing positions of picked fiducials
%(output files easily overlayed in eg imagej)

if exist('corr_init')==2
    corr_init();
elseif exist('corr_init_orig')==2
    corr_init_orig();
else 
    a=msgbox('No initialization script found!','Error','modal');uiwait(a);
    a=msgbox('Please update algorithms!','Error','modal');uiwait(a);
    return 
end 

if exist('hm_overlays')~=1 
    a=msgbox('Initialization script is not the newest version, please update!');uiwait(a); hm_overlays = 0; 
end
accuracy1 = accuracy;


if exist([outfileroot,'.appltfm.mat'],'file')
    load([outfileroot,'.appltfm.mat']);
    pathname=[];file=[];file1=[];
    namebase = outfileroot;
else

[filename, pathname] = uigetfile({'appltfm.mat'},'select correlation transform',loc_hmcoos);%,'/struct/briggs/wanda/DataLightMicroscopy');
load([pathname,filename]);

dpod=strfind(filename,'.appl');
namebase=filename(1:dpod-1);

%generate filename
s=strfind(file,'_');
file1=file(1:s(end)-1);
end
accuracy = accuracy1;
% 
% gm=imread([outfileroot,file,'_gm.tif']);
% rm=imread([outfileroot,file,'_rm.tif']);

% fluorsel1 = questdlg('What fluorescence signal are you interested in?','Signal Selector','GFP','RFP','Cancel');


slice=0;

lmf = [pathname,namebase];

if exist([lmf,'_em.tif'],'file')    
    lmf = [pathname,namebase,'_em.tif'];
elseif exist([lmf,'_em.jpg'],'file') 
    lmf = [pathname,namebase,'_em.jpg'];
else
    error('LM image file not found');
end

% read images and pick beads
lm=imread(lmf);




   
    
if ~isempty(strfind(smf,'tif'))
    sm=imread(smf);
else
    error('wrong file format, needs to be TIFF');
end
    
if strcmp(smf,hmf)
    hm = sm;
else
    if ~isempty(strfind(hmf,'tif'))
        hm=imread(hmf);
    else
        error('wrong file format, needs to be TIFF');
    end
end

if hm_overlays
    gm=imread([pathname,namebase,'_gm.tif']);
    rm=imread([pathname,namebase,'_rm.tif']);
end
smf1=smf;
lm=imadjust(lm);
 
% gm=conv8to16bit(gm);
% rm=conv8to16bit(rm);
% hm=conv8to16bit(hm);
% sm=conv8to16bit(sm);
% lm=conv8to16bit(lm);

hm=imadjust(uint16(hm));
sm=imadjust(uint16(sm));


if exist([outfileroot,file1,'.lmhmcoos.mat'],'file')
    load([outfileroot,file1,'.lmhmcoos.mat'])
    [ip,bp]=cpselect(lm,hm,ip,bp,'Wait',true);
else
    
    
    
    [filename1, pathname1] = uigetfile({'lmhmcoos.mat'},'select existing picked beads',loc_hmcoos);
    if ~isstr(filename1)
        [ip,bp]=cpselect(lm,hm,'Wait',true);
    elseif exist([pathname1,filename1])==2    
        load([pathname1,filename1]);
%         [ip,bp]=cpselect(lm,hm,ip,bp,'Wait',true);
    else
        [ip,bp]=cpselect(lm,hm,'Wait',true);
    end
end



while size(ip,1) <4
    k=msgbox('you need at least 4 pairs for fit','Error','modal');
    uiwait(k);
    [ip,bp]=cpselect(lm,hm,ip,bp,'Wait',true);
end







save([outfileroot,file1,'.lmhmcoos.mat'],'ip','bp', 'hmf','smf')

ip3=ip;bp3=bp;

thm=cp2tform(ip3,bp3,'linear conformal');
spotpos=tformfwd(thm,impos);
hmaccuracy=accuracy/pixelsize_hm;
hm_accuracy=hmaccuracy;
tfmcircle=martin_circle(sm,hmaccuracy,round(spotpos));
if exist('fluorsel','var')
file=[file,'_',fluorsel];
end
%convert to 16 bit
tfmcircle=conv8to16bit(tfmcircle);
% 
% if isa(lm2,'uint16')
%     lm2=lm2;
% else
%     lm2=256*uint16(em);
% end

file_2 = fopen([outfileroot,file,'_hm_transform.log'],'w');
fprintf(file_2,[outfileroot,file,'_hm_transform.log      ---   Logfile of LowMag2HighMag-transformation\n\n']);
fprintf(file_2,['lowmag file: ',outfileroot,file,'_em.tif  --  highmag file: ',hmf,'  --  slice of interest: ',smf1,'\n\n']);
fprintf(file_2,'coordinates of transformed fluorescence spot:');
fprintf(file_2,'%2.3f %2.3f',spotpos);
fprintf(file_2,['\n prediction circle radius (px): ',int2str(hmaccuracy)]);
fclose(file_2);

save([outfileroot,file,'.sliceinfo.mat'],'slice','hm_accuracy','spotpos');

tfmcircle=uint8(tfmcircle*255);
%writes output files
imwrite(uint8(hm/255),[outfileroot,file,'_hm.jpg']);
% imwrite(lm2,[outfileroot,file,'_lm2hm.tif'],'Compression','none');
imwrite(uint8(sm/255),[outfileroot,file,'_sm.tif'],'Compression','none');
imwrite(tfmcircle,[outfileroot,file,'_hm_prediction.tif'],'Compression','none');

if hm_overlays
 [gm2 xdata ydata]=imtransform(gm,thm,'FillValues',128,'XData', [1 size(hm,2)],'YData',[1 size(hm,1)],'Size',size(hm));
 [rm2 xdata ydata]=imtransform(rm,thm,'FillValues',128,'XData', [1 size(hm,2)],'YData',[1 size(hm,1)],'Size',size(hm));  
 imwrite(gm2,[outfileroot,file,'_hm_gm.tif'],'Compression','none');
 imwrite(rm2,[outfileroot,file,'_hm_rm.tif'],'Compression','none');
end

impred=tfmcircle+uint8(sm/255);
imwrite(impred,[outfileroot,file,'_hm_prd_overlay.jpg']);


% imshow(impred);


