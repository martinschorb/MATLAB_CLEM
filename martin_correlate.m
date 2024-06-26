function martin_correlate(varargin)

% % version MartinSchorb 160321
% % Copyright EMBL 2016, All rights reserved
% %
% =========================================================================
%         DO NOT MODIFY !!!!!  use init script to set up parameters!!!  
% =========================================================================    
% %
% % usage is martin_correlate('beadimage','emimage','spotimage','outputfileroot','fluorsel','otherimage (optional)','omfluor (optional)')
% %   minimal alternative:  - no parameters
% %                         - martin_correlate(outfileroot)  
% % designed for correlating light and em images using fluorescent electron
% % dense fiducials.
% % 
% % looks for previously picked fiducial coordinates
% % 
% % calls cpselect for control point registration and uses cp2tform
% % 
% % corrects for image shift between channels using bleed-thru beads
% % 
% % uses martin_tfm_beads to suggest optimal transformation according to
% % lowest error in predictions of a single "blind" bead
% % 
% % calls martin_corr_gui4 to check and select transformation
% % output overlayed transformations and predictions to screen to select
% % best transformation for the region of interest or to modify fiducials
% % 
% % transforms fluorescence images according to the selected transformation
% % 
% % outputs files in tif format containing (transformed) images
% % outputs files in tif format representing positions of picked fiducials
% % and predictions of the selected transform
% % (output files easily overlayed in eg imagej)

if exist('corr_init','file')==2
    corr_init();
elseif exist('corr_init_orig','file')==2
    corr_init_orig(); 
else 
    a=msgbox('No initialization script found!','Error','modal');uiwait(a);
    a=msgbox('Please update algorithms!','Error','modal');uiwait(a);
    return 
end 
% global status 
status=0;

if exist('init')~=1 
    a=msgbox('Initialization script is not the newest version, please update!');uiwait(a); gaussloc = 0;
end

if nargin>1
   fmf = varargin{1};
   emf = varargin{2};
   if nargin<3
     imf=fmf;
     outfileroot=[fmf,'_corr'];
     [init,emf,fmf,imf,omf,outfile,fluorsel,omfluor,slices,hmf]=martin_correlate_init(init,outfileroot,'GFP',emf,fmf,imf);
   elseif nargin==3
       outfileroot=varargin{3};
       imf=fmf;
[init,emf,fmf,imf,omf,outfile,fluorsel,omfluor,slices,hmf]=martin_correlate_init(init,outfileroot,'GFP',emf,fmf,imf);
  else
   imf = varargin{3};
   outfileroot=varargin{4};
   fluorsel='other';
    if nargin>4
        fluorsel = varargin{5};
    end
    if nargin>5
        omf = varargin{6};
        omfluor = varargin{7};
        [init,emf,fmf,imf,omf,outfile,fluorsel,omfluor_trash,slices,hmf]=martin_correlate_init(init,outfileroot,fluorsel,emf,fmf,imf,omf,omfluor);
    else
        [init,emf,fmf,imf,omf,outfile,fluorsel,omfluor,slices,hmf]=martin_correlate_init(init,outfileroot,fluorsel,emf,fmf,imf);
    end
   end
else
    if nargin==1
        outfileroot = varargin{1};
    else
        outfileroot = '';
    end
    [init,emf,fmf,imf,omf,outfile,fluorsel,omfluor,slices,hmf]=martin_correlate_init(init,outfileroot);
end

outfileroot = outfile;
accuracy = init.accuracy/init.pixelsize_lm;



% read images and pick beads
em=martin_loadim(emf,slices.em);
fm=martin_loadim(fmf,slices.fm);im=martin_loadim(imf,slices.im);

if ~isempty(omf)
    om=martin_loadim(omf,slices.om);
else
    om=zeros(4);
end
flip = init.flip;
% adjust flip/orientation of fluorescence images
if flip==1
    fm=fm';im=im';om=om';
end
% adjust contrast of images according to init values
% em=imadjust(em);


if contr_fid==0
    fm_view=(fm);
else
    fm_view=martin_contrast(fm);
end
if contr_poi==0
    im_view=(im);
else
    im_view=martin_contrast(im);
end


em2=em;
if isa(em,'uint16')
    em=uint8(em/256);
elseif isa(em,'int16')
    em=uint8(imadjust(uint16(em))/256);
end

s_em=size(em);
s_fm=size(fm);
s_im=size(im);
em=imadjust(em);
if length(s_fm)>2
if s_fm(3)>1
    fm=fm(:,:,slices.fm);
end
end
if length(s_im)>2
if s_im(3)>1
    im=im(:,:,slices.im);
end
end   
%generate filename
file='';
%check if already previously picked
filecheck=exist([outfileroot,file,'.pickspots1.mat'],'file');

if filecheck==0 
  %import previously clicked positions
pause(0.001)
  [filename, pathname] = uigetfile('*.pickspots1.mat','select previously picked beads',loc_pickspots);
        if isequal(filename,0)
            disp('No previously picked positions selected, looking for icy coordinate files.');
            clickskip=0;
            em_ext=strfind(emf,'.');
            em_base=emf(1:em_ext(end)-1);
            emxml = [em_base,'.xml'];
            ip = xml_ec_point_import(emxml);
            
            fm_ext=strfind(fmf,'.');
            fm_base=fmf(1:fm_ext(end)-1);
            fmxml = [fm_base,'.xml'];
            bp = xml_ec_point_import(fmxml);
            if (~isempty(ip))&(~isempty(bp))
            [ip,bp]=cpselect(em,fm_view,ip,bp,'Wait',true); 
            else
            [ip,bp]=cpselect(em,fm_view,'Wait',true); 
            end   
        end
else
    in1=load([outfileroot,file,'.pickspots1.mat']);
    clickskip=1;
    ip=in1.ip;
    bp=in1.bp;
    
end


status=0;
while status==0  
%     if exist('ip2','var')>0
%         [ip,bp]=cpselect(em,fm_view,ip,bp,'Wait',true) ;
%     else
%         ip2=ip;bp2=bp;
%     end
% 190

while size(ip,1) < init.minbeads
    
    k=msgbox(['you need at least ',num2str(init.minbeads),' pairs for this transformation'],'Error','modal');
    uiwait(k);
    [ip,bp]=cpselect(em,fm_view,ip,bp,'Wait',true);
end

% computation warning
if size(ip,1) >14
    k=msgbox('Accuracy estimation might take a while when choosing too many fiducials.');
    uiwait(k);
    [ip,bp]=cpselect(em,fm_view,ip,bp,'Wait',true);
end
    
fm2=fm;
[mlen,idx]=max(s_fm);

if clickskip==0

numfids=size(ip,1);
bp1=bp;
if gaussloc > 1
    imsir=floor(imboxsize/2);
for ispot=1:numfids
    sixf=double(fm(floor(bp(ispot,2))-imsir:floor(bp(ispot,2))+imsir , floor(bp(ispot,1))-imsir:floor(bp(ispot,1))+imsir));
    [mu,sig,Amp,check] = martin_2dgaussfit(sixf,1,fit_interactive);
    
    if isnan(mu)
        bp1(ispot,:)=[NaN,NaN];
    else
        
    if check
        bp1(ispot,:)=bp(ispot,:);
        
                   
   % 178     
    else
        
        bp(ispot,:)=floor(bp(ispot,:))+mu(1:2)-[1 1]-[imsir imsir];
    end
    end

    
end
    nanidx=find(isnan(bp1(:,1)));
    bp(nanidx,:)=[];
    ip(nanidx,:)=[];    
end

% 241
%reshows the control points so you can check them...
%  ip4=ip;bp4=bp;
    [ip,bp]=cpselect(em,fm_view,ip,bp,'Wait',true) ;
%      ip2=ip4;bp2=bp4;ip=ip4;bp=bp4;

end

%export pixel values
    pos_log=[ip,bp];
    save([outfileroot,file,'.pickspots1.mat'], 'ip','bp','emf','fmf','imf','omf','hmf','slices','fluorsel','omfluor'); 

    file_1 = fopen([outfileroot,file,'_picked1.txt'],'w');
    fprintf(file_1,['Picked pixel values of corresponding fluorospheres \n\n El. Tomogram:',emf,'  slice:',num2str(slices.em),'\n Fluorospheres: ',...
        fmf,'  slice:',num2str(slices.fm),'\n fluorescence image of interest:',imf,'  slice:',num2str(slices.im),'\n other image',omf,'  slice:',num2str(slices.om),...
        '\n-----------\n EM image -  FM image\n']);
    fprintf(file_1,'%4.2f,%4.2f  -   %4.2f, %4.2f \n',pos_log'); 
    fclose(file_1);

% asks for region of interest using the control points
numfids=size(ip,1);
if exist('ipint','var')==0

% fluorsel = questdlg('What fluorescence signal are you interested in?','Signal Selector','GFP','RFP','Cancel');
        k=msgbox(['Click one spot in both images to pick region of interest     --    ',fluorsel,' Image shown on the right']);
        uiwait(k);
% switch fluorsel
%     case ''
%         return
%     case 'GFP'
%         im=gm;im_view=gm_view;imtxt='gm';
%     case 'RFP'
%         im=rm;im_view=rm_view;imtxt='rm';
% end
numspots=1;
if multispot==1;
    numq='s';
    ipint=[];
    bpint=[];    
    while ~strcmp(numq,'Correct')
        [ipint,bpint]=cpselect(em,im_view,ip,bp,'Wait',true);
        if size(ipint,1)==numfids
                k=msgbox('No spot selected!');
                uiwait(k);
                continue
        end        
        ipint=ipint(numfids+1:end,:);
        bpint=bpint(numfids+1:end,:);
        numspots=size(bpint,1);
        numq = questdlg(['You have selected ' num2str(numspots) ' fluorescent spots of interest'],'Check number of spots','Correct','No select again','Cancel');
    end

else
    ipint=[0 0];
    bpint=[0 0];
    while ~(size(ipint,1)==numfids+1&(size(bpint,1)== numfids+1))
        k=msgbox(['Click one spot in both images to pick region of interest     --    ',fluorsel,' Image shown on the right']);
        uiwait(k);
        [ipint,bpint]=cpselect(em,im_view,ip,bp,'Wait',true) ;
    end
    ipint=ipint(end,:);
    bpint=bpint(end,:);
end

numspots=size(bpint,1);
if mod(gaussloc,2) == 1
    imsir=floor(imboxsize/2);
    for ispot=1:numspots
        sixg=double(im(floor(bpint(ispot,2))-imsir:floor(bpint(ispot,2))+imsir , floor(bpint(ispot,1))-imsir:floor(bpint(ispot,1))+imsir));
        [mu,sig,Amp,check] = martin_2dgaussfit(sixg,1,fit_interactive);
        if isnan(mu)
            bpint1(ispot,:)=bpint(ispot,:);
        else

            bpint1(ispot,:)=floor(bpint(ispot,:))+mu(1:2)-[1 1]-[imsir imsir];
        end
    end
else
    bpint1=bpint;
end



% 324 reshows point of interest

[ipint,bpint]=cpselect(em,im_view,ipint,bpint1,'Wait',true) ;



end





% 336
% runs fluorescence image drift correction
if ~shift_skip
    [bluespot,fluospot]=martin_chromaticshift_drift2(fm,fm_view,im,im_view,fmboxsize,imboxsize,fluorsel,loc_shiftcoos,outfileroot,0);
    if isnan(fluospot)
        k=msgbox(['No bleed through spots found! ',fluorsel,' Image...']);
        uiwait(k);
        medshift=[]
    else  

    sdiff=fluospot-bluespot

    fspot=find(abs(sdiff)>5);
    idspot=mod(fspot,length(sdiff));
    idspot=idspot+(idspot==0)*length(sdiff);
    sdiff(idspot,:)=[];

    n_shift=size(sdiff,1);
    medshift=median(sdiff,1);
    shifterr=std(sdiff)/sqrt(n_shift);
    
    if isnan(medshift)
        medshift=[];
    end

    disp(['median of Shift correction [px]: ', num2str(medshift),' deviation: ', num2str(shifterr),' number of points: ', num2str(n_shift)]);
    disp('---  ---  ---  ---  ---  ---  ---  ---  ---');
    % disp(['Shift correction in pixel: ', num2str(medshift)]);
    %corrects for median shift of bleed-thru beads
    bpint2=bpint;
    bpint=bpint-repmat(medshift,[numspots,1]);
    show=[bpint(:,2) bpint(:,1);bpint2(:,2) bpint2(:,1)];
    end
    
    switch fluorsel
        case 'GFP'
            medshift_GFP=medshift;
        case 'RFP'
            medshift_RFP=medshift;
        otherwise
            medshift_other=medshift;
    end    
else
    medshift_GFP=[];
    medshift_RFP=[];
end

























% 407
[output,pickedem]=martin_tfm_beads(ip,bp,ipint,bpint,em,3,accuracy,init.trafo,outfileroot);
% clear test

% test(1)=sum(sum((output.all.bptfm-ip4).^2))/length(ip4);
output.emsize=s_em;
output.fmsize=s_fm;
% ip=ip4;
% bp=bp4;

tfmselect=0;
    %      generate images
         tfmed=uint8(zeros(s_em));
         newColorImage=uint8(zeros([s_em,3]));
         bptfm=output.all.bptfm;
         bpr=round(bptfm);
         m_bp=max(bpr);
        sr=size(bpr,1);
        for n=1:sr
            if min(bpr(n,:))>6 & max(bpr(n,:))<(s_em-6)
            tfmed(bpr(n,2)-5:bpr(n,2)+5,bpr(n,1)-5:bpr(n,1)+5)=255;
            tfmed(bpr(n,2),bpr(n,1))=10;end
        end   
        output.all.circle=output.all.circle(1:s_em(1),1:s_em(2));
        newColorImage(:,:,1) =uint8(output.all.circle*255)+0.8*em;
        newColorImage(:,:,2) =tfmed+0.8*em;
        newColorImage(:,:,3) =pickedem/10*255+uint8(output.all.circle*255)+0.8*em;
        output.all.image=newColorImage;

 % shows the GUI to select the transformation
status= martin_corr_gui4(output,tfmselect,status);
    
end


% ip=ip4;
% bp=bp4;
save([outfileroot,file,'.pickspots1.mat'], 'ip','bp','emf','fmf','imf','omf','hmf','bpint','slices','outfileroot'); 


% 
clear global status

file=[file,'_',fluorsel];

if tfmselect==0 % transform using all beads has been selected
    appltfm=output.all.tfm;
%     tfmed=output.all.tfmed;
    rgb=output.all.image;
    beads='all beads.';
    file=[file,'_all'];
    prederrlist=[];
    a=((output.all.bptfm-ip)).^2;
    t=[];
   
else
    appltfm=output.blind(tfmselect).optimtfm;
%     tfmed=output.blind(tfmselect).tfmed;
    rgb=output.blind(tfmselect).image;
%     pred=output.blind(tfmselect).pred;
    file=[file,'_tfm',int2str(tfmselect)];
    imwrite(pred,[outfileroot,file,'_pred.tif'],'Compression','none');
    beads=['bead ',int2str(tfmselect),' as blind bead and beads ',output.blind(tfmselect).sel(output.blind(tfmselect).optimum(1),output.blind(tfmselect).optimum(2)).beads,' as transformation base.'];
    t='\n\n Errors of the predicted beads:\n\n pixel squared - pixels - Angstrom\n--------------------------------\n';

% output a list of errors for all predicted beads to estimate accuracy
    prederrlist=output.blind(tfmselect).preddev;
    prederrlist(2,:)=sqrt(prederrlist(1,:));
    
end


% transform the fluorescence microscopy images

% [fm2 xdata ydata]=imtransform(fm,appltfm,'FillValues',128,'XData', [1 s_em(2)],'YData',[1 s_em(1)],'Size',s_em);
% [im2 xdata ydata]=imtransform(im,appltfm,'FillValues',128,'XData', [1 s_em(2)],'YData',[1 s_em(1)],'Size',s_em);
% if ~isempty(omf)
%     [om2 xdata ydata]=imtransform(om,appltfm,'FillValues',128,'XData', [1 s_em(2)],'YData',[1 s_em(1)],'Size',s_em);
% end

% [picked2 xdata ydata]=imtransform(picked,tfm,'FillValues',128,'XData', [1 s_em(2)],'YData',[1 s_em(1)],'Size',s_em);


%convert to 8 bit
% % if isa(em,'uint8')
    
% else
%     em=uint8(em2/256);
% end


%generate accuracy map

impos=tformfwd(appltfm,bpint);
impos1=round(impos);
circle1=martin_circle(em,accuracy,impos1);
impred=uint8(circle1*255)+em;


imwrite(circle1,[outfileroot,file,'_prediction.tif'],'Compression','none');
imwrite(impred,[outfileroot,file,'_pred_overlay.jpg']);

% write output files
% imwrite(fm2,[outfileroot,file,'_fm.tif'],'Compression','none');
% imwrite(em,[outfileroot,file,'_em.tif'],'Compression','none');
% imwrite(im2,[outfileroot,file,'_im.tif'],'Compression','none');

% if ~isempty(omf)
%     imwrite(om2,[outfileroot,file,'_',omfluor,'_om.tif'],'Compression','none');
% end
% imwrite(tfmed,[outfileroot,file,'_tfmed.tif'],'Compression','none');
% imwrite(pickedem,[outfileroot,file,'_pickedem.tif'],'Compression','none');
imwrite(rgb,[outfileroot,file,'_predictions.jpg']);

save([outfileroot,file,'.appltfm.mat'],'appltfm','emf','imf','omf','slices','file','circle1','fluorsel','accuracy','impos');

% save([outfileroot,file,'_tfmaccuracy.mat'],'prederrlist','allerrlist');

file_2 = fopen([outfileroot,file,'_transform.log'],'w');
fprintf(file_2,[outfileroot,file,'_transform.log      ---   Logfile of transformation\n\n']);
fprintf(file_2,['Selected transformation used: ', beads,'\n\n EM Stack: ',emf,'\n Fluorospheres: ',fmf,'  slice:',num2str(slices.fm),...
    '\n fluorescence image of interest:',imf,'  slice:',num2str(slices.im),'\n other image',omf,'  slice:',num2str(slices.om),'\n\n-----\n']);
% fprintf(file_2,['Shift error (pixel):  ',int2str(shifterr),'   #of spots used for shift: ',int2str(n_shift),'\n\n-----\n']);
% fprintf(file_2,['lowmag tomogram: ',stfile,'   Pixel size: ']);
% fprintf(file_2,'%2.3g',psize);
fprintf(file_2,'\n coordinates of transformed fluorescence spot:\n');
fprintf(file_2,'%2.2f %2.2f\n',impos');
fprintf(file_2,['\n\n prediction circle radius (px): ',int2str(accuracy),'\n\n-----------------\n']);
% fprintf(file_2,[ptext,t]);
% if isequal(psize,[])
%     fprintf(file_2,'%3.3f ,  %2.3f\n',prederrlist);
% else
%     fprintf(file_2,'%3.3f ,  %2.3f ,  %4.3f\n',prederrlist);
% end
% 
% fprintf(file_2,['\n\n Errors of all beads:\n\n pixel squared - pixels - Angstrom\n--------------------------------\n']);
% if isequal(psize,[])
%     fprintf(file_2,'%3.3f ,  %2.3f\n',allerrlist);
% else
%     fprintf(file_2,'%3.3f ,  %2.3f ,  %4.3f\n',allerrlist);
% end
fclose(file_2);

% show selected transformation prediction
% imshow(rgb)
% imshow(impred)

if init.hmauto>0
   
   if isempty(hmf)
      hmf = [emf(1:end-8),'.mrc'];
   end
       
% checks EM image file names for MRC extension


ext={'mrc','st','rec'};
lm_check=0;
hm_check=0;

    for i=1:lenght(ext)
        lm_check=or(lm_check,~isempty(findstr(upper(ext{i},emf))));
        hm_check=or(hm_check,~isempty(findstr(upper(ext{i},hmf))));
    end

   if and(lm_check,hm_check)
   
   magx=0;
   
   else
       magx=init.pixelsize_lm/init.pixelsize_hm;
   end
   hm=martin_loadim(hmf,slices.hm);
   hm=uint8(imadjust(uint16(hm))/256);
   hmslices=[slices.em,slices.hm];
  

   switch init.hmauto
       case 1
           spotpos=martin_LM2HMauto(impos,emf,hmf,init.hmcrop,magx,0,hmslices); 
       case 2
           spotpos=martin_LM2HMauto(impos,emf,hmf,init.hmcrop,magx,1,hmslices);
   end
    
   hmaccuracy=init.accuracy/init.pixelsize_hm;
   hm_accuracy=hmaccuracy;
      
   tfmcircle=martin_circle(hm,hmaccuracy,round(spotpos));

%convert to 16 bit
% tfmcircle=conv8to16bit(tfmcircle);
% % 
% if isa(lm2,'uint16')
%     lm2=lm2;
% else
%     lm2=256*uint16(em);
% end

file_2 = fopen([outfileroot,file,'_hm_transform.log'],'w');
fprintf(file_2,[outfileroot,file,'_hm_transform.log      ---   Logfile of LowMag2HighMag-transformation\n\n']);
fprintf(file_2,['lowmag file: ',outfileroot,file,'_em.tif  --  highmag file: ',hmf,'  --  slice of interest: ',slices.hm'\n\n']);
fprintf(file_2,'coordinates of transformed fluorescence spot:');
fprintf(file_2,'%2.3f %2.3f',spotpos);
fprintf(file_2,['\n prediction circle radius (px): ',int2str(hmaccuracy)]);
fclose(file_2);

save([outfileroot,file,'.sliceinfo.mat'],'hmf','slices','hm_accuracy','spotpos');

tfmcircle=uint8(tfmcircle*255);
%writes output files
% imwrite(hm,[outfileroot,file,'_hm.jpg']);
% imwrite(lm2,[outfileroot,file,'_lm2hm.tif'],'Compression','none');
% imwrite(uint8(sm/255),[outfileroot,file,'_sm.tif'],'Compression','none');
imwrite(tfmcircle,[outfileroot,file,'_hm_prediction.tif'],'Compression','none');

impred=tfmcircle+hm;
imwrite(impred,[outfileroot,file,'_hm_prd_overlay.tif']);

figure
imshow(impred);

    
end

