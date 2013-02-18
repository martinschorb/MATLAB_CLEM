function martin_correlate(varargin)

% % version MartinSchorb 130107
% % 
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
   imf = varargin{3};
   outfileroot=varargin{4};
   fluorsel='other';
    if nargin>4
        fluorsel = varargin{5};
    end
    if nargin>5
        omf = varargin{6};
        omfluor = varargin{7};
        [init,emf,fmf,imf,omf,outfile,fluorsel,omfluor,slices]=martin_correlate_init(init,outfileroot,fluorsel,emf,fmf,imf,omf,omfluor);
    else
        [init,emf,fmf,imf,omf,outfile,fluorsel,omfluor,slices]=martin_correlate_init(init,outfileroot,fluorsel,emf,fmf,imf);
    end
else
    if nargin==1
        outfileroot = varargin{1};
    else
        outfileroot = '';
    end
    [init,emf,fmf,imf,omf,outfile,fluorsel,omfluor,slices]=martin_correlate_init(init,outfileroot);
end

outfileroot = outfile;

accuracy = init.accuracy/init.pixelsize_lm;
% read images and pick beads
em=imread(emf);
fm=imread(fmf,slices.fm);im=imread(imf,slices.im);
if ~isempty(omf)
    om=imread(omf,slices.om);
else
    om=zeros(4);
end

% adjust flip/orientation of fluorescence images
if flip==1
    fm=fm';im=im';om=om';
end
% adjust contrast of images according to init values
em=imadjust(em);
if contr_fid==0
    fm_view=imadjust(fm);
else
    fm_view=martin_contrast(fm);
end
if contr_poi==0
    im_view=imadjust(im);
else
    im_view=martin_contrast(im);
end

if isa(em,'uint16')
    em2=em;
    em=uint8(em/256);
else
    em2=em;
end
s_em=size(em);
s_fm=size(fm);

%generate filename
file='';

% get pixel size of em-tom
% pos1=strfind(emf,'/stac');
% stfile=emf(1:pos1);
% pos2=strfind(emf,'aphy/');
% stfile=[stfile,emf(pos2+5:pos2+10),'_',file,'lma.st'];
% 
% if exist(stfile,'file')
%     [s sz]=unix(['header ',stfile]);
%     pos3=strfind(sz,'Cell axes');
%     pos4=strfind(sz(pos3+42:end),'.');
%     sz=str2num(sz(pos3+42:pos3+43+pos4(1)));
%     psize=sz/(s_em(2));ptext='';
% else
%     psize=[];
%     ptext=['NO PIXEL INFORMATION FOUND FOR TOMOGRAM '];
%     disp(ptext);
% end

%check if already previously picked
filecheck=exist([outfileroot,file,'_pickspots1.mat'],'file');
filecheck2=exist([outfileroot,file,'.pickspots1.mat'],'file');

if filecheck==0 & filecheck2==0
  %import previously clicked positions
pause(0.001)
  [filename, pathname] = uigetfile('*.pickspots1.mat','select previously picked beads',loc_pickspots);
        if isequal(filename,0)
            disp('No previously picked positions selected');            
            [fm_view1,rotid] = martin_rotateimage(em,fm_view);
            [ip,bp]=cpselect(em,fm_view1,'Wait',true); 
            bp=martin_coordinate_sort(bp,rotid,s_fm);
        else          
               a=open([pathname,filename]);
               ip=a.ip;bp=a.bp;
                [ip,bp]=cpselect(em,fm_view,ip,bp,'Wait',true);
status=0;
        end
else
    load([outfileroot,file,'.pickspots1.mat']);
end

status=0;
while status==0
    
    if exist('ip2','var')>0
        [ip,bp]=cpselect(em,fm_view,ip2,bp2,'Wait',true) ;
    else
        ip2=ip;bp2=bp;
    end
% 135
while size(ip2,1) < init.minbeads
    k=msgbox(['you need at least ',num2str(init.minbeads),' pairs for this transformation'],'Error','modal');
    uiwait(k);
    [ip2,bp2]=cpselect(em,fm_view,ip2,bp2,'Wait',true);
end

% computation warning
if size(ip,1) >14
    k=msgbox('Accuracy estimation might take a while when choosing too many fiducials.');
    uiwait(k);
    [ip2,bp2]=cpselect(em,fm_view,ip2,bp2,'Wait',true);
end
    
fm2=fm;
[mlen,idx]=max(s_fm);

numfids=size(ip2,1);
bp1=bp2;
if gaussloc > 1
    imsir=floor(imboxsize/2);
for ispot=1:numfids
    sixf=double(fm(floor(bp2(ispot,2))-imsir:floor(bp2(ispot,2))+imsir , floor(bp2(ispot,1))-imsir:floor(bp2(ispot,1))+imsir));
    [mu,sig,Amp,check] = martin_2dgaussfit(sixf,1,fit_interactive);
    if isnan(mu)
        bp1(ispot,:)=[NaN,NaN];
    else
    if check
        bp1(ispot,:)=bp2(ispot,:);
        
            
        
   % 168     
    else
        
        bp2(ispot,:)=floor(bp2(ispot,:))+mu(1:2)-[1 1]-[imsir imsir];
    end
    end

    
end
    nanidx=find(isnan(bp1(:,1)));
    bp2(nanidx,:)=[];
    ip2(nanidx,:)=[];    
end













































% 226
%reshows the control points so you can check them...
%  ip4=ip;bp4=bp;
    [ip4,bp4]=cpselect(em,fm_view,ip2,bp2,'Wait',true) ;
     ip2=ip4;bp2=bp4;ip=ip4;bp=bp4;
     numfids=size(ip2,1);
%export pixel values
    output=[ip4,bp4];
% 234

    save([outfileroot,file,'.pickspots1.mat'], 'ip','bp','emf','fmf','imf','omf','slices','fluorsel','omfluor'); 

    file_1 = fopen([outfileroot,file,'_picked1.txt'],'w');
    fprintf(file_1,['Picked pixel values of corresponding fluorospheres \n\n El. Tomogram:',emf,'  slice:',num2str(slices.em),'\n Fluorospheres: ',...
        fmf,'  slice:',num2str(slices.fm),'\n fluorescence image of interest:',imf,'  slice:',num2str(slices.im),'\n other image',omf,'  slice:',num2str(slices.om),...
        '\n-----------\n EM image -  FM image\n']);
    fprintf(file_1,'%4.2f,%4.2f  -   %4.2f, %4.2f \n',output'); 
    fclose(file_1);
    
    
    
   
   
% asks for region of interest using the control points

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
        [ipint,bpint]=cpselect(em,im_view,ip4,bp4,'Wait',true);
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
        [ipint,bpint]=cpselect(em,im_view,ip4,bp4,'Wait',true) ;
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
        if check | isnan(mu)
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
    [bluespot,fluospot]=martin_chromaticshift_drift2(fm,fm_view,im,im_view,fmboxsize,imboxsize,fluorsel,loc_shiftcoos,outfileroot,fit_interactive);
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
    end
    
else
    medshift_GFP=[];
    medshift_RFP=[];
end


























% 407
[output,pickedem]=martin_tfm_beads(ip4,bp4,ipint,bpint,em,3,accuracy,init.trafo,outfileroot);
% clear test

% test(1)=sum(sum((output.all.bptfm-ip4).^2))/length(ip4);
output.emsize=s_em;
output.fmsize=s_fm;
ip=ip4;
bp=bp4;

%     calculates the accuracy for the predicted region:
% for i=1:size(output.blind,2)
%     rowmin=output.blind(i).rowmin;
%     colmin=output.blind(i).colmin;
%     test(i+1)=output.blind(i).minimum;
% %     output.blind(i).predacc=2*((0.08*output.blind(i).sel(rowmin,colmin).ls)/median(d.ls)+15*d.optcloserr/median(d.optcloserr)+13*d.optclosdist/median(d.optclosdist)).^0.6-1+d.optprederr;
%     
% end    

% [mmum select]=min(test);

% the index of selected transformation
% tfmselect=select-1;
tfmselect = 0; % temporarily always use tfm using all beads until accuracy issue is solved!
if tfmselect>0
%     nblind=select;
%     %generate images
%     tfmed=uint8(zeros(s_em));
%     bpot2=tformfwd(output.blind(nblind).optimtfm,output.blind(nblind).sel(rowmin,colmin).bp);
%     bpr=round(bpot2);
%     bpotfm=output.blind(nblind).bpotfm;
%     pred=uint8(zeros(size_em));
%     ppr=round([bpotfm;output.blind(nblind).sel(rowmin,colmin).blindtfm]);
%     m_bp=max(bpr);m_all=max([bpr;ppr]);
%     if (m_all(1)<s_em(1)-6) && (m_all(2)<s_em(2)) && (min(min([bpr;ppr]))>6)
% 
%         for n=1:size(bpr,1)
%             tfmed(bpr(n,2)-5:bpr(n,2)+5,bpr(n,1)-5:bpr(n,1)+5)=10;
%             tfmed(bpr(n,2),bpr(n,1))=10;
%         end
% 
% 
% 
%         for n=1:size(ppr,1)
%             pred(ppr(n,2)-5:ppr(n,2)+5,ppr(n,1)-5:ppr(n,1)+5)=10;
%             pred(ppr(n,2),ppr(n,1))=10;
%         end
% 
% 
%         newColorImage(:,:,1) =pred/10*255+uint8(output.blind(nblind).circle*255)+em;
%         newColorImage(:,:,2) =tfmed/10*255+em;
%         newColorImage(:,:,3) =pickedem/10*255+uint8(output.blind(nblind).circle*255)+em;
%         output.blind(nblind).image=newColorImage;
%     else
%         output.blind(nblind).image=['blind: ',int2str(nblind),' -- ',output.blind(nblind).sel(rowmin,colmin).beads];
%         disp( ['Bad transformation using bead ',int2str(nblind),' as blind bead']);disp(['and beads ',output.blind(nblind).sel(rowmin,colmin).beads,' as transformation base.']);
%     end
%     
    
else
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
%     else
%        output.all.image=['Bad transformation using all beads as transformation base.'];
%        disp( ['Bad transformation using all beads as transformation base.']);
%     end
    
    
end
%%%%%%%%%%%%%%%%%%% end of pasted ls_blind


% clear test
% test(1)=sum(sum((output.all.bptfm-ip4).^2))/length(ip4);
% ip=ip4;
% bp=bp4;
% 
% %     
% %     calculates the accuracy for the predicted region:
% for i=1:size(output.blind,2)
%     rowmin=output.blind(i).rowmin;
%     colmin=output.blind(i).colmin;
%     test(i+1)=output.blind(i).minimum;
% %     output.blind(i).predacc=2*((0.08*output.blind(i).sel(rowmin,colmin).ls)/median(d.ls)+15*d.optcloserr/median(d.optcloserr)+13*d.optclosdist/median(d.optclosdist)).^0.6-1+d.optprederr;
%     
% end    
% 
% [mmum select]=min(test);
% 
% % the index of selected transformation
% tfmselect=select-1;
% 
% 
%  % shows the GUI to select the transformation
% status= martin_corr_gui3(output,tfmselect,status);
%     
% end
% 
% 
% ip=ip4;
% bp=bp4;
% save([outfileroot,file,'.pickspots1.mat'], 'ip','bp','emf','fmf','gmf','rmf',['medshift_',fluorsel],'bpint'); 
% 
% 
% % 
% clear global status

% 
% if exist([outfileroot,file,'_transforms.mat'],'file')==0
%     save([outfileroot,file,'_transforms.mat'], 'output');
% end



 % shows the GUI to select the transformation
status= martin_corr_gui4(output,tfmselect,status);
    
end


ip=ip4;
bp=bp4;
save([outfileroot,file,'.pickspots1.mat'], 'ip','bp','emf','fmf','imf','omf',['medshift_',fluorsel],'bpint','slices'); 


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
    a=((output.all.bptfm-ip4)).^2;
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
%     if ~isempty(psize)
%         prederrlist(3,:)=prederrlist(2,:)*psize;
%     end
    % output a list of errors for all beads to estimate accuracy
%     a=(output.blind(tfmselect).sel(output.blind(tfmselect).optimum(1),output.blind(tfmselect).optimum(2)).bptfm(:,1:2)-output.blind(tfmselect).ip).^2;

   
    
end

% allerrlist(1,:)=a(:,1)+a(:,2);
% allerrlist(2,:)=sqrt(allerrlist(1,:));
% 
% if ~isempty(psize) 
%     
%     allerrlist(3,:)=allerrlist(2,:)*psize;
% end





% transform the fluorescence microscopy images

[fm2 xdata ydata]=imtransform(fm,appltfm,'FillValues',128,'XData', [1 s_em(2)],'YData',[1 s_em(1)],'Size',s_em);
[gm2 xdata ydata]=imtransform(im,appltfm,'FillValues',128,'XData', [1 s_em(2)],'YData',[1 s_em(1)],'Size',s_em);
if ~isempty(omf)
    [om2 xdata ydata]=imtransform(om,appltfm,'FillValues',128,'XData', [1 s_em(2)],'YData',[1 s_em(1)],'Size',s_em);
end

% [picked2 xdata ydata]=imtransform(picked,tfm,'FillValues',128,'XData', [1 s_em(2)],'YData',[1 s_em(1)],'Size',s_em);


%convert to 8 bit
if isa(em,'uint8')
    
else
    em=uint8(em2/256);
end


%generate accuracy map

impos=tformfwd(appltfm,bpint);
impos1=round(impos);
circle1=martin_circle(em,accuracy,impos1);
impred=uint8(circle1*255)+em;



imwrite(circle1,[outfileroot,file,'_prediction.tif'],'Compression','none');
imwrite(impred,[outfileroot,file,'_pred_overlay.jpg']);



% write output files
imwrite(fm2,[outfileroot,file,'_fm.tif'],'Compression','none');
imwrite(em,[outfileroot,file,'_em.tif'],'Compression','none');
imwrite(gm2,[outfileroot,file,'_gm.tif'],'Compression','none');

if ~isempty(omf)
    imwrite(om2,[outfileroot,file,'_',omfluor,'.tif'],'Compression','none');
end
imwrite(tfmed,[outfileroot,file,'_tfmed.tif'],'Compression','none');
imwrite(pickedem,[outfileroot,file,'_pickedem.tif'],'Compression','none');
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
fprintf(file_2,'%2.2f %2.2f\n',impos);
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
imshow(impred)

