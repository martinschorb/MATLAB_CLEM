% function filename=martin_correlate_batch()

% % version MartinSchorb 101122
% %
% % usage is john_manualregister_beads('beadimage','emimage','gfpimage','rfpimage','outputfileroot')
% %
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
 clear all

if exist('corr_init')==2
    corr_init();
elseif exist('corr_init_orig')==2
    corr_init_orig(); 
else 
    a=msgbox('No initialization script found!','Error','modal');uiwait(a);
    a=msgbox('Please update algorithms!','Error','modal');uiwait(a);
    return 
end 





[filename, pathname] = uigetfile({'*.pickspots1.mat'},'select Fiducial coordinate files','MultiSelect', 'on','/struct/briggs/schorb/beadcoord/');

  if iscell(filename)
        indend = length(filename);
    else
        file1=filename;
        indend=1;
  end

  keep=[who;'keep'];  
  
for i=1:indend
    status=0;
    
    load([pathname,filename{i}]);
    
      file = filename{i}(1:end-15);
% 
%       dirfind = strfind(file,'/');
%       file = file(dirfind(end)+1:end);
      
    file1 = file;
    
    lmdate = martin_dbread(file,8,'/struct/briggs/schorb/_Endo-Data/outliers.csv');
    
    if isempty(lmdate)
    
    if ~isempty(str2num(file(1))) 
       if ~isempty(str2num(file(2))) 
        file = [file(3),file(1:2),file(4:end)];
       else
        file = [file(2),file(1),file(3:end)];
       end
        
    end
        lmdate = martin_dbread(file,8,'/struct/briggs/schorb/_Endo-Data/outliers.csv');

    end
    
    
    if ~isempty(str2num(file(end))) 
        filebase = file;
    else
        filebase = file(1:end-1);
    end
    
    
    
    
  


    if isempty(lmdate)
        
        disp('emergency stop');
        pause
    end
    
    if iscell(lmdate)
        lmdate = lmdate{1}(2:end-1);
    end
    
    lmdate=num2str(lmdate);
       
    while isempty(str2num(lmdate(1)))
        lmdate = lmdate(2:end);
    end
    
    if str2num(lmdate(1))>1
        lmdate = ['0' lmdate];
    end

    cd(['/struct/briggs/wanda/DataLightMicroscopy/',lmdate]);
    
    fmd = dir([filebase(1:3),'*']);
    
    if isempty(fmd)
        fmd = dir([file1(1:3),'*']);
    end
    
    if isempty(fmd)
        fmd = dir(['*',file1(1:2),'*']);
    end
    
     if isempty(fmd)
        fmd = dir(['*',file(1:2),'*']);
    end
    
    if isempty(fmd)
        fmd = dir(['*',file1(2),file1(1),'*']);
    end
    
    
    if isempty(fmd)
        disp('emergency stop');
        pause
    end
    
    if isempty(strfind(fmd(1).name,'tif'))        
        cd(fmd(1).name);
    end
    

%     
%     f=figure('visible','off');
%     imagesc(em)
%     axis off
%     colormap('gray')
%     axis equal
%     hold on
% %     scatter(ip2(nblind,1),ip2(nblind,2),4404,'rx')
%     scatter(ip(:,1),ip(:,2),442,'co')
%     set(f,'Position',[16000 1600 1600 1200])
%     
%     sp1 = strfind(file1,'_');
%     file1 = [file1(1:sp1),'{}-',file1(sp1+1:end)];
%     
%    tt = text(0,0,[file1]);%'Minimum error: ',num2str(round(testdist*5.068)),' nm   -   ',file1]);
%      
%     cd  /struct/briggs/schorb/accuracy/outliers;
%     
% %     saveas(f,['good-',num2str(oo),'.jpg']);
%     saveas(f,[file1,'.jpg']);
% end
%     oo=oo+1;
%     close all
%     
% end




if fmf(1:2)=='..'
    cd('corr')
end

if exist(gmf)==0
    
[gmf, pathname2] = uigetfile({'*.tif'},['select GFP image (',gmf,')']);

cd(pathname2);

end

outfileroot=['/struct/briggs/schorb/batchcorr/',file1,'_bcorr1'];

% read images and pick beads
% em=imread(emf);
%em=em';
fm=imread(fmf);gm=imread(gmf);rm=imread(rmf);
fm=fm';gm=gm';rm=rm';
% fm=fm.*(65535./(0.2*max(max(fm))));
    
    cd corr;
    
    emd = dir([file1,'*','_em.tif']);
    
     if isempty(emd)
       emd = dir([filebase,'*','_em.tif']);
     end
       
     
    if isempty(emd)
        disp('emergency stop');
        pause
    end
    
    
    
%     
%     
    emf = (emd(1).name);
    
    em=imread(emf);
    
% adjust contrast of images according to init values
em=imadjust(em);
if contr_b==0
    fm_view=imadjust(fm);
else
    fm_view=martin_contrast(fm);
end
if contr_g==0
    gm_view=imadjust(gm);
else
    gm_view=martin_contrast(gm);
end
if contr_r==0
    rm_view=imadjust(rm);
else
    rm_view=martin_contrast(rm);
end
if isa(em,'uint16')
    em2=em;
    em=uint8(em/256);
else
    em2=em;
end


s_em=size(em);
s_fm=size(fm);



% %check if output folder exists
% strp1=findstr(outfileroot, 'corr');
% folder=outfileroot(1:strp1-1);
% foldercheck=exist([folder,'corr']);
% 
% if foldercheck<7
%     a=msgbox('please ensure the folder ./corr is created !','Error','modal');
%     uiwait(a);
%     
% end

%generate filename
file='';
% if size(outfileroot,2)<strp1+9
%     pos=findstr(emf,'/lowmag/');
%     if pos>0
%         s=emf(1:pos-1);
%         p2=findstr(s,'/');
%         file=[file,emf(max(p2)+1:pos-1)];
%     end
% end

% 
% % get pixel size of em-tom
% pos1=strfind(emf,'/stac');
% stfile=emf(1:pos1);
% pos2=strfind(emf,'aphy/');
% stfile=[stfile,emf(pos2+5:pos2+10),'_',file,'lma.st'];
% 
% if exist(stfile)
%     [s sz]=unix(['header ',stfile]);
%     pos3=strfind(sz,'Cell axes');
%     pos4=strfind(sz(pos3+42:end),'.');
%     sz=str2num(sz(pos3+42:pos3+43+pos4(1)));
%     psize=sz/(size(em,2));ptext='';
% else
%     psize=[];
%     ptext=['NO PIXEL INFORMATION FOUND FOR TOMOGRAM '];
%     disp(ptext);
% end
% 
% %check if already previously picked
% filecheck=exist([outfileroot,file,'_pickspots1.mat']);
% filecheck2=exist([outfileroot,file,'.pickspots1.mat']);
% 
% if filecheck==0 & filecheck2==0
%   %import previously clicked positions
%         [filename, pathname] = uigetfile({'pickspots1.mat'},'select previously picked beads',loc_pickspots);
%         if isequal(filename,0)ile
%             disp('No previously picked positions selected');
%             [ip,bp]=cpselect(em,fm,'Wait',true);  
%         else          
%                a=open([pathname,filename]);
%                ip=a.ip;bp=a.bp;
%                 [ip,bp]=cpselect(em,fm,ip,bp,'Wait',true);
% status=0;
%         end
% 
% else
%     %import previously clicked positions
%     if filecheck2==2
%        load([outfileroot,file,'.pickspots1.mat']);
% %        [ip,bp]=cpselect(em,fm,ip,bp,'Wait',true);
%     else
%         load([outfileroot,file,'_pickspots1.mat']);
% %         [ip,bp]=cpselect(em,fm,ip,bp,'Wait',true);
%     end
%     
% end
% 
% 
% 
% while size(ip,1) <5
%     k=msgbox('you need at least 5 pairs for fit','Error','modal');
%     uiwait(k);
%     [ip,bp]=cpselect(em,fm,ip,bp,'Wait',true);
% end
% 
% % computation warning
% 
% if size(ip,1) >14
%     k=msgbox('Accuracy estimation might take a while when choosing too many fiducials.');
%     uiwait(k);
%     [ip,bp]=cpselect(em,fm,ip,bp,'Wait',true);
% end

    
fm2=fm;
[mlen,idx]=max(s_fm);

numfids=size(ip,1);
bp=floor(bp);ip2=ip;bp2=bp;


[ip2,bp2]=cpselect(em,fm_view,ip2,bp2,'Wait',true) ;





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
        
            
        
   % 168     
    else
        
        bp1(ispot,:)=floor(bp(ispot,:))+mu(1:2)-[1 1]-[imsir imsir];
    end
    end

    
end

    bp1(find(isnan(bp1(:,1))),:)=[];

end










































% 223
status=0;
while status==0

%reshows the control points so you can check them...
% ip4=ip;bp4=bp;
    [ip4,bp4]=cpselect(em,fm_view,ip2,bp2,'Wait',true) ;
     ip2=ip4;bp2=bp4;
     numfids=size(ip2,1);
%export pixel values
    output=[ip2,bp2];
% 234
    file_1 = fopen([outfileroot,file,'_picked1.txt'],'w');
    fprintf(file_1,['Picked pixel values of corresponding fluorospheres \n\n El. Tomogram:',emf,'\n Fluorospheres: ',fmf,'\n GFP-Image:',gmf,'\n RFP-Image',rmf,'\n-----------\n EM image -  FM image\n']);
    fprintf(file_1,'%4.0f,%4.0f -  %4.0f, %4.0f \n',output'); 
    fclose(file_1);
    
    
   save([outfileroot,file,'.pickspots1.mat'], 'ip','bp','emf','fmf','gmf','rmf'); 
    
   
   
% asks for region of interest using the control points

if exist('ipint')==0

fluorsel = questdlg('What fluorescence signal are you interested in?','Signal Selector','GFP','RFP','Cancel');

switch fluorsel
    case ''
        return
    case 'GFP'
        im=gm;im_view=gm_view;imtxt='gm';
    case 'RFP'
        im=rm;im_view=rm_view;imtxt='rm';
end
numspots=1;
if multispot==1;
    numq='s';
    ipint=[0 0];
    bpint=[0 0];    
    while ~strcmp(numq,'Correct')
        [ipint,bpint]=cpselect(em,im_view,ip4,bp4,'Wait',true);
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
    if check
        bpint1(ispot,:)=bpint(ispot,:);
    else
        
        bpint1(ispot,:)=floor(bpint(ispot,:))+mu(1:2)-[1 1]-[imsir imsir];
    end
end
end






















% 324 reshows point of interest

[ipint,bpint]=cpselect(em,im_view,ipint,bpint1,'Wait',true) ;

gm1(:,:,2)=gm(:,:,1);
gm1(:,:,3)=gm(:,:,1);
    for ispot=1:numspots
        gm1(bpint(ispot,1)-1:bpint(ispot,1)+1,bpint(ispot,2)-1:bpint(ispot,2)+1,2:3)=65000;
    end
end
% if exist(['medshift_',fluorsel]) == 0

% runs fluorescence image drift correction
if ~shift_skip
    [bluespot,fluospot]=martin_chromaticshift_drift2(fm',fm2',im',im_filtered',fmboxsize,imboxsize,fluorsel,loc_shiftcoos,outfileroot);
    if isequal(fluospot,ones(2))
        k=msgbox(['No bleed through spots found! ',fluorsel,' Image...']);
        uiwait(k);
    %     [bluespot,fluospot]=martin_chromaticshift_drift(fm',gm',gfl,fmboxsize,imboxsize,outfileroot);
    end    

    sdiff=fluospot-bluespot


    %%%%%%%%%% reversing x and y coordinates of difference because they were
    %%%%%%%%%% calculated for transposed images JB
    sdiff=circshift(sdiff,[0 1]);
    %%%%%%%%%% ---------------------------

    fspot=find(abs(sdiff)>5);
    idspot=mod(fspot,length(sdiff));
    idspot=idspot+(idspot==0)*length(sdiff);
    sdiff(idspot,:)=[];

    %n_shift=length(diff);
    %medshift=median(diff);
    %shifterr=std(diff)/sqrt(n_shift);

    %%%%%%%%% making dimension specific else fails for 1 bead JB

    n_shift=size(sdiff,1);
    medshift=median(sdiff,1);
    shifterr=std(sdiff)/sqrt(n_shift);
    if isnan(medshift) medshift=[0 0];end

    disp(['median of Shift correction [px]: ', num2str(medshift),' deviation: ', num2str(shifterr),' number of points: ', num2str(n_shift)]);
    switch fluorsel
        case 'GFP'
            medshift_GFP=medshift;
        case 'RFP'
            medshift_RFP=medshift;
    end


    % else
    % 
    % medshift=eval(genvarname(['medshift_',fluorsel]));
    % 
    % end

    disp('---  ---  ---  ---  ---  ---  ---  ---  ---');
    % disp(['Shift correction in pixel: ', num2str(medshift)]);
    %corrects for median shift of bleed-thru beads
    bpint2=bpint;
    bpint=bpint-repmat(medshift,[numspots,1]);
    show=[bpint(:,2) bpint(:,1);bpint2(:,2) bpint2(:,1)];
    
else
    medshift_GFP=[];
    medshift_RFP=[];
end

% [ipint7,bpint7]=cpselect(em,im',[ipint;ipint],show,'Wait',true) ;

%runs the transformation accuracy prediction

% [output,pickedem]=martin_ls_blind5(ip4,bp4,ipint,bpint,em,3,accuracy,outfileroot);


%%%%%%%%%%%%%%%%%%%%%%%%%

% % pasted new martin_ls_blind algorithm


% [output,pickedem]=martin_tfm_beads(ip4,bp4,ipint,bpint,em,3,accuracy,trafo,outfileroot);
% clear test

alltfm = cp2tfm(bp4,ip4,trafo);
bptfm = tformfwd(alltfm,bp4);

% test(1)=sum(sum((output.all.bptfm-ip4).^2))/length(ip4);
output.emsize=s_em;
output.fmsize=s_fm;
ip=ip4;
bp=bp4;



tfmselect = 0; % temporarily always use tfm using all beads until accuracy issue is solved!
if tfmselect>0

else  
    
    
    %      generate images
         tfmed=uint8(zeros(s_em));
         newColorImage=uint8(zeros([s_em,3]));
%          bptfm=output.all.bptfm;
         bpr=round(bptfm);
         m_bp=max(bpr);
        sr=size(bpr,1);
        
        
        %generate accuracy map

impos=tformfwd(appltfm,bpint);
impos1=round(impos);
circle1=martin_circle(em,accuracy,impos1);
impred=uint8(circle1*255)+em;
        
        
        
        for n=1:sr
            if min(bpr(n,:))>6 & max(bpr(n,:))<(s_em-6)
            tfmed(bpr(n,2)-5:bpr(n,2)+5,bpr(n,1)-5:bpr(n,1)+5)=255;
            tfmed(bpr(n,2),bpr(n,1))=10;end
        end   
%         output.all.circle=output.all.circle(1:s_em(1),1:s_em(2));
        newColorImage(:,:,1) = uint8(circle1*255)+0.8*em;
        newColorImage(:,:,2) = tfmed+0.8*em;
        newColorImage(:,:,3) = pickedem/10*255+uint8(circle1*255)+0.8*em;
        output.all.image = newColorImage;
%     else
%        output.all.image=['Bad transformation using all beads as transformation base.'];
%        disp( ['Bad transformation using all beads as transformation base.']);
%     end
    
    
end

 % shows the GUI to select the transformation
status= martin_corr_gui4(output,tfmselect,status);
    
end


ip=ip4;
bp=bp4;
save([outfileroot,file,'.pickspots1.mat'], 'ip','bp','emf','fmf','gmf','rmf',['medshift_',fluorsel],'bpint'); 


% 
clear global status

file=[file,'_',fluorsel];

if tfmselect==0 % transform using all beads has been selected
%     appltfm=output.all.tfm;
%     tfmed=output.all.tfmed;
    rgb=output.all.image;
    beads='all beads.';
    file=[file,'_all'];
%     prederrlist=[];
%     a=((output.all.bptfm-ip4)).^2;
%     t=[];
   
else
   
    
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
[gm2 xdata ydata]=imtransform(gm,appltfm,'FillValues',128,'XData', [1 s_em(2)],'YData',[1 s_em(1)],'Size',s_em);
[rm2 xdata ydata]=imtransform(rm,appltfm,'FillValues',128,'XData', [1 s_em(2)],'YData',[1 s_em(1)],'Size',s_em);
% [picked2 xdata ydata]=imtransform(picked,tfm,'FillValues',128,'XData', [1 s_em(2)],'YData',[1 s_em(1)],'Size',s_em);


%convert to 8 bit
if isa(em,'uint8')
    
else
    em=uint8(em2/256);
end






imwrite(circle1,[outfileroot,file,'_prediction.tif'],'Compression','none');
imwrite(impred,[outfileroot,file,'_pred_overlay.tif'],'Compression','none');



% write output files
imwrite(fm2,[outfileroot,file,'_fm.tif'],'Compression','none');
imwrite(em,[outfileroot,file,'_em.tif'],'Compression','none');
imwrite(gm2,[outfileroot,file,'_gm.tif'],'Compression','none');
imwrite(rm2,[outfileroot,file,'_rm.tif'],'Compression','none');
imwrite(tfmed,[outfileroot,file,'_tfmed.tif'],'Compression','none');
imwrite(pickedem,[outfileroot,file,'_pickedem.tif'],'Compression','none');
imwrite(rgb,[outfileroot,file,'_predictions.tif'],'Compression','none');

save([outfileroot,file,'.appltfm.mat'],'appltfm','emf','file','circle1','fluorsel','accuracy','impos');

% save([outfileroot,file,'_tfmaccuracy.mat'],'prederrlist','allerrlist');

file_2 = fopen([outfileroot,file,'_transform.log'],'w');
fprintf(file_2,[outfileroot,file,'_transform.log      ---   Logfile of transformation\n\n']);
fprintf(file_2,['Selected transformation used: ', beads,'\n\n EM Stack: ',emf,'\n Fluorospheres: ',fmf,'\n GFP Image: ',gmf, '\n RFP Image: ',rmf,'\n\n-----\n']);
% fprintf(file_2,['Shift error (pixel):  ',int2str(shifterr),'   #of spots used for shift: ',int2str(n_shift),'\n\n-----\n']);
fprintf(file_2,['lowmag tomogram: ',stfile,'   Pixel size: ']);
fprintf(file_2,'%2.3g',psize);
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



a=who;
excl=find(ismember(a,[keep;'a']));
a(excl)=[];
for ii=1:length(a)
clear(a{ii});
end
end

