function martin_imcorrelate(lmf,hmf,imf,outfileroot)

%version MartinSchorb 100222
%
%usage is john_manualregister_LMtoHMtomo3('lowmag','highmag','any highmag slice', 'outputfileroot');
%
%designed for correlating light and em images using fluorescent electron
%dense fiducials.
%
%calls cpselect for control point registration and uses cp2tform to
%calculate the transform
%
%outputs overlayer registered image to screen
%outputs files in tif format containing (transformed) images
%outputs files in tif format representing positions of picked fiducials
%(output files easily overlayed in eg imagej)


file=[];


%  if exist([outfileroot,file,'.appltfm.mat'])
%      load([outfileroot,file,'.appltfm.mat']);
%  else
%      [filename, pathname] = uigetfile({'appltfm.mat'},'select correlation transform');%,'/struct/briggs/wanda/DataLightMicroscopy');
%      load([pathname,filename]);
%  end

% 
% gm=imread([outfileroot,file,'_gm.tif']);
% rm=imread([outfileroot,file,'_rm.tif']);
% 
% fluorsel1 = questdlg('What fluorescence signal are you interested in?','Signal Selector','GFP','RFP','Cancel');
% 
% 
% dotpos1=strfind(smf,'.tif');
% getslice=[smf(1:dotpos1-1)];
% dotpos2=strfind(getslice,'.');
% slice=getslice(dotpos2(end)+1:length(getslice));
% 
% 
% %generate filename
% s=strfind(file,'_');
% file1=file(1:s(end)-1);
% 
% 
% 
% 
% 
% switch fluorsel1
%     case ''
%         return
%     case 'GFP'
% %         im=gm;
%         imtxt='gm';
%     case 'RFP'
% %         im=rm;
%         imtxt='rm';
% end

% 
% if ~isequal(fluorsel1,fluorsel)
%     k=msgbox('different fluorescence channel than in source transformation selected!','Error','modal');
%     uiwait(k);
%     fluorsel1 = questdlg('confirm fluorescence signal selection','Signal Selector','GFP','RFP','Cancel');
%     circle1=imread([outfileroot,file,'_',fluorsel,'_prediction.tif']);
% end
% 
%     
 lm=imread(lmf);
 hm=imread(hmf);
 
 hm=hm(:,:,1);
 lm=lm(:,:,1);
 
%  lm=imread([outfileroot,file,'_em.tif']);
if exist([outfileroot,file,'.lmhmcoos.mat'])
    load([outfileroot,file,'.lmhmcoos.mat'])
%     [ip,bp]=cpselect(lm,hm,ip,bp,'Wait',true);
else
    [filename1, pathname1] = uigetfile({'lmhmcoos.mat'},[file,' - select existing picked beads']);
    if ~isstr(filename1)
        [ip,bp]=cpselect(imadjust(lm),imadjust(hm),'Wait',true);
    elseif exist([pathname1,filename1])==2    
        load([pathname1,filename1]);
        [ip,bp]=cpselect(imadjust(lm),imadjust(hm),ip,bp,'Wait',true);
    else
        [ip,bp]=cpselect(imadjust(lm),imadjust(hm),'Wait',true);
    end
    
end

% % read images and pick beads

sm=imread(imf);
% 
% 
% 
% lm=imadjust(lm);
% 
% 
% % 
% % gm=conv8to16bit(gm);
% % rm=conv8to16bit(rm);
hm=conv8to16bit(hm);
% sm=conv8to16bit(sm);
% % lm=conv8to16bit(lm);
% 
% hm=imadjust(hm);
% sm=imadjust(sm);


% fm=im;
% 
% 
%     %,'/struct/briggs/wanda/DataLightMicroscopy');
%     if ~isstr(filename1)
%         [ip,bp]=cpselect(lm,hm,'Wait',true);
%     elseif exist([pathname1,filename1])==2    
%         load([pathname1,filename1]);
%         [ip,bp]=cpselect(lm,hm,ip,bp,'Wait',true);
%     else
%         [ip,bp]=cpselect(lm,hm,'Wait',true);
%     end
% end


% 
% while size(ip,1) <4
%     k=msgbox('you need at least 4 pairs for fit','Error','modal');
%     uiwait(k);
%     [ip,bp]=cpselect(lm,hm,ip,bp,'Wait',true);
% end
file1=file;

%fit beads to get subpixel centres

emboxsize=25; % must be odd number

fmboxsize=49; % must be odd number

%fit beads to get subpixel centres
ip=round(ip); bp=round(bp);

emsir=(emboxsize-1)/2;
fmsir=(fmboxsize-1)/2;
% 
% for si=1:size(ip,1)
%     sixe=lm(ip(si,2)-emsir:ip(si,2)+emsir,ip(si,1)-emsir:ip(si,1)+emsir);
%     sixe=double(sixe);
%     sixe=(sixe.*-1)+max(max(sixe));
% 
%    
%     sixf=hm(round(bp(si,2))-fmsir:bp(si,2)+fmsir,bp(si,1)-fmsir:bp(si,1)+fmsir);
%     sixf=double(sixf);
%     sixf=(sixe.*-1)+max(max(sixf));
%     
%     % last argument enables interactive mode to score sub pixel fitting...
% %      a=cntrd1(sixe,[emsir emsir],emboxsize-12,0);
% 
% a=[0 0];
% [a(1),a(2),sx,sy,peak0D]= Gaussian2D(sixe,gfl,.8*emboxsize);
% 
% 
% %     [xpeak,ypeak,junk]=john_findpeak(sixe,1);
% 
% if min(a)>0 & max(a)<emboxsize
%      ip2(si,1)=ip(si,1)+a(1)-1-emsir; ip2(si,2)=ip(si,2)+a(2)-1-emsir;
% else
%      ip2(si,:)=ip(si,:);
% end
% %     [xpeak,ypeak,junk]=john_findpeak(sixf,1);
% %      b=cntrd1(sixf,[fmsir fmsir],fmboxsize-6,0);
% 
% b=[0 0];
% [b(1),b(2),sx,sy,peak0D]= Gaussian2D(sixf,gfl,.8*fmboxsize);
% 
% if min(b)>0 & max(b)<fmboxsize
%     bp2(si,1)=bp(si,1)+b(1)-1-fmsir; bp2(si,2)=bp(si,2)+b(2)-1-fmsir;
% else
%     bp2(si,:)=bp(si,:);
%    
% end
% 
% end
% ip=ip2;
% bp=bp2;
save([outfileroot,file1,'.lmhmcoos.mat'],'ip','bp');%, 'hmf','smf')
ip3=ip;
bp3=bp;
% reshows the control points so you can check them...
% [ip3,bp3]=cpselect(lm,hm,ip2,bp2,'Wait',true);

% 
pickedhm=zeros(size(hm));
bp3r=round(bp3);
for n=1:size(bp3,1)
picked(bp3r(n,2),bp3r(n,1))=10;
end
pickedlm=zeros(size(lm));
ip3r=round(ip3);
for n=1:size(ip3,1)
pickedlm(ip3r(n,2)-1:ip3r(n,2)+1,ip3r(n,1)-1:ip3r(n,1)+1)=10;
pickedlm(ip3r(n,2),ip3r(n,1))=10;
end

thm=cp2tform(ip3,bp3,'linear conformal');
% spotpos=tformfwd(thm,impos);
% hm_accuracy=mean([norm(tformfwd(thm,impos+sqrt(.5)*[accuracy,accuracy])-spotpos),norm(tformfwd(thm,impos+sqrt(.5)*[accuracy,-accuracy])-spotpos)]);


[lm2 xdata ydata]=imtransform(sm,thm,'FillValues',128,'XData', [1 size(hm,2)],'YData',[1 size(hm,1)],'Size',size(hm));
% [gm2 xdata ydata]=imtransform(gm,thm,'FillValues',128,'XData', [1 size(hm,2)],'YData',[1 size(hm,1)],'Size',size(hm));
% [rm2 xdata ydata]=imtransform(rm,thm,'FillValues',128,'XData', [1 size(hm,2)],'YData',[1 size(hm,1)],'Size',size(hm));
% [picked2 xdata ydata]=imtransform(pickedlm,thm,'FillValues',128,'XData', [1 size(hm,2)],'YData',[1 size(hm,1)],'Size',size(hm));
% [tfmcircle xdata ydata]=imtransform(circle1,thm,'FillValues',128,'XData', [1 size(hm,2)],'YData',[1 size(hm,1)],'Size',size(hm));

% file=[file,'_',fluorsel];

%convert to 16 bit
% tfmcircle=conv8to16bit(tfmcircle);
% 
% if isa(lm2,'uint16')
%     lm2=lm2;
% else
%     lm2=256*uint16(em);
% end

% file_2 = fopen([outfileroot,file,'_hm_transform.log'],'w');
% fprintf(file_2,[outfileroot,file,'_hm_transform.log      ---   Logfile of LowMag2HighMag-transformation\n\n']);
% fprintf(file_2,['lowmag file: ',outfileroot,file,'_em.tif  --  highmag file: ',hmf,'  --  slice of interest: ',smf,'\n\n']);
% % fprintf(file_2,'coordinates of transformed fluorescence spot:');
% % fprintf(file_2,'%2.3f %2.3f',spotpos);
% % fprintf(file_2,['\n prediction circle radius (px): ',int2str(hm_accuracy)]);
% fclose(file_2);

% save([outfileroot,file,'.sliceinfo.mat'],'smf');

% tfmcircle=uint16(tfmcircle/max(max(tfmcircle))*65535);
%writes output files
% imwrite(hm,[outfileroot,file,'_hm.tif'],'Compression','none');
imwrite(lm2,[outfileroot,file,'_lm2hm.tif'],'Compression','none');
% imwrite(sm,[outfileroot,file,'_sm.tif'],'Compression','none');
% imwrite(tfmcircle,[outfileroot,file,'_hm_prediction.tif'],'Compression','none');
% imwrite(gm2,[outfileroot,file,'_hm_gm.tif'],'Compression','none');
% imwrite(rm2,[outfileroot,file,'_hm_rm.tif'],'Compression','none');
% 
% impred=tfmcircle+sm;
% imwrite(impred,[outfileroot,file,'_hm_prd_overlay.tif'],'Compression','none');

% figure
% imshow(impred);

% imwrite(picked2,[outfileroot,'_pickedhm.tif'],'Compression','none');
% imwrite(pickedhm,[outfileroot,'_pickedlm.tif'],'Compression','none');
% hmc = ind2rgb(hm2,jet(55525));
% figure; imshow(hmc,'XData', xdata, 'YData', ydata);
% hold on
% h = imshow(lm, gray(45536));
% set(h, 'AlphaData', 0.5);
%ylim = get(gca, 'YLim');
%set(gca, 'YLim', [0.5 ylim(2)]) 
