function john_manualregister_LMtoHMtomo_batch()

%version MartinSchorb 100222
%
%usage is john_manualregister_LMtoHMtomo3('highmaggold_image','any highmag slice', 'outputfileroot');
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





[filename, pathname] = uigetfile({'appltfm.mat'},'select correlation transform','/struct/briggs/schorb/101118/pombe');%,'/struct/briggs/wanda/DataLightMicroscopy');
load([pathname,filename]);

% 
% gm=imread([outfileroot,file,'_gm.tif']);
% rm=imread([outfileroot,file,'_rm.tif']);

% fluorsel1 = questdlg('What fluorescence signal are you interested in?','Signal Selector','GFP','RFP','Cancel');


slice=0;

%generate filename
s=strfind(file,'_');
file1=file(1:s(end)-1);


dpod=strfind(filename,'.appl');
namebase=filename(1:dpod-1);

fluorsel='RFP';
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

   


% if exist([outfileroot,file1,'.lmhmcoos.mat'])
%     load([outfileroot,file1,'.lmhmcoos.mat'])
%     [ip,bp]=cpselect(lm,hm,ip,bp,'Wait',true);
% else
    [filename1, pathname1] = uigetfile({'lmhmcoos.mat'},'select existing picked beads','/struct/briggs/wanda/DataLightMicroscopy/100119');%,'/struct/briggs/wanda/DataLightMicroscopy');
    if ~isstr(filename1)
%         [ip,bp]=cpselect(lm,hm,'Wait',true);
    elseif exist([pathname1,filename1])==2    
        load([pathname1,filename1]);
%         [ip,bp]=cpselect(lm,hm,ip,bp,'Wait',true);
    else
%         [ip,bp]=cpselect(lm,hm,'Wait',true);
    end
% end
outfileroot=[pathname,namebase];

 lm=imread([pathname,namebase,'_em.tif']);

% read images and pick beads
%  lm=imread(lmf);
hm=imread(hmf);sm=imread(smf);


smf1=smf;
lm=imadjust(lm);


% 
% gm=conv8to16bit(gm);
% rm=conv8to16bit(rm);
hm=conv8to16bit(hm);
sm=conv8to16bit(sm);
% lm=conv8to16bit(lm);

hm=imadjust(hm);
sm=imadjust(sm);


% fm=im;


while size(ip,1) <4
    k=msgbox('you need at least 4 pairs for fit','Error','modal');
    uiwait(k);
%     [ip,bp]=cpselect(lm,hm,ip,bp,'Wait',true);
end
save([outfileroot,file1,'.lmhmcoos.mat'],'ip','bp', 'hmf','smf')
%fit beads to get subpixel centres

emboxsize=49; % must be odd number
emsir=(emboxsize-1)/2;
fmboxsize=55; % must be odd number
fmsir=(fmboxsize-1)/2;

emF=fspecial('gaussian',emsir-5,emsir/5);
fmF=fspecial('gaussian',fmsir-5,fmsir/5);
ip=round(ip);
bp=round(bp);

% for si=1:size(ip,1)
%     sixe=lm(ip(si,2)-emsir:ip(si,2)+emsir,ip(si,1)-emsir:ip(si,1)+emsir);
%     sixe=double(sixe);
%    
% %     imshow(imadjust(sixe));pause
%     sixe=(imfilter(sixe,emF)); 
%     [max_f, imax] = max(abs(sixe(:)));
%     [ypeak, xpeak] = ind2sub(size(sixe),imax(1));
%     
%     
%     
%     sixf=hm(bp(si,2)-fmsir:bp(si,2)+fmsir,bp(si,1)-fmsir:bp(si,1)+fmsir);
%     sixf=(sixf.*-1)+max(max(sixf));
% %     imshow(imadjust(sixf));pause 
%     sixf=double(imfilter(sixf,fmF));   
% 
%     
%     
%     % last argument enables interactive mode to score sub pixel fitting...
%     if (ypeak>round(emsir/2)*2-5 & ypeak<emsir+5)&(xpeak>round(emsir/2)*2-3 & xpeak<emsir+5)
%     a=cntrd1(sixe,[ypeak, xpeak],round(emsir/2)*2-5,0);
%     else
%     a=cntrd1(sixe,[emsir, emsir],round(emsir/2)*2-5,0);
%     end
% %     [xpeak,ypeak,junk]=john_findpeak(sixe,1);
%      ip2(si,1)=ip(si,1)+a(1)-1-emsir; ip2(si,2)=ip(si,2)+a(2)-1-emsir;
%     
%      
% %     [xpeak,ypeak,junk]=john_findpeak(sixf,1);
% 
%     [max_f, imax] = max(abs(sixf(:)));
%     [ypeak, xpeak] =ind2sub(size(sixe),imax(1));
%     if (ypeak>round(fmsir/2)*2-5 & ypeak<fmsir+8)&(xpeak>round(fmsir/2)*2-5 & xpeak<fmsir+8)
%     b=cntrd1(sixf,[ypeak, xpeak],round(fmsir/2)*2+7,0);
%     else
%     b=cntrd1(sixf,[fmsir, fmsir],round(fmsir/2)*2+7,0);
%     end
%      
%     bp2(si,1)=bp(si,1)+b(1)-1-fmsir; bp2(si,2)=bp(si,2)+b(2)-1-fmsir;
% end


% reshows the control points so you can check them...
% [ip3,bp3]=cpselect(lm,hm,ip2,bp2,'Wait',true);
ip3=ip;
bp3=bp;

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
spotpos=tformfwd(thm,impos);
hm_accuracy=mean([norm(tformfwd(thm,impos+sqrt(.5)*[accuracy,accuracy])-spotpos),norm(tformfwd(thm,impos+sqrt(.5)*[accuracy,-accuracy])-spotpos)]);


% [lm2 xdata ydata]=imtransform(lm,thm,'FillValues',128,'XData', [1 size(hm,2)],'YData',[1 size(hm,1)],'Size',size(hm));
% [gm2 xdata ydata]=imtransform(gm,thm,'FillValues',128,'XData', [1 size(hm,2)],'YData',[1 size(hm,1)],'Size',size(hm));
% [rm2 xdata ydata]=imtransform(rm,thm,'FillValues',128,'XData', [1 size(hm,2)],'YData',[1 size(hm,1)],'Size',size(hm));
% [picked2 xdata ydata]=imtransform(pickedlm,thm,'FillValues',128,'XData', [1 size(hm,2)],'YData',[1 size(hm,1)],'Size',size(hm));
[tfmcircle xdata ydata]=imtransform(circle1,thm,'FillValues',128,'XData', [1 size(hm,2)],'YData',[1 size(hm,1)],'Size',size(hm));

file=[file,'_',fluorsel];

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
fprintf(file_2,['\n prediction circle radius (px): ',int2str(hm_accuracy)]);
fclose(file_2);

save([outfileroot,file,'.sliceinfo.mat'],'slice','hm_accuracy','spotpos');

tfmcircle=uint16(tfmcircle/max(max(tfmcircle))*65535);
%writes output files
imwrite(hm,[outfileroot,file,'_hm.tif'],'Compression','none');
% imwrite(lm2,[outfileroot,file,'_lm2hm.tif'],'Compression','none');
imwrite(sm,[outfileroot,file,'_sm.tif'],'Compression','none');
imwrite(tfmcircle,[outfileroot,file,'_hm_prediction.tif'],'Compression','none');
% imwrite(gm2,[outfileroot,file,'_hm_gm.tif'],'Compression','none');
% imwrite(rm2,[outfileroot,file,'_hm_rm.tif'],'Compression','none');
% 
impred=tfmcircle+sm;
imwrite(impred,[outfileroot,file,'_hm_prd_overlay.tif'],'Compression','none');


imshow(impred);

% imwrite(picked2,[outfileroot,'_pickedhm.tif'],'Compression','none');
% imwrite(pickedhm,[outfileroot,'_pickedlm.tif'],'Compression','none');
% hmc = ind2rgb(hm2,jet(55525));
% figure; imshow(hmc,'XData', xdata, 'YData', ydata);
% hold on
% h = imshow(lm, gray(45536));
% set(h, 'AlphaData', 0.5);
%ylim = get(gca, 'YLim');
%set(gca, 'YLim', [0.5 ylim(2)]) 
