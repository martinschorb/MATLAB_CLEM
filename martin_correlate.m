function martin_correlate(fmf,emf,gmf,rmf,outfileroot)

% % version MartinSchorb 101122
% %
% % usage is john_manualregister_beads('beadimage','emimage','gfpimage','rfpimage','outputfileroot');
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


if exist('corr_init')==2
    corr_init();
elseif exist('corr_init_orig')==2
    corr_init_orig() 
else 
    a=msgbox('No initialization script found!','Error','modal');uiwait(a);
    a=msgbox('Please update algorithms!','Error','modal');uiwait(a);
    return 
end 
% global status 
status=0;

% read images and pick beads
em=imread(emf);
%em=em';
fm=imread(fmf);gm=imread(gmf);rm=imread(rmf);
fm=fm';gm=gm';rm=rm';
% fm=fm.*(65535./(0.2*max(max(fm))));
em=imadjust(em);
fm=imadjust(fm);
gm=imadjust(gm);
rm=imadjust(rm);


if isa(em,'uint16')
    em2=em;
    em=uint8(em/256);
else
    em2='';
end




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


% get pixel size of em-tom
pos1=strfind(emf,'/stac');
stfile=emf(1:pos1);
pos2=strfind(emf,'aphy/');
stfile=[stfile,emf(pos2+5:pos2+10),'_',file,'lma.st'];

if exist(stfile)
    [s sz]=unix(['header ',stfile]);
    pos3=strfind(sz,'Cell axes');
    pos4=strfind(sz(pos3+42:end),'.');
    sz=str2num(sz(pos3+42:pos3+43+pos4(1)));
    psize=sz/(size(em,2));ptext='';
else
    psize=[];
    ptext=['NO PIXEL INFORMATION FOUND FOR TOMOGRAM '];
    disp(ptext);
end

%check if already previously picked
filecheck=exist([outfileroot,file,'_pickspots1.mat']);
filecheck2=exist([outfileroot,file,'.pickspots1.mat']);

if filecheck==0 & filecheck2==0
  %import previously clicked positions
        [filename, pathname] = uigetfile({'pickspots1.mat'},'select previously picked beads',loc_pickspots);
        if isequal(filename,0)
            disp('No previously picked positions selected');
            [ip,bp]=cpselect(em,fm,'Wait',true);  
        else          
               a=open([pathname,filename]);
               ip=a.ip;bp=a.bp;
                [ip,bp]=cpselect(em,fm,ip,bp,'Wait',true);
status=0;
        end

else
    %import previously clicked positions
    if filecheck2==2
       load([outfileroot,file,'.pickspots1.mat']);
%        [ip,bp]=cpselect(em,fm,ip,bp,'Wait',true);
    else
        load([outfileroot,file,'_pickspots1.mat']);
%         [ip,bp]=cpselect(em,fm,ip,bp,'Wait',true);
    end
    
end



while size(ip,1) <5
    k=msgbox('you need at least 5 pairs for fit','Error','modal');
    uiwait(k);
    [ip,bp]=cpselect(em,fm,ip,bp,'Wait',true);
end

% computation warning

if size(ip,1) >14
    k=msgbox('Accuracy estimation might take a while when choosing too many fiducials.');
    uiwait(k);
    [ip,bp]=cpselect(em,fm,ip,bp,'Wait',true);
end

    
fm2=fm;
[mlen,idx]=max(size(fm2));
fm2(:,(end+1):mlen)=fm2(:,1:(mlen-size(fm2,2)));
fm_filtered=tom_bandpass1(double(fm2),70,1344,2);
[fmean, fmax, fmin, fstd, fvariance] = tom_dev(fm_filtered);
fm_filtered=double(uint16(fm_filtered));

%fit beads to get subpixel centres
ip=floor(ip); bp=floor(bp);

% emsir=(emboxsize-1)/2;
fmsir=(fmboxsize-1)/2;
% % 
% emF=fspecial('gaussian',emsir-3,2);
%fmF=fspecial('gaussian',3,2);
sss=[];
for si=1:size(ip,1)
%     sixe=em(ip(si,2)-emsir:ip(si,2)+emsir,ip(si,1)-emsir:ip(si,1)+emsir);
%     sixe=double(sixe);
%     sixe1=(sixe.*-1)+max(max(sixe));
%     sixe2=(imfilter(sixe1,emF)); 
bp2(si,:)=bp(si,:);

% % to let it converge
for iii=1:4 
    bp(si,:)=floor(bp2(si,:));
    sixf=fm_filtered(bp(si,2)-fmsir:bp(si,2)+fmsir,bp(si,1)-fmsir:bp(si,1)+fmsir);
%     sixf1=double(ideal_high(sixf,1));
     
%     sixf2=double(imfilter(sixf,fmF));  
    
%     [C,rows]=max(sixf2);
%     [maximum,colmax]=max(C);
% 	rowmax=rows(colmax);
    
    % last argument enables interactive mode to score sub pixel fitting...
%       a=cntrd1(sixe2,[emsir+1 emsir+1],floor(.5*emboxsize),0);

% a=[0 0];
% [a(1),a(2),sx,sy,peak0D]= Gaussian2D_1(sixe,gfl,.75*emboxsize);


%     [xpeak,ypeak,junk]=john_findpeak(sixe,1);
% 
% if min(a(1:2))>0 & max(a(1:2))<emboxsize
%      ip2(si,1)=ip(si,1)+a(1)-1-emsir; ip2(si,2)=ip(si,2)+a(2)-1-emsir;
% else
%      ip2(si,:)=ip(si,:);
% end
%     [xpeak,ypeak,junk]=john_findpeak(sixf,1);
% cent1=[rowmax colmax];
% if min(cent1)<floor(.5*fmboxsize)/2 | max(cent1)>fmboxsize-floor(.5*fmboxsize)/2
%     b=cent1;
% else
    b=cntrd1(sixf,[fmsir+1 fmsir+1],floor(5),0);
% end   
% b=[0 0];
% [b(1),b(2),sx,sy,peak0D]= Gaussian2D_1(sixf,gfl,.75*fmboxsize);

if min(b(1:2))>0 & max(b(1:2))<fmboxsize
    bp2(si,1)=bp(si,1)+b(1)-1-fmsir; bp2(si,2)=bp(si,2)+b(2)-1-fmsir;
else
    bp2(si,:)=bp(si,:);
  
end
 sss(si,iii)=bp2(si,1);


end

end

% sss

status=0;
while status==0

%reshows the control points so you can check them...
% ip4=ip;bp4=bp;
    [ip4,bp4]=cpselect(em,fm,ip,bp2,'Wait',true) ;
     ip=ip4;bp=bp4;
%export pixel values
    output=[ip,bp];
        file_1 = fopen([outfileroot,file,'_picked1.txt'],'w');
%     file_2 = fopen([outfileroot,'_pickspots1.txt'],'w');
%     fprintf(file_2,'%4.0f,%4.0f,%4.0f, %4.0f \n',output);
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
        im=gm;imtxt='gm';
    case 'RFP'
        im=rm;imtxt='rm';
end


%  k=msgbox(['Click fluorescent spot to pick region of interest  - ',fluorsel,'-image shown']);
%  uiwait(k);
gm1(:,:)=gm(:,:,1);
gm1(:,:,2)=gm(:,:,1);
gm1(:,:,3)=gm(:,:,1);

% for ii=1:length(bp)
%     gm1(bp(ii,2)-1:bp(ii,2)+1,bp(ii,1)-1:bp(ii,1)+1,2:3)=65000;
% end
% imshow(gm1);
% axis([min(bp(:,1))-50 max(bp(:,1))+50 min(bp(:,2))-50 max(bp(:,2))+50]);
% [xx,yy] = ginput(1);
% close(gcf);
% bpint=floor([xx,yy]);
% ipint=[1 1];


ipint=[0 0];
bpint=[0 0];
while ~(size(ipint,1)==size(ip4,1)+1&(size(bpint,1)== size(bp4,1)+1))
% k=msgbox(['Click one spot in both images to pick region of interest     --    ',fluorsel,' Image shown on the right']);
%     uiwait(k);
    [ipint,bpint]=cpselect(em,im,ip4,bp4,'Wait',true) ;
end


ipint=ipint(end,:);
bpint=bpint(end,:);

%apply highpass-filter to eliminate cellular autofluorescence and fit intensity peak to get subpixel centre
% bb=3;
im2=im;
[mlen,idx]=max(size(im2));
im2(:,(end+1):mlen)=im2(:,1:(mlen-size(im2,2)));
im_filtered=tom_bandpass1(double(im2),70,1344,2);
im_filtered=double(uint16(im_filtered));
imsir=(imboxsize-1)/2;

for iii=1:4

    sixg=double(im_filtered(floor(bpint(2))-imsir:floor(bpint(2))+imsir,floor(bpint(1))-imsir:floor(bpint(1))+imsir));
    % sixg=ideal_high(sixg,1);
    % sixg=imfilter(sixg,fmF);
    % 
    % [C,rows]=max(sixg);
    % [maximum,colmax]=max(C);
    % rowmax=rows(colmax);

    % c=[0 0];
    % [c(1),c(2),sx,sy,peak0D]= Gaussian2D_1(sixg,gfl,.75*fmboxsize);


     c=cntrd1(sixg,[imsir imsir]+[1 1],7,0);


    if min(c(1:2))>0 & max(c(1:2))<imboxsize
        bpint=floor(bpint)+c(1:2)-[1 1]-[imsir imsir];
end

end
% reshows point of interest

[ipint,bpint]=cpselect(em,im,ipint,bpint,'Wait',true) ;


gm1(:,:,2)=gm(:,:,1);
gm1(:,:,3)=gm(:,:,1);
gm1(bpint(1)-1:bpint(1)+1,bpint(2)-1:bpint(2)+1,2:3)=65000;

end
% if exist(['medshift_',fluorsel]) == 0

% runs fluorescence image drift correction

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
bpint=bpint-medshift;

show=[bpint(2) bpint(1);bpint2(2) bpint2(1)];


% [ipint7,bpint7]=cpselect(em,im',[ipint;ipint],show,'Wait',true) ;

%runs the transformation accuracy prediction

% [output,pickedem]=martin_ls_blind5(ip4,bp4,ipint,bpint,em,3,accuracy,outfileroot);


%%%%%%%%%%%%%%%%%%%%%%%%%

% % pasted new martin_ls_blind algorithm


[output,pickedem]=martin_tfm_beads(ip4,bp4,ipint,bpint,em,3,accuracy,outfileroot);

clear test
test(1)=sum(sum((output.all.bptfm-ip4).^2))/length(ip4);
ip=ip4;
bp=bp4;

%     
%     calculates the accuracy for the predicted region:
for i=1:size(output.blind,2)
    rowmin=output.blind(i).rowmin;
    colmin=output.blind(i).colmin;
    test(i+1)=output.blind(i).minimum;
%     output.blind(i).predacc=2*((0.08*output.blind(i).sel(rowmin,colmin).ls)/median(d.ls)+15*d.optcloserr/median(d.optcloserr)+13*d.optclosdist/median(d.optclosdist)).^0.6-1+d.optprederr;
    
end    

[mmum select]=min(test);

% the index of selected transformation
tfmselect=select-1;

if tfmselect>0
    nblind=select;
    %generate images
    tfmed=uint8(zeros(size(em)));
    bpot2=tformfwd(output.blind(nblind).optimtfm,output.blind(nblind).sel(rowmin,colmin).bp);
    bpr=round(bpot2);
    bpotfm=output.blind(nblind).bpotfm;
    pred=uint8(zeros(size(em)));
    ppr=round([bpotfm;output.blind(nblind).sel(rowmin,colmin).blindtfm]);

    if (max(max([bpr;ppr]))<2042) && (min(min([bpr;ppr]))>6)

        for n=1:size(bpr,1)
            tfmed(bpr(n,2)-5:bpr(n,2)+5,bpr(n,1)-5:bpr(n,1)+5)=10;
            tfmed(bpr(n,2),bpr(n,1))=10;
        end



        for n=1:size(ppr,1)
            pred(ppr(n,2)-5:ppr(n,2)+5,ppr(n,1)-5:ppr(n,1)+5)=10;
            pred(ppr(n,2),ppr(n,1))=10;
        end


        newColorImage(:,:,1) =pred/10*255+uint8(output.blind(nblind).circle*255)+em;
        newColorImage(:,:,2) =tfmed/10*255+em;
        newColorImage(:,:,3) =pickedem/10*255+uint8(output.blind(nblind).circle*255)+em;
        output.blind(nblind).image=newColorImage;
    else
        output.blind(nblind).image=['blind: ',int2str(nblind),' -- ',output.blind(nblind).sel(rowmin,colmin).beads];
        disp( ['Bad transformation using bead ',int2str(nblind),' as blind bead']);disp(['and beads ',output.blind(nblind).sel(rowmin,colmin).beads,' as transformation base.']);
    end
    
    
else
    %      generate images
         tfmed=uint8(zeros(size(em)));
         newColorImage=uint8(zeros([size(em),3]));
         bptfm=output.all.bptfm;
         bpr=round(bptfm);
    if (max(max(bpr))<size(em,1)-6) && (min(min(bpr))>6)
        for n=1:size(bpr,1)
            tfmed(bpr(n,2)-5:bpr(n,2)+5,bpr(n,1)-5:bpr(n,1)+5)=10;
            tfmed(bpr(n,2),bpr(n,1))=10;
        end   
        output.all.circle=output.all.circle(1:size(em,1),1:size(em,2));
        newColorImage(:,:,1) =uint8(output.all.circle*255)+em;
        newColorImage(:,:,2) =tfmed/10*255+em;
        newColorImage(:,:,3) =pickedem/10*255+uint8(output.all.circle*255)+em;
        output.all.image=newColorImage;
    else
       output.all.image=['Bad transformation using all beads as transformation base.'];
       disp( ['Bad transformation using all beads as transformation base.']);
    end
    
    
    
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
% if exist([outfileroot,file,'_transforms.mat'])==0
%     save([outfileroot,file,'_transforms.mat'], 'output');
% end



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

[fm2 xdata ydata]=imtransform(fm,appltfm,'FillValues',128,'XData', [1 size(em,2)],'YData',[1 size(em,1)],'Size',size(em));
[gm2 xdata ydata]=imtransform(gm,appltfm,'FillValues',128,'XData', [1 size(em,2)],'YData',[1 size(em,1)],'Size',size(em));
[rm2 xdata ydata]=imtransform(rm,appltfm,'FillValues',128,'XData', [1 size(em,2)],'YData',[1 size(em,1)],'Size',size(em));
% [picked2 xdata ydata]=imtransform(picked,tfm,'FillValues',128,'XData', [1 size(em,2)],'YData',[1 size(em,1)],'Size',size(em));


%convert to 16 bit
if isa(em2,'uint16')
    em=em2;
else
    em=256*uint16(em);
end


%generate accuracy map

impos=tformfwd(appltfm,bpint);
impos1=round(impos);
circle1=martin_circle(em,accuracy,impos1);
impred=uint16(circle1*655535)+em;



imwrite(circle1,[outfileroot,file,'_prediction.tif'],'Compression','none');
imwrite(impred,[outfileroot,file,'_pred_overlay.tif'],'Compression','none');



% write output files
imwrite(fm2,[outfileroot,file,'_fm.tif'],'Compression','none');
imwrite(em,[outfileroot,file,'_em.tif'],'Compression','none');
imwrite(gm2,[outfileroot,file,'_gm.tif'],'Compression','none');
imwrite(rm2,[outfileroot,file,'_rm.tif'],'Compression','none');
imwrite(tfmed,[outfileroot,file,'_tfmed.tif'],'Compression','none');
imwrite(pickedem,[outfileroot,file,'_pickedem.tif'],'Compression','none');
% imwrite(rgb,[outfileroot,file,'_predictions.tif'],'Compression','none');

save([outfileroot,file,'.appltfm.mat'],'appltfm','emf','file','circle1','fluorsel','accuracy','impos');

% save([outfileroot,file,'_tfmaccuracy.mat'],'prederrlist','allerrlist');

file_2 = fopen([outfileroot,file,'_transform.log'],'w');
fprintf(file_2,[outfileroot,file,'_transform.log      ---   Logfile of transformation\n\n']);
fprintf(file_2,['Selected transformation used: ', beads,'\n\n EM Stack: ',emf,'\n Fluorospheres: ',fmf,'\n GFP Image: ',gmf, '\n RFP Image: ',rmf,'\n\n-----\n']);
% fprintf(file_2,['Shift error (pixel):  ',int2str(shifterr),'   #of spots used for shift: ',int2str(n_shift),'\n\n-----\n']);
fprintf(file_2,['lowmag tomogram: ',stfile,'   Pixel size: ']);
fprintf(file_2,'%2.3g',psize);
fprintf(file_2,'\n coordinates of transformed fluorescence spot:');
fprintf(file_2,'%2.3f %2.3f',impos);
fprintf(file_2,['\n prediction circle radius (px): ',int2str(accuracy),'\n\n-----------------\n']);
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

