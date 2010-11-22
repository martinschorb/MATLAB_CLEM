function [bluespot,fluospot]=martin_chromaticshift_drift2(bm,bm2,im,im2,bmboxsize,imboxsize,fluorsel,outfileroot)

% % version MartinSchorb 100303
% %
% % usage is martin_chromaticshift_drift2(beadimage , fluorimage ,  boxsize for beads , boxsize for fluo-image , file root);
% %
% % calls digitize07 to get picks of bleed-thru fluospheres. (optional
% % input of existing coordinates)
% % Uses Gaussian fitting for getting each channels coordinates and outputs these.
% % 
% %
% % 
% % 
% % 

bmboxsize=bmboxsize+2;

gfl=0.1;

% rm=imadjust(imread(rmf));gm=imadjust(imread(gmf));%bm=imadjust(imread(bmf));
% 

[filename, pathname] = uigetfile({'shiftcoos.mat'},'select previously picked fluoshift coordinates','/struct/briggs/schorb/101108');
        if isequal(filename,0)
            disp('No previously picked positions selected');
            rgb=imadjust(bm);
            rgb(:,:,2)=imadjust(im);
            % rgb(:,:,3)=imadjust(bm);
            rgb(1,2,3)=0;
            imwrite(rgb,[outfileroot,'rgbtmp.tif'],'Compression','none');
            martin_digitize([outfileroot,'rgbtmp.tif']);
            disp('Press a key to continue...')
            pause;
            XY=evalin('base', 'XY');
            
        else          
            load([pathname,filename]);
            
        end

save([outfileroot,'_',fluorsel,'_fluoshift.shiftcoos.mat'],'XY')





% bm=rm;
%fit intensity peak to get subpixel centre

% % bmboxsize=15; % must be odd number
bmsir=(bmboxsize-1)/2;
% bmF=fspecial('gaussian',3,2);

% % imboxsize=15; % must be odd number
imsir=(imboxsize-1)/2;

% gmF=fspecial('gaussian',3,2);

mm=max([bmsir imsir]);
if ~isempty(XY)
l=length(XY);

t=find(XY<mm+1);
ii=(mod(t,l)==0)*l+mod(t,l);
XY(ii,:)=[];

s=size(im);

t=find(XY(:,2)>s(1)-mm-1);
ii=(mod(t,l)==0)*l+mod(t,l);
XY(ii,:)=[];

t=find(XY(:,1)>s(2)-mm-1);
ii=(mod(t,l)==0)*l+mod(t,l);
XY(ii,:)=[];

XY=round(XY);

bm_filtered=tom_bandpass(double(bm2),35,1344,2);
bm_filtered=double(uint16(bm_filtered));

% im_filtered=tom_bandpass(double(im2),70,1344,2);
% im_filtered=double(uint16(im_filtered));
im_filtered=im2;

 bluespot=XY;
 fluospot=XY;

    for si=1:size(XY,1)
        for iii=1:4
        sixb=bm_filtered(floor(bluespot(si,2))-bmsir:floor(bluespot(si,2))+bmsir,floor(bluespot(si,1))-bmsir:floor(bluespot(si,1))+bmsir);
        
%         sixb=double(ideal_high(sixb,.15));
%         sixb=imfilter(sixb,bmF);
        

        sixg=im_filtered(floor(fluospot(si,2))-imsir:floor(fluospot(si,2))+imsir,floor(fluospot(si,1))-imsir:floor(fluospot(si,1))+imsir);
%         sixg=double(ideal_high(sixg,.15));
%         sixg=imfilter(sixg,gmF);   
        
      

        % last argument enables interactive mode to score sub pixel fitting...
           a=cntrd1(sixb,[bmsir bmsir]+[1 1],bmboxsize-12,0);


%      a=[0 0];
%      [a(1),a(2),sx,sy,peak0D]= Gaussian2D(sixb,gfl,.75*bmboxsize);


    bluespot(si,1)=floor(bluespot(si,1))+a(1)-1-bmsir; bluespot(si,2)=floor(bluespot(si,2))+a(2)-1-bmsir;
    % %     [xpeak,ypeak,junk]=john_findpeak(sixf,1);


    b=cntrd1(sixg,[imsir imsir]+[1 1],imboxsize-12,0);

%      
%     b=[0 0];
%     [b(1),b(2),sx,sy,peak0D]= Gaussian2D(sixf,gfl,.75*imboxsize); 


    fluospot(si,1)=floor(fluospot(si,1))+b(1)-1-imsir; fluospot(si,2)=floor(fluospot(si,2))+b(2)-1-imsir;
        end


    end
else

fluospot=ones(2);
bluespot=ones(2);
end

