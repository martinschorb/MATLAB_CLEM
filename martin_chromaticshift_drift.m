function [bluespot,fluospot]=martin_chromaticshift_drift(bm,im,gfl,bmboxsize,imboxsize,outfileroot)

% % version MartinSchorb 100818
% %
% % usage is martin_chromaticshift_drift(beadimage , fluorimage , gauss-fit parameter , boxsize for beads , boxsize for fluo-image , file root);
% %
% % calls digitize07 to get picks of bleed-thru fluospheres. (optional
% % input of existing coordinates)
% % Uses Gaussian fitting for getting each channels coordinates and outputs these.
% % 
% %
% % 
% % 
% % 



gfl=0.1;

% rm=imadjust(imread(rmf));gm=imadjust(imread(gmf));%bm=imadjust(imread(bmf));
% 

[filename, pathname] = uigetfile({'shiftcoos.mat'},'select previously picked fluoshift coordinates','/struct/briggs/wanda/DataLightMicroscopy');
        if isequal(filename,0)
            disp('No previously picked positions selected');
            rgb=imadjust(bm);
            rgb(:,:,2)=imadjust(im);
            % rgb(:,:,3)=imadjust(bm);
            rgb(1,2,3)=0;
%             imwrite(rgb,[outfileroot,'rgbtmp.tif'],'Compression','none');
%             martin_digitize([outfileroot,'rgbtmp.tif']);
            martin_digitize2(rgb);
            disp('Press a key to continue...')
            pause;
            XY=evalin('base', 'XY');            
        else          
            load([pathname,filename]);
            
        end

save([outfileroot,'_fluoshift.shiftcoos.mat'],'XY')





% bm=rm;
%fit intensity peak to get subpixel centre

bmboxsize=bmboxsize+2; % must be odd number
bmsir=(bmboxsize-1)/2;
% bmF=fspecial('gaussian',bmboxsize-3,2);

% imboxsize=15; % must be odd number
imsir=(imboxsize-1)/2;

% gmF=fspecial('gaussian',gmboxsize-3,2);

% rmboxsize=15; % must be odd number
% rmsir=(rmboxsize-1)/2;
% rmF=fspecial('gaussian',rmboxsize-3,2);


mm=max([bmsir imsir]);
if ~isempty(XY)&(XY~=[0 0])
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

    for si=1:size(XY,1)
        sixb=bm(XY(si,2)-bmsir:XY(si,2)+bmsir,XY(si,1)-bmsir:XY(si,1)+bmsir);
    %     sixb=double(imfilter(sixb,bmF));
       sixb=double(sixb);

        sixf=im(XY(si,2)-imsir:XY(si,2)+imsir,XY(si,1)-imsir:XY(si,1)+imsir);
    %     sixg=double(imfilter(sixg,gmF));   
        sixf=double(sixf);

    %     sixr=rm(XY1(si,2)-rmsir:XY1(si,2)+rmsir,XY1(si,1)-rmsir:XY1(si,1)+rmsir);
    % %     sixr=double(imfilter(sixr,rmF));
    %    sixr=double(sixr);


        % last argument enables interactive mode to score sub pixel fitting...
    %       a=cntrd1(sixb,[bmsir bmsir],bmboxsize-12,0);


     a=[0 0];
     [a(1),a(2),sx,sy,peak0D]= Gaussian2D(sixb,gfl,.75*bmboxsize);


    bluespot(si,1)=XY(si,1)+a(1)-1-bmsir; bluespot(si,2)=XY(si,2)+a(2)-1-bmsir;
    % %     [xpeak,ypeak,junk]=john_findpeak(sixf,1);


    %      b=cntrd1(sixg,[gmsir gmsir],gmboxsize-6,0);

    b=[0 0];
    [b(1),b(2),sx,sy,peak0D]= Gaussian2D(sixf,gfl,.75*imboxsize); 


    fluospot(si,1)=XY(si,1)+b(1)-1-imsir; fluospot(si,2)=XY(si,2)+b(2)-1-imsir;


    %      c=cntrd1(sixr,[rmsir rmsir],rmboxsize-6,0);

    %      
    % c=[0 0];
    % [c(1),c(2),sx,sy,peak0D]= Gaussian2D(sixr,gfl,.75*rmboxsize);


    %  redspot(si,1)=XY(si,1)+c(1)-1-rmsir; redspot(si,2)=XY(si,2)+c(2)-1-rmsir;



    end
else

fluospot=ones(2);
bluespot=ones(2);
end


% 
% rgb=imadjust(rm);
% rgb(:,:,2)=imadjust(gm);
% rgb(:,:,3)=imadjust(bm);


% allspots=[bluespot;greenspot;redspot];
% coordinates=[bluespot,fluospot];
% reshows point of interest

% cpselect(imadjust(gm),rgb,allspots,allspots);