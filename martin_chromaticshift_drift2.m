function [bluespot,fluospot]=martin_chromaticshift_drift2(bm,bm_view,im,im_view,bmboxsize,imboxsize,fluorsel,loc_shiftcoos,outfileroot,fit_interactive)

% % version MartinSchorb 101122
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

bmboxsize = bmboxsize + 2;
imboxsize = imboxsize + 2;

gfl=0.1;


[filename, pathname] = uigetfile({'shiftcoos.mat'},'select previously picked fluoshift coordinates',loc_shiftcoos);
 
disp(['filename: ',filename]);

if isequal(filename,0)
            disp('No previously picked positions selected');
            rgb=bm_view;
            rgb(:,:,2)=im_view;
            % rgb(:,:,3)=imadjust(bm);
            rgb(1,2,3)=0;

            martin_digitize2(rgb);
            
            disp('Press a key to continue...')
            pause;
            XY=evalin('base', 'XY');
            evalin('base','clear XY');
        else          
            load([pathname,filename]);
            
        end

save([outfileroot,'_',fluorsel,'_fluoshift.shiftcoos.mat'],'XY')


%fit intensity peak to get subpixel centre

% % bmboxsize=15; % must be odd number
bmsir=(bmboxsize-1)/2;

% % imboxsize=15; % must be odd number
imsir=(imboxsize-1)/2;

mm=max([bmsir imsir]);

if ~isempty(XY)
%     clean points too close to edge of images
l=length(XY);
s=size(im);

t=find(XY(:)<mm+1);
[row,col]=ind2sub(size(XY),t);
XY(row,:)=[];

t=find(XY(:,2)>s(1)-mm-1);
XY(t,:)=[];

t=find(XY(:,1)>s(2)-mm-1);
XY(t,:)=[];

XY=round(XY);

 bluespot=XY;
 fluospot=XY;
 
% fit 2D Gaussian

    for si=1:size(XY,1)

        sixb=bm(floor(bluespot(si,2))-bmsir:floor(bluespot(si,2))+bmsir,floor(bluespot(si,1))-bmsir:floor(bluespot(si,1))+bmsir);

        [mu,sig,Amp,check] = martin_2dgaussfit(sixb,1,fit_interactive);
        
        if ~check
             bluespot(si,:)=floor(bluespot(si,:))+mu(1:2)-[1 1]-[bmsir bmsir];
        end

            
        sixg=im(floor(fluospot(si,2))-imsir:floor(fluospot(si,2))+imsir,floor(fluospot(si,1))-imsir:floor(fluospot(si,1))+imsir);

        [mu,sig,Amp,check] = martin_2dgaussfit(sixg,1,fit_interactive);
        
        if ~check
             fluospot(si,:)=floor(fluospot(si,:))+mu(1:2)-[1 1]-[imsir imsir];
        end

    end
else

fluospot=NaN;
bluespot=NaN;
end

