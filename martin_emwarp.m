function [TT,hm3,scale,hm2,lm1] = martin_emwarp(inlm,inhm,crop,magx,showimage,slices)

% % version MartinSchorb 160113
% % Copyright EMBL 2016, All rights reserved
%
%usage is martin_LM2HMauto('lowmag_image','highmag_image', crop option, mag_factor, interactive)
%
% designed for correlating em images of different magnifications.
%
%
% inputs are: 
% - LM and HM images, either as Matlab arrays or as filenames for
% tif or mrc
%
% options: 
% 
% - manual cropping (1/0, default is 0, so HM is centered in LM)
% - magnification ratio between images (if not provided as mrc)
% - interactive display (1/0, default is 1 );
%
%
% outputs image transform
% outputs transformed images, crop scale (optional)

if nargin<6
    slices =[0 0];
end


if nargin<5
    showimage=1;
end


[lm,lmheader]=martin_loadim(inlm,slices(1));
[hm,hmheader]=martin_loadim(inhm,slices(2));


if isempty(lmheader) | isempty(hmheader)
   if  nargin<4 | ~isnumeric(magx) | magx==0
       a=msgbox('Please provide magnification difference or a mrc image with valid header!','Error','modal');
       error('Please provide magnification difference or a mrc image with valid header!','Error','modal');
   end
else
    
   lmpx=lmheader.MRC.xlen/lmheader.MRC.nx;
   hmpx=hmheader.MRC.xlen/hmheader.MRC.nx;
   
   magx= lmpx/hmpx;
end



if nargin<3
    crop = 0;
    % region is in the center
end




sz_l=size(lm);
sz_h=size(hm);

target_sz = sz_h/magx;

l_sz=round(target_sz);

if crop==0
    cvec = [(sz_l-l_sz)/2 l_sz-1];    
    lm1=imcrop(lm,cvec);    
else
    [lm1,cvec]=imcrop(imadjust(lm));
    close gcf
end

   
hm1=hm;

scale=1;

while scale>1/magx
    hm1=binning(hm1);
    scale=scale/2;
end

T = imregcorr(hm1, lm1);

% 
% ir_small=imref2d(size(lm1));
% hm2=imwarp(hm1,T,'OutputView',ir_small);

T1=T.T+[zeros(2,3); cvec(1:2) 0];
T1(1:2,1:2)=T1(1:2,1:2)*scale;


TT=T;
TT.T=T1;

ir=imref2d(sz_l);

hm3=imwarp(hm,TT,'OutputView',ir);


% figure, imshowpair(lm1,hm2);
if showimage == 1
    figure, imshowpair(lm,hm3);
end


