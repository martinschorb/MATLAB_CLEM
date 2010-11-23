function image = tom_bandpass1(image,low,hi,smooth)
%TOM_BANDPASS performs bandpass filtering of image or volume
%
%   image = tom_bandpass1(image,low,hi,smooth)
%   
%   An IMAGE (1, 2, or 3D) is bandpass filtered using FFT. LOW and HI
%   specifies the low- and highpass. Additionally a smoothing parameter
%   SMOOTH can be specified to avoid hard cropping (as in TOM_SPHEREMASK).
%
%PARAMETERS 
%  IN
%   image         iamge or volume to be filtered
%   low           lowest frequ (in pixels)
%   hi            highest frequ
%   smooth        smoothing (optional) - if low = 0 no smoothing around
%                   zero frequency
%
%  OUT
%   image         filtered image or volume
%
%SEE ALSO
%   TOM_FILTER
%
%    Copyright (c) 2004
%    TOM toolbox for Electron Tomography
%    Max-Planck-Institute for Biochemistry
%    Dept. Molecular Structural Biology
%    82152 Martinsried, Germany
%    http://www.biochem.mpg.de/tom
%
%   02/14/03 FF
%  last change 03/31/05 FF - docu updated
error(nargchk(3,4,nargin));
scf = 1/(size(image,1)*size(image,2)*size(image,3));
if nargin<4 % fast version
    [x,y,z]=ndgrid( -floor(size(image,1)/2):-floor(size(image,1)/2)+(size(image,1)-1),...
        -floor(size(image,2)/2):-floor(size(image,2)/2)+size(image,2)-1, ...
        -floor(size(image,3)/2):-floor(size(image,3)/2)+size(image,3)-1);
    len = sqrt(x.^2 +y.^2+z.^2);clear y; clear z;clear x;
    lowp = len <= hi;
    highp=len >= low;
    image = fftshift(tom_fourier1(image));
    image = scf*real(tom_ifourier1(ifftshift(highp.*lowp.*image)));
else % smoothened version
    image = fftshift(tom_fourier1(image));
    if low > 0
        image= real(tom_ifourier1(ifftshift(tom_spheremask1(image,hi,smooth) - tom_spheremask1(image,low,smooth))));
    else
        image= real(tom_ifourier1(ifftshift(tom_spheremask1(image,hi,smooth))));
    end;
end;
    
