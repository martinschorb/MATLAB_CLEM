function im=conv8to16bit(im)

%converts an image from 8bit to 16 bit
% Version 20130313 Martin Schorb

s=size(im,3);

im1=uint16(im);
im2=uint16(zeros(size(im)));

if isa(im,'uint8')
    for i=1:s
        im2(:,:,i)=256*im1(:,:,i);
    end
    im=im2;
else
    disp('ERROR: wrong image format')
end
