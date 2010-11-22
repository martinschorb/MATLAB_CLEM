function c = tom_fourier(a)
%TOM_FOURIER Calculates the Fourier Transformation of Matrix A
%   C=TOM_FOURIER(A) This function calculates the Fourier Transformation of the input matrix A.
%   Matrix A can be a 1-D, 2-D or 3-D matrix. the result is stored in the matrix C.
%
%   Example
%  ---------
%
%          11    3   -6  -14             16      3+19i   -6      3+19i
%       a = 1  -13    9    6        c = -24-2i   50+9i    8-64i  1+21i
%           0    8   13   -3             8       5-37i    54     5+37i
%          -8   17  -15    7            -24+2i   1-21i    8+64i  59-9i
%
%
%   see also TOM_FILTER, TOM_ORCD, TOM_RTSYM
%
%   08/07/02   AL
%
%
[s1 s2 s3]=size(a);

if s1 == 1 | s2 == 1
    c=fft(a);
elseif s3 == 1
    c=fft2(a);
else
    c=fftn(a);
end

