function c = tom_ifourier(a)

%TOM_IFOURIER Calculates the inverse Fourier Transformation of Matrix A
%   C=TOM_FOURIER(A) This function calculates the inverse Fourier Transformation of the
%   input matrix A.
%   Matrix A can be a 1-D, 2-D or 3-D matrix. the result is stored in the matrix C.
%
%   Example ... to come
%  ---------
%
%
%   see also TOM_FOURIER, TOM_FILTER
%
%   09/11/02   AL
%
%   Copyright (c) 2004
%    TOM toolbox for Electron Tomography
%    Max-Planck-Institute for Biochemistry
%    Dept. Molecular Structural Biology
%    82152 Martinsried, Germany
%    http://www.biochem.mpg.de/tom

[s1 s2 s3]=size(a);

if s1 == 1 | s2 == 1
    c=ifft(a);
elseif s3 == 1
    c=ifft2(a);
else
    c=ifftn(a);
end

