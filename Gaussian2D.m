function [cx,cy,sx,sy,PeakOD] = Gaussian2D(m,tol,radius);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% To fit a 2-D gaussian
%% m = Image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [cx,cy,sx,sy,PeakOD] = Gaussian2D(m,tol);
% 
% [sizey sizex] = size(m);
% [x,y] = MeshGrid(1:sizex,1:sizey);
% fit = abs(PeakOD)*(exp(-0.5*(x-cx).^2./(sx^2)-0.5*(y-cy).^2./(sy^2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% a function to fit a thermal cloud 2-D

%% m = image
%% tol = fitting tolerance

options = optimset('Display','off','TolFun',tol,'LargeScale','off');
m=double(m);
[sizey sizex] = size(m);
[cx,cy,sx,sy] = centerofmass(m);

sy=radius;sx=radius;

pOD = max(max(m));

mx = double(m(ceil(cy),:));
x1D = double(1:sizex);
ip1D = double([cx,sx,pOD]);
fp1D = fminunc(@fitgaussian1D,ip1D,options,mx,x1D);

if fp1D(1)>0 & fp1D(1)<sizex
    cx = fp1D(1);
end

sx = fp1D(2);
PeakOD = fp1D(3);

my = double(m(:,ceil(cx))');
y1D = double(1:sizey);
ip1D = double([cy,sy,pOD]);
fp1D = fminunc(@fitgaussian1D,ip1D,options,my,y1D);

if fp1D(1)>0 & fp1D(1)<sizey
    cy = fp1D(1);
end

sy = fp1D(2);
PeakOD = fp1D(3);
[X,Y] = meshgrid(1:sizex,1:sizey);

initpar = [cx,cy,sx,sy,PeakOD];
fp = fminunc(@fitgaussian2D,initpar,options,m,X,Y);
cx = fp(1);
cy = fp(2);
sx = fp(3);
sy = fp(4);
PeakOD = fp(5);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: find c of m of distribution
function [cx,cy,sx,sy] = centerofmass(m);

[sizey sizex] = size(m);
vx = sum(m);
vy = sum(m');

vx = vx.*(vx>0);
vy = vy.*(vy>0);

x = [1:sizex];
y = [1:sizey];

cx = sum(vx.*x)/sum(vx);
cy = sum(vy.*y)/sum(vy);

sx = sqrt(sum(vx.*(abs(x-cx).^2))/sum(vx));
sy = sqrt(sum(vy.*(abs(y-cy).^2))/sum(vy));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z] = fitgaussian1D(p,v,x);

%cx = p(1);
%wx = p(2);
%amp = p(3);

zx = p(3)*exp(-0.5*(x-p(1)).^2./(p(2)^2)) - v;

z = sum(zx.^2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z] = fitgaussian2D(p,m,X,Y);

%cx = p(1);
%cy = p(2);
%wx = p(3);
%wy = p(4);
%amp = p(5);

ztmp = p(5)*(exp(-0.5*(X-p(1)).^2./(p(3)^2)-0.5*(Y-p(2)).^2./(p(4)^2))) - m;

z = sum(sum(ztmp.^2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


