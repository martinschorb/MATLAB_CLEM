function [mu,sig,A,check] = martin_2dgaussfit(im,linear,interactive,c0)

% 2D-Gauss Fit (i.e. of one peak to a PSF-convolved image)
% 
% usage is [mu,sig,check] = martin_2dgaussfit(image,linear_background_enable,interactive_mode_enable,initial_guess_of_center)
% 
% version MartinSchorb 130107
% 
% with the inputs: image (has to be square) , indicators whether linear
% background subtraction and interactive mode should be enabled, and
% optional initial coordinates of the center.
%
% returns: mu - the fitted centre of the gaussian, sig - its standard dev.,
%   A - its Amplitude, check - an indicator whether the outputs are
%   feasible.



if nargin<3
    interactive = 0;
end

if nargin<2
    linear = 0;
end

check = 0;


im=double(im);
sz=size(im);
im1=im;

if sz(1)~=sz(2)
    error('expect square image');
end



% determine background from edge pixels
l=im(1,1:end-1);r=im(end,2:end);u=im(2:end,1);d=im(1:end-1,end);

ln=sz(1)-1;
back = [l(:);r(:);u(:);d(:)];
backcoos = [ones(1,ln),(ln+1)*ones(1,ln),2:ln+1,1:ln ; 1:ln,2:ln+1,ones(1,ln),(ln+1)*ones(1,ln)]';

back_s = std(back);
backg1=median(back);

% check for and remove outlier peaks at the edge
outpeaks = find(back > (backg1 + 2.5*back_s));

[row,col] = ind2sub([length(back)/4 4],outpeaks);

const = ones(4*ln,1);

const(outpeaks)=[];
back(outpeaks)=[];
backcoos(outpeaks,:)=[];

[X,Y]=meshgrid(1:sz(1),1:sz(2));

% substract linear gradient background if wanted
switch linear
    case 0
        linstr = 'off';
        IM1 = mean(back(:));
    case 1
        linstr = 'on';
        coeff = [fliplr(backcoos),const]\back;
        IM1 = coeff(1)*X + coeff(2)*Y + coeff(3);
end

im1 = im - IM1;

% im = im - mean(back)

for i=1:length(row)
    switch num2str(col(i))
        case '1'
            im1(1,row(i)) = 0;
        case '2'
            im1(end,row(i)+1) = 0;
        case '3'
            im1(row(i)+1,1) = 0;
        case '4'
            im1(row(i),end) = 0;
    end
end

%  initial parameters


    
    centerim = im1(2:end-1,2:end-1);
    maxpeak = max(centerim(:));
    
if nargin<4
    
    middle=ceil(sz(1:2)/2);

else

    middle = c0;
end



%Define a function which returns the residual between image and your fitted gaussian
gauss2D = @(params) params(1)/params(4)*exp(-(1/(2*params(4))*((X(:)-params(2)).^2+(Y(:)-params(3)).^2))) - im1(:);

%Define initial guesses for parameters
params0 = [maxpeak,middle,mean(sz)/4];

%not displaying debugging info
opts = optimset('Display','Off');

%Fit
fitparams = lsqnonlin(gauss2D,params0,[],[],opts) ;

A = fitparams(1);
mu = fitparams(2:3);
sig = fitparams(4);





% check feasability of parameters

if abs(log(A)-log(maxpeak)) > 1;
    check = 1;
end

if sig > sz(1)^2
    check = 1;
end

if sig < 0
    check =1;
end

if max(mu) > sz(1)-1
    check = 1;
end

if min(mu) < 2
    check = 1;   
end
check = 1;
    
if (interactive + check) > 1
    [X1,Y1]=meshgrid(1:0.2:sz,1:0.2:sz);
    Z=A/sig*exp(-(1/(2*sig)*((X1-mu(1)).^2+(Y1-mu(2)).^2)));
    
    hFig = figure('MenuBar','none','Toolbar','figure','NumberTitle','off','Visible','off','Position',[0,120,900,900],'Name','Check fit');
%     clean toolbar
    tbh = findall(hFig,'Type','uipushtool');
    delete(tbh);
    tbh = findall(hFig,'Type','uitoggletool');
    delete(tbh([1:2,end]));
    rotate3d on;
%     


%     buttons
    h_keepbutton = uicontrol('Parent',hFig,'Style','PushButton','Units','Normalized','Position',[0.15 0.02 .1 0.03],'Callback',@keepbutton,'String','keep fit');
    h_oldbutton = uicontrol('Parent',hFig,'Style','PushButton','Units','Normalized','Position',[0.3 0.02 .1 0.03],'Callback',@oldbutton,'String','keep original','ForegroundColor','m');
    h_manualbutton = uicontrol('Parent',hFig,'Style','PushButton','Units','Normalized','Position',[0.45 0.02 .1 0.03],'Callback',@manualbutton,'String','select manually');
    h_redobutton = uicontrol('Parent',hFig,'Style','PushButton','Units','Normalized','Position',[0.6 0.02 .1 0.03],'Callback',@redobutton,'String','redo fit');
    h_kickbutton = uicontrol('Parent',hFig,'Style','PushButton','Units','Normalized','Position',[0.75 0.02 .1 0.03],'Callback',@kickbutton,'String','kick out point');
    
    
    
    display1 = uicontrol('Style','text','String',['centre of fit: ',num2str(mu),'   linear background subtraction: ',linstr],'Units','Normalized','Position',[0.15 0.97 .6 .024]);
    
%     graphics
    im2 = im-mean(im(:));
    surf(im2);
    hold all
    aa=mesh(X1,Y1,Z);
    set(aa,'LineWidth',2)
    set(aa,'EdgeColor','k');
    set(aa,'FaceAlpha',0);
    bb = plot3(repmat(middle(1),[2 1]),repmat(middle(2),[2 1]),[0 max([im1(:);Z(:)])]);
    set(bb,'LineWidth',3);
    set(bb,'Color','m');
    
    set(hFig,'Visible','on');
    
    uiwait
    uiresume
    close all


end


% uiwait(hFig)
% -----------------------------------

function keepbutton(h_keepbutton,event)
    close all
    check = 0;
end

% -----------------------------------

function oldbutton(h_slider2,event)
    close all
    mu = middle;
    A=NaN; sig=NaN;
    check = 0;
end


% -----------------------------------

function manualbutton(h_manualbutton,event)
    close all
    [xxx,mu] = cpselect(uint16(im1),imadjust(uint16(im1)),mu,mu,'Wait',true);
    
%     reduce box size
    cs = round(log2(sz));
    newsize = max(sz - 2*cs,[3 3]);
    dis1 = ceil(newsize/2)-1;
    corner = min(max(round(mu)-dis1,1),sz-newsize+1);
    im1 = imcrop(im,[corner,newsize-1]);
    c1 = max(mu-corner+dis1-1,[1 1]);
    [mu1,sig,A,check] = martin_2dgaussfit(im1,linear,2,c1);
    if c1 == [1 1]
        mu = mu1;
    else
        mu = mu1-dis1+corner;
    end
end


% -----------------------------------

function redobutton(h_redobutton,event)
    close all
    [mu,sig,A,check] = martin_2dgaussfit(im,~(linear),2);
    
    
end


% -----------------------------------

function kickbutton(h_redobutton,event)
    close all
    mu = NaN;
    A=NaN; sig=NaN;
    check = 0;
    
end

end