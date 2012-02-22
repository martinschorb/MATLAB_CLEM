  function outcircle = martin_circle(em,error,fluortfm)
 
% %version MartinSchorb 120222
% %
% % martin_circle(EM-image, accuracy of transformation, coordinate of spot
% %    of interest)
% 
% creates a circle around a given location; if multiple coordinates are
% given, the spots are numbered

r=round(error);
m = 2*r;
x = 0:(m-1) ;
cent = (m-1)/2;
x2 = (x-cent).^2;
dst=zeros(m,m);
for i=1:m
    dst(i,:)=sqrt((i-cent)^2+x2);
end

ind=find(dst < r & dst > r-7);

msk=zeros([2*r,2*r]);
msk(ind)=1.0;

outcircle=zeros(size(em));
numspots=size(fluortfm,1);
ping=0;
for jj=1:numspots

if fluortfm(jj,1)-error+1>0 & fluortfm(jj,2)-error+1>0 &fluortfm(jj,1)+error<size(em,2) &fluortfm(jj,2)+error<size(em,1)
    outcircle(round(fluortfm(jj,2)-error+1):round(fluortfm(jj,2)+error),round(fluortfm(jj,1)-error+1):round(fluortfm(jj,1)+error))=msk;
%     outcircle=outcircle(1:size(em,1),1:size(em,2));
    ping=1;
end

if numspots>1
    hf = figure('color','white','units','pixel','Visible', 'off');
    text('units','pixels','position',[100 100],'fontsize',30,'string',num2str(jj));
    set(gca,'units','pixels','position',[1 1 200 200],'visible','off');
    tim = getframe(gca);
    close all;
    tim2 = tim.cdata;
    tim3 = imcrop(tim2,[90 80 56 32]);
    tim4 = ~tim3(:,:,1);
    outcircle(round(fluortfm(jj,2)+error+12-16):round(fluortfm(jj,2)+error+12+16),round(fluortfm(jj,1))-28:round(fluortfm(jj,1))+28)=tim4;
end

end

if ~ping
    outcircle=zeros(size(em));
end



end
% outcircle=outcircle(:,:,1)';