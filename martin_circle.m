  function outcircle = martin_circle(em,error,fluortfm)
 
% % version MartinSchorb 160627
% % Copyright EMBL 2011-2016, 
% %
%      This file is part of EMBL script for high-accuracy CLEM.
%  
%      EMBL script for high-accuracy CLEM is free software: you can redistribute it and/or modify
%      it under the terms of the GNU General Public License as published by
%      the Free Software Foundation, either version 3 of the License, or
%      (at your option) any later version.
%  
%      EMBL script for high-accuracy CLEM is distributed in the hope that it will be useful,
%      but WITHOUT ANY WARRANTY; without even the implied warranty of
%      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%      GNU General Public License for more details.
%  
%      You should have received a copy of the GNU General Public License
%      along with EMBL script for high-accuracy CLEM.  If not, see <http://www.gnu.org/licenses/>.
%  

% %
% % martin_circle(EM-image, accuracy of transformation, coordinate of spot
% %    of interest)
% 
% creates a circle around a given location; if multiple coordinates are
% given, the spots are numbered

r=round(error);
fluortfm = round(fluortfm);
m = 2*r;
x = 0:(m-1) ;
cent = (m-1)/2;
x2 = (x-cent).^2;
dst=zeros(m,m);
for i=1:m
    dst(i,:)=sqrt((i-cent)^2+x2);
end

ind=find(dst < r & dst > r-7);

msk=logical(zeros([2*r,2*r]));
msk(ind)=1;
s_o=size(em);
outcircle=logical(zeros(s_o));
numspots=size(fluortfm,1);
ping=0;
for jj=1:numspots

if fluortfm(jj,1)-error+1>0 & fluortfm(jj,2)-error+1>0 &fluortfm(jj,1)+error<s_o(2) &fluortfm(jj,2)+error<s_o(1)
    ar1=outcircle(round(fluortfm(jj,2)-error+1):round(fluortfm(jj,2)+error),round(fluortfm(jj,1)-error+1):round(fluortfm(jj,1)+error));
    outcircle(round(fluortfm(jj,2)-error+1):round(fluortfm(jj,2)+error),round(fluortfm(jj,1)-error+1):round(fluortfm(jj,1)+error))=bitor(ar1,msk);
%     outcircle=outcircle(1:size(em,1),1:size(em,2));
    ping=1;
end

if numspots>1
    hf = figure('color','white','units','pixel','Visible', 'off');
    text('units','pixels','position',[100 100],'HorizontalAlignment','center','fontsize',80,'string',num2str(jj));
    set(gca,'units','pixels','position',[1 1 200 200],'visible','off');
    orig_mode = get(hf, 'PaperPositionMode');
    set(hf, 'PaperPositionMode', 'auto');
    cdata = hardcopy(hf, '-Dzbuffer', '-r0');
    
    
%     tim = getframe(gca,[80 85 57 41]);
    close all;
%     tim2 = tim.cdata;
    tim3=~cdata(:,:,1);
    tim4 = imcrop(tim3,[79 310 42 30]);

 
if fluortfm(jj,1)>21 & fluortfm(jj,2)+error+13-15>0 &fluortfm(jj,1)+21<size(em,2) & fluortfm(jj,2)+error+13+15<size(em,1)
    area=outcircle(round(fluortfm(jj,2)+error+13-15):round(fluortfm(jj,2)+error+13+15),round(fluortfm(jj,1))-21:round(fluortfm(jj,1))+21);
    outcircle(round(fluortfm(jj,2)+error+13-15):round(fluortfm(jj,2)+error+13+15),round(fluortfm(jj,1))-21:round(fluortfm(jj,1))+21)=bitor(area,tim4);
end

end

if ~ping
    outcircle=zeros(size(em));
end



end
% outcircle=outcircle(:,:,1)';