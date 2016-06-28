 function [output,pickedem] = martin_tfm_beads(ip,bp,pos,bpint,em,kmin,accuracy,trafo,outfileroot)

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
% %usage is martin_tfm_beads(bead coordinates,position of interest,EM image,minimum number of beads used for transformations,acc,transformation type, output folder);
% %
% % 
% %
% %designed for estimating the accuracy of correlating light and em images using fluorescent electron
% %dense fiducials. 
% % Modified to always use the transform given by all beads (the choice of
% transform is sample dependent, this does perform reliably in any case.


%create timestamp

d=clock;
d1='';
for i=1:size(d,2)
    if d(i)<10
        d1=[d1,'0',int2str(d(i))];
    else
        d1=[d1,int2str(d(i))];
    end
end
stamp=d1(3:12);

%saves initial coordinates
ip2=ip;
bp2=bp;
ntot=size(ip2,1);
pos=mean(pos,1);
%generate picked image
sz_em=size(em);
pickedem=uint8(zeros(sz_em));
ip2r=round(ip2);
for n=1:size(ip2,1)
pickedem(max(ip2r(n,2)-5,1):min(ip2r(n,2)+5,sz_em(1)),max(ip2r(n,1)-5,1):min(ip2r(n,1)+5,sz_em(2)))=10;
pickedem(ip2r(n,2),ip2r(n,1))=10;
end

% output.pickedem=pickedem;

devall=zeros(1,ntot);
output=struct;
nblind=1;
output.blind=struct;
output.blind.sel=struct;

% get bead closest to site of interest


   dist1=[0 0];        
    for cl=1:ntot
        dist1(cl)=norm(pos-ip2(cl,:));
    end

     tfm=cp2tform(bp2,ip2,trafo);
     output.all.tfm=tfm;           
     output.all.circle=martin_circle(em,accuracy,tformfwd(tfm,bpint));
     %transform coordinates 
     bptfm=tformfwd(tfm,bp2);
     output.all.bptfm=bptfm;


end 
