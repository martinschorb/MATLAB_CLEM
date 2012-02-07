 function [stats] = martin_beads_analysis2(pos,allpos,index)

%version MartinSchorb 091229
%
%usage is martin_beads_analysis(bead coordinates used for tfm, all bead coord, blind bead index);
%
%designed for calculating statistics on fiducial distribution.
%
% calculates the estimated position of every bead according to the
% transformation obtained using the other fiducials
%
% output is a vector with number of beads,  standard deviation, mean
% ellipticity , mean coord. of cloud.
%
%  stats=[numbeads , SD from centre, ellipticity(ratio of the major axis lengths), mean_x,mean_y, dist of blind bead from major axis, sd of rotated in minor axis direction]

backup=allpos;
% pos(index,:)=[];

m1=mean(pos(:,1));
m2=mean(pos(:,2));
sd=std([pos(:,1);pos(:,2)]);
pos2=pos;
pos2(:,1)=pos(:,1)-m1;
pos2(:,2)=pos(:,2)-m2;


alpha=atan2(m1,m2);
rot=[cos(alpha) sin(alpha);-sin(alpha) cos(alpha)];
rpos=pos;
rpos=rot*pos2';
rpos=rpos';
sdmx1=std(rpos(:,1));
sdmx2=std(rpos(:,2));


backup2(:,1)=backup(:,1)-m1;
backup2(:,2)=backup(:,2)-m2;
rall=rot*backup2';
stats(6)=rall(1,index)^2;
stats(7)=sdmx1;
%output;

numbeads=size(pos,1);
ellipt=max(sdmx1,sdmx2)/min(sdmx1,sdmx2);


stats(1)=numbeads;
stats(2)=sd;
stats(3)=ellipt;
stats(4:5)=[m1 m2];