  function outcircle = martin_circle(em,error,fluortfm)
 
% %version MartinSchorb 100201
% %
% % martin_circle(EM-image, accuracy of transformation, coordinate of spot
% %    of interest)

%create circle around trial location

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
if fluortfm(1)-error+1>0 & fluortfm(2)-error+1>0 &fluortfm(1)+error<size(em,1) &fluortfm(2)+error<size(em,2)
    outcircle(round(fluortfm(1)-error+1):round(fluortfm(1)+error),round(fluortfm(2)-error+1):round(fluortfm(2)+error))=msk;
    outcircle=outcircle(1:size(em,1),1:size(em,2));
else
    outcircle=zeros(size(em));
end
outcircle=outcircle(:,:,1)';