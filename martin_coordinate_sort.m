function out=martin_coordinate_sort(coos,rotid,imsize)

% % version MartinSchorb 130312
% % Copyright EMBL 2013, All rights reserved
% 
% manages coordinates after #of 90degree rotations (rotid)

rotid=rem(rotid,4);
l=size(coos,1);

if ~rem(rotid,2)
    imsize=fliplr(imsize);
end

switch rotid
    case 0
        out=coos;
    case 1
        out(:,1)=repmat(imsize(2),[l,1])-coos(:,2)+1;
        out(:,2)=coos(:,1);
    case 2
        out(:,2)=repmat(imsize(2),[l,1])-coos(:,2)+1;
        out(:,1)=repmat(imsize(1),[l,1])-coos(:,1)+1;
    case 3
        out(:,1)=coos(:,2);
        out(:,2)=repmat(imsize(1),[l,1])-coos(:,1)+1;
end
      
        