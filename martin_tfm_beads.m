 function [output,pickedem] = martin_tfm_beads(ip,bp,pos,bpint,em,kmin,accuracy,outfileroot)
 
% %version MartinSchorb 120203
% %
% %usage is martin_tfm_beads(bead coordinates,position of interest,EM image,minimum number of beads used for transformations,acc,output folder);
% %
% % 
% %
% %designed for estimating the accuracy of correlating light and em images using fluorescent electron
% %dense fiducials.
% %
% % 
% output is a struct including
%     for each blind fiducial a set of:  
%     - sel(for each transform):
%       * beads: the picked beads (string)
%       * points: same thing but as vector
%       * ip and bp: their coordinates
%       * tfm: the transformation
%       * stat_used: the distribution of the used beads
%       * bptfm: the predicted coordinates using the transform and
%           indication if point was used for transform.
%       * dist: DISTANCE OF THE ESTIMATED POINT FROM center of used
%           beads
%       * ls: sum of squared residuals (DISTANCE OF THE ESTIMATED positions
%          FROM picked spots) for least squares analysis.
%       * blindtfm: prediction of the blind bead's position
%       * blinddev: deviation of this prediction from picked position
%     - all: the results from the transformation using all beads (as above)
%     -  
% 
% 
%----------------------


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
    [mincl dix]=min(dist1);
    output.closdist=mincl;
    output.dix=dix;
% 
%     % include transformation using all beads
% 
     %calculate current transformation
     tfm=cp2tform(bp2,ip2,'linear conformal');
     output.all.tfm=tfm;           %(*)martin_ls_blind3.m
     output.all.circle=martin_circle(em,accuracy,tformfwd(tfm,round(bpint)));
     %transform coordinates 
     bptfm=tformfwd(tfm,bp2);
     output.all.bptfm=bptfm;
     %calculate squared residue
%      output.blind(nblind).all.ls=0;
%      lsall=0;
%      for r=1:n
%          lsall=lsall+(bptfm(r,1)-ip(r,1))^2+(bptfm(r,2)-ip(r,2))^2;
%      end
%      output.blind(nblind).all.ls=lsall;
     output.all.stat_used=martin_beads_analysis(ip2);
     
%      generate images
%      tfmed=uint8(zeros(size(em)));
%      newColorImage=uint8([size(em),3]));
%      bpr=round(bptfm);
% if (max(max(bpr))<size(em,1)-6) && (min(min(bpr))>6)
%     for n=1:size(bpr,1)
%         tfmed(bpr(n,2)-5:bpr(n,2)+5,bpr(n,1)-5:bpr(n,1)+5)=10;
%         tfmed(bpr(n,2),bpr(n,1))=10;
%     end   
%   
%     newColorImage(:,:,1) =uint8(output.all.circle*255)+em;
%     newColorImage(:,:,2) =tfmed/10*255+em;
%     newColorImage(:,:,3) =pickedem/10*255+uint8(output.all.circle*255)+em;
%     output.all.image=newColorImage;
% else
%    output.all.image=['Bad transformation using all beads as transformation base.'];
%    disp( ['Bad transformation using all beads as transformation base.']);
% end
output.all.closerr=norm(ip2(dix)-output.all.bptfm(dix));
% output.all.tfmed=tfmed;


% prediction loops -------------------------------------------------

for nblind=1:ntot  % blind bead index

    %ip3=[ip2(1:-1,:);999999 9999999;ip2(nblind+1:end,:)];
    ip=ip2;
    bp=bp2;
    
 
   
    ip(nblind,:)=[];
    bp(nblind,:)=[];
    n=size(ip,1); %total number of picked beads
    output.blind(nblind).ip=ip;
    output.blind(nblind).bp=bp;
    m4=150000000*ones(n-kmin+1,max(nchoosek(n+1,kmin+2),nchoosek(n,kmin)));
    
    for k=kmin:(n-1)  % # of trafo base index
        
        tsize=nchoosek(n,k);
        permidx=combnk(1:n,k);
        
        for cnt=1:tsize %index to go through the possible transformations
            output.blind(nblind).sel(k-kmin+1,cnt).beads=permidx(cnt,:);
            
            tip=ip(permidx(cnt,:),:);
            tbp=bp(permidx(cnt,:),:);
            output.blind(nblind).sel(k-kmin+1,cnt).beads=permidx(cnt,:);
           
            output.blind(nblind).sel(k-kmin+1,cnt).ip=tip;             %(*)
            output.blind(nblind).sel(k-kmin+1,cnt).bp=tbp;             %(*)
            
            %analysis of the distribution
            output.blind(nblind).sel(k-kmin+1,cnt).stat_used=martin_beads_analysis2(tip,ip2,nblind);       %(*)
            output.blind(nblind).sel(k-kmin+1,cnt).point=permidx(cnt,:);


             %calculate current transformation
              tfm=cp2tform(tbp,tip,'linear conformal');
             output.blind(nblind).sel(k-kmin+1,cnt).tfm=tfm;           %(*)
              %transform coordinates and estimate
             bptfm=tformfwd(tfm,bp2);

            %mark whether predicted position used to generate transform
            for q=1:k
                bptfm(output.blind(nblind).sel(k-kmin+1,cnt).point(q),3)=1;
            end
              output.blind(nblind).sel(k-kmin+1,cnt).bptfm=bptfm;
%               output.blind(nblind).sel(k-kmin+1).bselect=bselect;

            
            
            %calculate internal property
            ls=0;

                
            ls=sum(sum((bptfm(:,1:2)-ip2).^2));
                
                
            m4(k-kmin+1,cnt)=ls/n;
            
                      
%             calculate accuracy of prediction
            output.blind(nblind).sel(k-kmin+1,cnt).blindtfm=tformfwd(tfm,bp2(nblind,:));
            output.blind(nblind).sel(k-kmin+1,cnt).blinddev=(output.blind(nblind).sel(k-kmin+1,cnt).blindtfm(1)-ip2(nblind,1))^2+(output.blind(nblind).sel(k-kmin+1,cnt).blindtfm(2)-ip2(nblind,2))^2;

   

    end
    
    end  
    

    



% find minimum of squared residues and the corresponding transformation

[C,rows]=min(m4);
if size(rows,2)>1
    [minimum,colmin]=min(C);
    rowmin=rows(colmin);
else
    minimum=C;
    rowmin=1;
    colmin=rows;
end
output.blind(nblind).optimum=[rowmin,colmin];
output.blind(nblind).minimum=minimum;
output.blind(nblind).rowmin=rowmin;
output.blind(nblind).colmin=colmin;
output.blind(nblind).closerr=norm(ip2(dix)-output.blind(nblind).sel(rowmin,colmin).bptfm(dix));
output.blind(nblind).optimtfm=output.blind(nblind).sel(rowmin,colmin).tfm;

% % apply optimal transform to predict all other beads.

% select beads
for ind=1:size(output.blind(nblind).sel(rowmin,colmin).point,2)
    bp(output.blind(nblind).sel(rowmin,colmin).point(end-ind+1),:)=[];
    ip(output.blind(nblind).sel(rowmin,colmin).point(end-ind+1),:)=[];
end

    bpotfm=tformfwd(output.blind(nblind).optimtfm,bp);
    output.blind(nblind).circle=martin_circle(em,accuracy,round(tformfwd(output.blind(nblind).optimtfm,bpint)));
%     
%     error estimation
    output.blind(nblind).devall=minimum;
    output.blind(nblind).preddev=minimum;
    output.blind(nblind).bpotfm=bpotfm;
    
        for eix=1:size(bpotfm,1)
            output.blind(nblind).devall=output.blind(nblind).devall+(norm(bpotfm(eix,:)-ip(eix,:)))^2;
            output.blind(nblind).preddev(eix+1)=(norm(bpotfm(eix,:)-ip(eix,:)))^2;
        end


% devall(nblind)= output.blind(nblind).devall;


%generate images
% tfmed=uint8(zeros(size(em)));
% bpot2=tformfwd(output.blind(nblind).optimtfm,output.blind(nblind).sel(rowmin,colmin).bp);
% bpr=round(bpot2);
% 
% pred=uint8(zeros(size(em)));
% ppr=round([bpotfm;output.blind(nblind).sel(rowmin,colmin).blindtfm]);
% 
% if (max(max([bpr;ppr]))<2042) && (min(min([bpr;ppr]))>6)
% 
%     for n=1:size(bpr,1)
%         tfmed(bpr(n,2)-5:bpr(n,2)+5,bpr(n,1)-5:bpr(n,1)+5)=10;
%         tfmed(bpr(n,2),bpr(n,1))=10;
%     end
% 
% 
% 
%     for n=1:size(ppr,1)
%         pred(ppr(n,2)-5:ppr(n,2)+5,ppr(n,1)-5:ppr(n,1)+5)=10;
%         pred(ppr(n,2),ppr(n,1))=10;
%     end
% 
% 
%     newColorImage(:,:,1) =pred/10*255+uint8(output.blind(nblind).circle*255)+em;
%     newColorImage(:,:,2) =tfmed/10*255+em;
%     newColorImage(:,:,3) =pickedem/10*255+uint8(output.blind(nblind).circle*255)+em;
%     output.blind(nblind).image=newColorImage;
% else
%    output.blind(nblind).image=['blind: ',int2str(nblind),' -- ',output.blind(nblind).sel(rowmin,colmin).beads];
% %    disp( ['Bad transformation using bead ',int2str(nblind),' as blind bead']);disp(['and beads ',output.blind(nblind).sel(rowmin,colmin).beads,' as transformation base.']);
% end

% output.blind(nblind).pred=pred;
% output.blind(nblind).tfmed=tfmed;


% figure(nblind);
% imshow(newColorImage);
% 
% imwrite(pred,[outfileroot,'_',stamp,'_',int2str(ntot),'bds_optfm_',int2str(nblind),'_pred.tif'],'Compression','none');
% imwrite(tfmed,[outfileroot,'_',stamp,'_',int2str(ntot),'bds_optfm_',int2str(nblind),'_tfmed.tif'],'Compression','none');
% imwrite(pickedem,[outfileroot,'_',stamp,'_',int2str(ntot),'bds_optfm_',int2str(nblind),'_pickedem.tif'],'Compression','none');



end 
