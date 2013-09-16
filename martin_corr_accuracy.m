function data = martin_corr_accuracy(kmin,files)

% %version MartinSchorb 120912
% %


% accthr = 50/(5.068/2);

if nargin==0
    kmin=3;
end

if exist('corr_init')==2
    corr_init();
elseif exist('corr_init_orig')==2
    corr_init_orig(); 
else 
    a=msgbox('No initialization script found!','Error','modal');uiwait(a);
    a=msgbox('Please update algorithms!','Error','modal');uiwait(a);
    return 
end 


% file_1 = fopen('/struct/briggs/schorb/batchcorr/images_good.csv','w');
% fprintf(file_1,['Correlation accuracy for certain images.\n\n threshold is set to ',num2str(accthr*5.068),' nm\n\n File    ;   number of total beads    ;  outimage numbers  ;  number of instances above threshold  \n------------------------\n']);
    
  

%  initialization and GUI

[outfile,in_dir,pxs,trafo,minbeads,maxbeads,maxdist]=martin_corr_accuracy_init(loc_pickspots,pixelsize_lm,trafo,kmin);

if isempty(outfile) | isempty(minbeads) | isempty(pxs)
    error('Please provide initial values');
end


x=1;
y=1;
z=1;
oo=1;
mult = 0;

kmin=minbeads;
data=struct;
stop=0;

while x>0 & stop==0
pause(0.001)
if nargin ~= 2
%import previously clicked positions
[filename, pathname] = uigetfile({'*.pickspots1.mat'},'select previously picked beads','MultiSelect', 'on',in_dir);
else
   pathname=[];filename=files;stop=1; 
end

% data.numcorr=zeros(6,1);

if isequal(filename,0)
    disp('Correlation accuracy estimation finished');
    x=0;    
else
    
    
    
for fileidx=1:size(filename,2) %file list index
        disp(['processing ',filename{fileidx},' ... ',num2str(fileidx),' of ',num2str(size(filename,2))]);
                data.files{fileidx} = cell2mat(filename(fileidx));
                
        a=open([pathname,cell2mat(filename(fileidx))]);
% 
% %%  analyse blind bead accuracy for the current file (image&bead coordinates)
        ip=a.ip;
        bp=a.bp;
        ip2=ip;
        bp2=bp;
        ntot=size(ip,1);
        totnum(fileidx)=ntot;
        
%         maximum number of beads
   
        
%         radius cleaning

alldists=squareform(pdist(ip));

        
        
        
if ntot<12 & ntot > kmin
%         nblind=1;
%         alli=y;
%         


for nblind=1:ntot  % blind bead index

    ip=ip2;
    bp=bp2;
    
    
   if isempty(maxdist)
       selbeads = [];
   else
    dists=alldists(:,nblind);
    selbeads = find(dists>(maxdist/pxs*1000));
   end
       
    
    
    
    
    
   selection=[nblind;selbeads];
    ip(selection,:)=[];
    bp(selection,:)=[]; 
    
     n = size(ip,1); %total number of picked beads
     if n<kmin
         continue
     end
     
     
%         if isempty(maxbeads)
            currmax = n;
%         else
%             currmax = min(maxbeads,n);
%         end
        
% pt of interest....
% 
pt1 = ip2(nblind,:);    
% pt2 = bp2(nblind,:); 


% find fiducial closest to POI

% dist = sum((ip-repmat(pt1,[n,1])).^2,2);
% 
% [mincl dix]=min(dist);   
% 
% dist2=dist;
% dist2(dix)=[];
% 
% [mincl2 dix2]=min(dist2);  

% dist_em = sqrt(sum((ip-repmat(center_em,[n 1])).^2,2));
% dist_fm = sqrt(sum((bp-repmat(center_fm,[n 1])).^2,2));

% scale = dist_em./dist_fm;

 ip1=ip;
 bp1=bp;

% alltfm=cp2tform(bp,ip,trafo);
% allbptfm=(tformfwd(alltfm,bp2));
% data.allprederr(y)=norm(allbptfm(nblind,:)-ip2(nblind,:));
% allprederr(nblind)=data.allprederr(y);
 
 corrindex=1;
%  data.corr(z).blind(nblind).ipall=ip;
%  data.corr(z).blind(nblind).bpall=bp;

        for k=kmin:currmax  % # of trafo base index
            tsize=nchoosek(currmax,k);
            permidx=combnk(1:currmax,k);

            
            for cnt=1:tsize %index to go through the possible transformations
                
             
                 tip=ip(permidx(cnt,:),:);
                 tbp=bp(permidx(cnt,:),:);
                 st = size(tip,1);

                  
                 tip2=tip;
                 tbp2=tbp;
                 
  
                
                 %calculate current transformation
                  tfm=cp2tform(tbp2,tip2,trafo);
%                   data.corr(z).blind(nblind).trafo{corrindex} = tfm;
%                   if isfield(tfm.tdata,'T')
%                   data.corr(z).blind(nblind).T(:,:,corrindex) = tfm.tdata.T;
%                   else
%                      data.corr(z).blind(nblind).T(:,:,corrindex) = tfm.tdata.tshifted.tdata.T;disp('shift thing...');
%                   end
%                   %transform coordinates and estimate
%                   bptfm1=tformfwd(tfm,bp1);      
                  bptfm2=tformfwd(tfm,bp2);                 
%                   bptfm=tformfwd(tfm,tbp);
% 
%      optional additional exponent for weighting
%         expo = 1;
%         expo_r = 1;
%                  ls = sum(sum(((bptfm(:,1:2)-tip).^2).^expo,2)); % least squares deviation for all beads
%                  lsr = sum((1./dist_em1).^expo_r.*sum(((bptfm(:,1:2)-tip).^2).^expo,2));

                usedbeads(corrindex) = st;
%                 m4(k-kmin+1,cnt)=ls/usedbeads;
                
                prederr = norm(ip2(nblind,:)-bptfm2(nblind,:));
                data.blinddev(x,:)=ip2(nblind,:)-bptfm2(nblind,:);
                data.prederr(x,1)=prederr;
                
                data.k_all(x)=k;
                
                prdr(corrindex)=prederr;
                
              %counter for data file
               x=x+1;
               corrindex = corrindex+1;
            end
            if cnt==1
                data.allbeads_err(y) = prederr;
                data.allbeaddiffv(y,:) = ip2(nblind,:)-bptfm2(nblind,:);
            end

        end

%        trafoind=1:corrindex-1;
       
[testdist,tix] = min(prdr);

data.n_used(y) = usedbeads(tix);


data.minprederr(y)=testdist;

meandist = mean(prdr);
data.meanprederr(y)=meandist;
data.n(y)=n;
clear prdr

% 
% if testdist < accthr
      file = filename{fileidx}(1:end-15);
      
%       disp(num2str(testdist));
%       dirfind = strfind(file,'/');
%       file = file(dirfind(end)+1:end);
      
%     file1 = file;
%     
%     lmdate = martin_dbread(file,8,'/struct/briggs/schorb/_Endo-Data/outliers.csv');
%     
%     if isempty(lmdate)
%     
%     if ~isempty(str2num(file(1))) 
%        if ~isempty(str2num(file(2))) 
%         file = [file(3),file(1:2),file(4:end)];
%        else
%         file = [file(2),file(1),file(3:end)];
%        end
%         
%     end
%         lmdate = martin_dbread(file,8,'/struct/briggs/schorb/_Endo-Data/outliers.csv');
% 
%     end
%     
%     
%     if ~isempty(str2num(file(end))) 
%         filebase = file;
%     else
%         filebase = file(1:end-1);
%     end
%     
    
    
%     check for multiplicity
%     if mult==0
%             
%             fprintf(file_1,[file1,'  ;  ',num2str(ntot),'  ;  (',num2str(oo)]);        
%     else
%             fprintf(file_1,[', ',num2str(oo)]);
%     end
%         
%     mult = mult+1;
    
    
%     hard work starts here....



%     if exist(['/struct/briggs/schorb/accuracy/outliers/good-',num2str(oo),'.jpg'])<8
%     
%   
%     
%     
%     if isempty(lmdate)
%         disp('emergency stop');
%         pause
%     end
%     
%     if iscell(lmdate)
%         lmdate = lmdate{1}(2:end-1);
%     end
%     
%     lmdate=num2str(lmdate);
%        
%     while isempty(str2num(lmdate(1)))
%         lmdate = lmdate(2:end);
%     end
%     
%     if str2num(lmdate(1))>1
%         lmdate = ['0' lmdate];
%     end
% 
%     cd(['/struct/briggs/wanda/DataLightMicroscopy/',lmdate]);
%     
%     fmd = dir([filebase(1:3),'*']);
%     
%     if isempty(fmd)
%         fmd = dir([file1(1:3),'*']);
%     end
%     
%     if isempty(fmd)
%         fmd = dir(['*',file1(1:2),'*']);
%     end
%     
%      if isempty(fmd)
%         fmd = dir(['*',file(1:2),'*']);
%     end
%     
%     if isempty(fmd)
%         fmd = dir(['*',file1(2),file1(1),'*']);
%     end
%     
%     
%     if isempty(fmd)
%         disp('emergency stop');
%         pause
%     end
%     
%     if isempty(strfind(fmd(1).name,'tif'))        
%         cd(fmd(1).name);
%     end
%     
% %     fm = imread(a.fmf);
%     
%     cd corr;
%     
%     emd = dir([file1,'*','_em.tif']);
%     
%      if isempty(emd)
%        emd = dir([filebase,'*','_em.tif']);
%      end
%        
%      
%     if isempty(emd)
%         disp('emergency stop');
%         pause
%     end
%     
%     
%     em = imread(emd(1).name);
%     
%     f=figure('visible','off');
%     imagesc(em)
%     axis off
%     colormap('gray')
%     axis equal
%     hold on
% %     scatter(ip2(nblind,1),ip2(nblind,2),4404,'rx')
%     scatter(ip(:,1),ip(:,2),442,'co')
%     set(f,'Position',[16000 1600 1600 1200])
%     
%     sp1 = strfind(file1,'_');
%     file1 = [file1(1:sp1),'{}-',file1(sp1+1:end)];
%     
%    tt = text(0,0,[file1]);%'Minimum error: ',num2str(round(testdist*5.068)),' nm   -   ',file1]);
%      
%     cd  /struct/briggs/schorb/accuracy/outliers;
%     
% %     saveas(f,['good-',num2str(oo),'.jpg']);
%     saveas(f,[file1,'.jpg']);
% end
%     oo=oo+1;
%     close all
%     
% end



        y=y+1;
    end
      z=z+1;  

 
end
%     if mult>0
%     fprintf(file_1,[')  ;  ',num2str(mult),'   ;   \n']);
%     mult = 0;
%     end

    end
    
     end
end

data.allerrors = sortrows(data.prederr*pxs)
 
data.besterrors = sortrows(data.minprederr'*pxs);

 lpa = 1/length(data.allerrors);
 data.allpercentage = lpa:lpa:1;
 
 lpb = 1/length(data.besterrors);
 data.bestpercentage = lpb:lpb:1;

%     fclose(file_1);


