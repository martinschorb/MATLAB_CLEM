function [prdr,diffvec] = martin_accuracy(kmin,files)

% %version MartinSchorb 130620
% %
% % Copyright EMBL 2013, All rights reserved



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
corrindex = 1;
while x>0 & stop==0
pause(0.001)
if nargin ~= 2
%import previously clicked positions
[filename, pathname] = uigetfile({'*.pickspots1.mat'},'select previously picked beads','MultiSelect', 'on',in_dir);
else
   pathname=[];filename=files;stop=1; 
end


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

            currmax = n;
        
% pt of interest....

pt1 = ip2(nblind,:);    

                 tip=ip;
                 tbp=bp;
                 st = size(tip,1);
                 
                 tip2=tip;
                 tbp2=tbp;
                
                 %calculate current transformation
                  tfm=cp2tform(tbp2,tip2,trafo);
                  bptfm2=tformfwd(tfm,bp2);                 

                usedbeads(corrindex) = st;
                
                prederr = norm(ip2(nblind,:)-bptfm2(nblind,:))*pxs;
                
                             
                prdr(corrindex)=prederr;
                
                diffvec(corrindex,:) = (ip2(nblind,:)-bptfm2(nblind,:))*pxs;
                
              %counter for data file
               corrindex = corrindex+1;

        end
  end
  end
end


