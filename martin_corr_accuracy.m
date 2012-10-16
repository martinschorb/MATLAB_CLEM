function data=martin_corr_accuracy(kmin,files)

% %version MartinSchorb 120912
% %

if exist('corr_init')==2
    corr_init();
elseif exist('corr_init_orig')==2
    corr_init_orig(); 
else 
    a=msgbox('No initialization script found!','Error','modal');uiwait(a);
    a=msgbox('Please update algorithms!','Error','modal');uiwait(a);
    return 
end 


    
   

%  initialization and GUI

[outfile,in_dir,pxs,trafo,minbeads]=martin_corr_accuracy_init(loc_pickspots,pixelsize_lm,trafo);

if isempty(outfile) | isempty(minbeads) | isempty(pxs)
    error('Please provide initial values');
end


x=1;
y=1;

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

data.numcorr=zeros(6,1);

if isequal(filename,0)
    disp('Correlation accuracy estimation finished');
    x=0;    
else
    
    
for fileidx=1:size(filename,2) %file list index
        disp(['processing ',filename{fileidx},' ... ',num2str(fileidx),' of ',num2str(size(filename,2))]);
        a=open([pathname,cell2mat(filename(fileidx))]);

%%  analyse blind bead accuracy for the current file (image&bead coordinates)
        ip=a.ip;
        bp=a.bp;
        ip2=ip;
        bp2=bp;
        ntot=size(ip,1);
        
if ntot>minbeads & ntot<12
        nblind=1;
        alli=y;
        
for nblind=1:ntot  % blind bead index
    ip3=[ip2(1:nblind-1,:);999999 9999999;ip2(nblind+1:end,:)];
    ip=ip2;
    bp=bp2;
    
    
% find fiducial closest to the feature of interest (here: to blind bead)

%     dist1=[0 0];        
%     for cl=1:ntot
%         dist1(cl)=norm(ip2(nblind,:)-ip3(cl,:));
%     end
%     [mincl dix]=min(dist1);
%     data.closdist(y)=mincl;

    selection=[nblind];
    ip(selection,:)=[];
    bp(selection,:)=[]; 
    
     n = ntot-1; %total number of picked beads
     
     
% distance measurements:    
    
center_em = mean(ip);    
center_fm = mean(bp); 


% pt of interest....
% 
pt1 = ip2(nblind,:);    
% pt2 = bp2(nblind,:); 

% find fiducial closest to POI

dist = sum((ip-repmat(pt1,[n,1])).^2,2);

[mincl dix]=min(dist);   

% dist_em = sqrt(sum((ip-repmat(center_em,[n 1])).^2,2));
% dist_fm = sqrt(sum((bp-repmat(center_fm,[n 1])).^2,2));

% scale = dist_em./dist_fm;
%     
%     
    


 data.file(y)=filename(fileidx);
 alltfm=cp2tform(bp,ip,trafo);
 allbptfm=(tformfwd(alltfm,bp2));
 allclosetfm=tformfwd(alltfm,bp);
 allclosetfm2=tformfwd(alltfm,bp2);
 
 data.allprederr(y)=norm(allbptfm(nblind,:)-ip2(nblind,:));   
       
 data.allcloseerr(y) = norm(allclosetfm(dix,:)-ip(dix,:)); 
 data.allcloseprederr(y) = norm(allclosetfm2(nblind,:)-ip2(nblind,:));  
%  n=n-1;
 ip1=ip;
 bp1=bp;
%  
%  ip(dix,:)=[];
%  bp(dix,:)=[];
 
%         if n>kmin
% %         m1=150000000*ones(n-kmin+1,max(nchoosek(n+1,kmin+2),nchoosek(n,kmin)));
% %         m2=150000000*ones(n-kmin+1,max(nchoosek(n+1,kmin+2),nchoosek(n,kmin)));
% %         m3=150000000*ones(n-kmin+1,max(nchoosek(n+1,kmin+2),nchoosek(n,kmin)));
%         m4=150000000*ones(n-kmin+1,max(nchoosek(n+1,kmin+2),nchoosek(n,kmin)));
% %         m5=150000000*ones(n-kmin+1,max(nchoosek(n+1,kmin+2),nchoosek(n,km
% %         in)));
%         else
%             m4=1500000000;
%         end

        for k=kmin:(n)  % # of trafo base index
            tsize=nchoosek(n,k);
            permidx=combnk(1:n,k);
            

            for cnt=1:tsize %index to go through the possible transformations
                
             
                 tip=ip(permidx(cnt,:),:);
                 tbp=bp(permidx(cnt,:),:);
                 
                 output.blind(nblind).sel(k-kmin+1,cnt).ip=tip;             %(*)
                 output.blind(nblind).sel(k-kmin+1,cnt).bp=tbp;             %(*)
                 
                 
           
                 st = size(tip,1);
%                  dist_em1 = sqrt(sum((tip-repmat(center_em,[st 1])).^2,2));
%                  [sorted,sidx] = sort(dist_em1);
%                  weight = st-sidx+1;
                  
                 tip2=tip;
                 tbp2=tbp;
                 
%                  tip2=[];
%                  tbp2=[];
%                  
%                  for i=1:st
%                      tip2 = [tip2;repmat(tip(i,:),[weight(i)^2,1])];
%                      tbp2 = [tbp2;repmat(tbp(i,:),[weight(i)^2,1])];
%                  end
%                                   
                
                 %calculate current transformation
                  tfm=cp2tform(tbp2,tip2,trafo);
                  data.corr(y).tfm(cnt) = tfm;
                  %transform coordinates and estimate
                  bptfm1=tformfwd(tfm,bp1);      
                  bptfm2=tformfwd(tfm,bp2);                 
                  bptfm=tformfwd(tfm,tbp);
                 
                 data.closerr(x)=norm(bptfm1(dix,:)-ip1(dix,:));
                

     %calculate internal deviations

%      optional additional exponent for weighting
        expo = 1;
%         expo_r = 1;
                 ls = sum(sum(((bptfm(:,1:2)-tip).^2).^expo,2)); % least squares deviation for all beads
%                  lsr = sum((1./dist_em1).^expo_r.*sum(((bptfm(:,1:2)-tip).^2).^expo,2));
%                  
%                  
%                  scale_error = 2;
%                  
% 
%                 if ls/n<2000                 
%                   m1(k-kmin+1,cnt)=data.closerr(x)*(ls/n)^(0.5);
%                 end
%                 m2(k-kmin+1,cnt)=data.closdist(x);
%                 m3(k-kmin+1,cnt)=data.closerr(x)*ls;

                usedbeads = st;
                m4(k-kmin+1,cnt)=ls/usedbeads;

%                 m4(k-kmin+1,cnt)=data.closerr(x);

                data.prederr(x,1)=norm(ip2(nblind,:)-bptfm2(nblind,:)); 
                prdr(cnt)= data.prederr(x,1);
                
                data.prederr(x,2)=n;
                if n>length(data.numcorr)
                    data.numcorr(n)=1;
                else
                    data.numcorr(n)=data.numcorr(n)+1;
                end
                
%                m1(k-kmin+1,cnt)=data.prederr(x,1);
               
               %counter for data file
               x=x+1;
            end
            


        end

        
        % find minimum of ...
        
        [C,rows]=min(m4);
        if size(rows,2)>1
            [minimum,colmin]=min(C);
            rowmin=rows(colmin);
        else
            minimum=C;
            rowmin=1;
            colmin=rows;
        end
%         
        sbp=output.blind(nblind).sel(rowmin,colmin).bp;
        sip=output.blind(nblind).sel(rowmin,colmin).ip;
%         
%         tfm=cp2tform(sbp,sip,trafo);
        
        data.opttfm = output.blind(nblind).sel(rowmin,colmin).tfm;
        bptfm2=tformfwd(tfm,bp2);
        
        data.avtfm(y) = cp2tform([bp1;sbp],[ip1;sip],'linear conformal');
        data.averr(y) = norm(ip2(nblind,:)-tformfwd(data.avtfm(y),bp2(nblind,:)));
        
%         data.optimum(fileidx,nblind)={[rowmin,colmin]};

        data.tfmbeads(y)=size(sbp,1);
        
        data.imbeads(y)=ntot;
        data.analysis(y).ip=ip;
        data.analysis(y).bp=bp;
        data.analysis(y).blindpos=ip2(nblind,:);
        
%         data.optimum(y)=min(min(m1));
        
        data.ls(y)=m4(rowmin,colmin);
        
        data.optprederr(y)=norm(ip2(nblind,:)-bptfm2(nblind,:));
        y=y+1;
    end
        
 

    end
    end
    
    end
end

data.errors = sortrows(data.prederr'*pxs);

data.deviation = sortrows(data.optprederr'*pxs);
lp = 1/length(data.deviation);
data.percentage = lp:lp:1;





