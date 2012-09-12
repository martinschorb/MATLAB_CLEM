function data=martin_corr_accuracy(kmin)

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


while x>0
pause(0.001)

%import previously clicked positions
[filename, pathname] = uigetfile({'*.pickspots1.mat'},'select previously picked beads','MultiSelect', 'on',in_dir);

if isequal(filename,0)
    disp('Correlation accuracy estimation finished');
    x=0;    
else
    
    
for fileidx=1:size(filename,2) %file list index

        a=open([pathname,cell2mat(filename(fileidx))]);

 %  analyse blind bead accuracy for the current file (image&bead coordinates)
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

    data.file(y)=filename(fileidx);
    
        n = ntot-1; %total number of picked beads
        
        if n>kmin
%         m1=150000000*ones(n-kmin+1,max(nchoosek(n+1,kmin+2),nchoosek(n,kmin)));
%         m2=150000000*ones(n-kmin+1,max(nchoosek(n+1,kmin+2),nchoosek(n,kmin)));
%         m3=150000000*ones(n-kmin+1,max(nchoosek(n+1,kmin+2),nchoosek(n,kmin)));
        m4=150000000*ones(n-kmin+1,max(nchoosek(n+1,kmin+2),nchoosek(n,kmin)));
%         m5=150000000*ones(n-kmin+1,max(nchoosek(n+1,kmin+2),nchoosek(n,km
%         in)));
        else
            m4=1500000000;
        end
            
        for k=kmin:(n)  % # of trafo base index
            tsize=nchoosek(n,k);
            permidx=combnk(1:n,k);


            for cnt=1:tsize %index to go through the possible transformations
                
             
                 tip=ip(permidx(cnt,:),:);
                 tbp=bp(permidx(cnt,:),:);

                 output.blind(nblind).sel(k-kmin+1,cnt).ip=tip;             %(*)
                 output.blind(nblind).sel(k-kmin+1,cnt).bp=tbp;             %(*)


                 %calculate current transformation
                  tfm=cp2tform(tbp,tip,trafo);

                  %transform coordinates and estimate
                 bptfm2=tformfwd(tfm,bp2);
                 bptfm=tformfwd(tfm,bp);
                 
%                 data.closerr(x)=norm(bptfm2(dix)-ip2(dix));
                

     %calculate internal deviations

%      optional additional exponent for weighting
        expo=1;

                 ls=sum(sum(((bptfm(:,1:2)-ip).^2).^expo,2));% least squares deviation for all beads

% 
%                 if ls/n<2000                 
%                   m1(k-kmin+1,cnt)=data.closerr(x)*(ls/n)^(0.5);
%                 end
%                 m2(k-kmin+1,cnt)=data.closdist(x);
%                 m3(k-kmin+1,cnt)=data.closerr(x)*ls;
                m4(k-kmin+1,cnt)=ls/n;

%                data.prederr(x)=norm(ip2(nblind,:)-bptfm2(nblind,:)); 
   
%                m5(k-kmin+1,cnt)=data.prederr(x);
               
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
        
        sbp=output.blind(nblind).sel(rowmin,colmin).bp;
        sip=output.blind(nblind).sel(rowmin,colmin).ip;
        
        tfm=cp2tform(sbp,sip,trafo);
        bptfm2=tformfwd(tfm,bp2);
        
        
%         data.optimum(fileidx,nblind)={[rowmin,colmin]};

        data.tfmbeads(y)=size(sbp,1);
        
        data.imbeads(y)=ntot;
        data.analysis(y).ip=ip;
        data.analysis(y).bp=bp;
        data.analysis(y).blindpos=ip2(nblind,:);
   
        data.ls(y)=m4(rowmin,colmin);
        
        data.optprederr(y)=norm(ip2(nblind,:)-bptfm2(nblind,:));
        y=y+1;
    end
        
 

    end
    end
    
    end
end

data.deviation = sortrows(data.optprederr'*pxs);
lp = 1/length(data.deviation);
data.percentage = lp:lp:1;





