function [data,output]=martin_blind_errorpred3(kmin,outfileroot)

% %version MartinSchorb 100129
% %
% %usage is martin_blindstat('outputfileroot');
% %
% %
% % uses martin_ls_blind.m for calculating estimations and ratings of
% %     transformations
% % designed for estimating the accuracy of correlating light and em images using fluorescent electron
% % dense fiducials.
% % 
% % 
% % 
% % output is a struct including:
% % 
% % xx - filename pickspots
% % 1 - index of blind bead
% % 2 - index of optimal transformation 
% % 3 - error of blind prediction
% % 4 - error of closest bead
% % 5 - relative error for prediction


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

%save data file output
% 
% file_1 = fopen([outfileroot,'blind_errorpred_',stamp,'.txt'],'w');
% fprintf(file_1,'Statistics for blind bead prediction errors \n\n');   




x=1;
y=1;
data=struct;
data.file=cell(2);
data.optimum=cell(2);


while x>0
%import previously clicked positions
[filename, pathname] = uigetfile({'*.mat'},'select previously picked beads','MultiSelect', 'on','/struct/briggs/schorb/');
if isequal(filename,0)
    disp('thanks...');
    x=0;    
else
%     fprintf(file_1,['Filenames: \n', pathname,'\n-------\n']);
    for fileidx=1:size(filename,2) %file list index
%         fprintf(file_1,[cell2mat(filename(fileidx)),'\n']);
        a=open([pathname,cell2mat(filename(fileidx))]);
        n=size(a.ip,1);         % # of beads

    %          %  analyse blind bead accuracy for the current file (image&bead coordinates)
        ip=a.ip;
        bp=a.bp;

        ip2=ip;
        bp2=bp;
        ntot=size(ip,1);
if ntot>4 & ntot<12
        output=struct;
        nblind=1;
    
        alli=y;
        
for nblind=1:ntot  % blind bead index
    ip3=[ip2(1:nblind-1,:);999999 9999999;ip2(nblind+1:end,:)];
    ip=ip2;
    bp=bp2;
    
    dist1=[0 0];        
    for cl=1:ntot
        dist1(cl)=norm(ip2(nblind,:)-ip3(cl,:));
    end
    [mincl dix]=min(dist1);
    data.closdist(x)=mincl;
    selection=[nblind]
    ip(selection,:)=[];
    bp(selection,:)=[]; 

    data.file(y)=filename(fileidx);
    
    alltfm=cp2tform(bp,ip,'linear conformal');
    allbptfm=(tformfwd(alltfm,bp2));
    
    data.allprederr(y)=norm(allbptfm(dix)-ip2(dix));
    allprederr(nblind)=data.allprederr(y);
    
    data.optstat(y,:)=martin_beads_analysis2(ip,ip2,nblind);
    
        n=size(ip,1); %total number of picked beads
        m1=150000000*ones(n-kmin+1,max(nchoosek(n+1,kmin+2),nchoosek(n,kmin)));
        m2=150000000*ones(n-kmin+1,max(nchoosek(n+1,kmin+2),nchoosek(n,kmin)));
        m3=150000000*ones(n-kmin+1,max(nchoosek(n+1,kmin+2),nchoosek(n,kmin)));
        m4=150000000*ones(n-kmin+1,max(nchoosek(n+1,kmin+2),nchoosek(n,kmin)));
        m5=150000000*ones(n-kmin+1,max(nchoosek(n+1,kmin+2),nchoosek(n,kmin)));
        
        for k=kmin:(n)  % # of trafo base index

            %generate binary selector(permutation) for beads

            a=[2^k-1:2^(n)-1];
            bm=dec2bin(a);  %binary matrix
            sel=zeros(size(bm)); %selection matrix
            lim=size(bm,1);

            for i=1:lim
                for j=1:size(bm,2)
                    sel(i,j)=str2num(bm(i,j));
                end
            end

            %remove all combinations involving less than k beads
            bselect=zeros(n,n);
            m=1;
            for l=1:lim
                if sum(sel(l,:))==k
                bselect(m,:)=sel(l,:);  %final selection matrix
                m=m+1;
                end
                l=l+1;
            end

            %calculate transforms and predicted positions
            tbp=zeros(k,2);
            tip=zeros(k,2);

            for cnt=1:size(bselect,1) %index to go through the possible transformations
                
             
                data.closdist(x)=mincl;
                
                output.blind(nblind).sel(k-kmin+1,cnt).beads='';
                p=1;
                for o=1:n
                    if bselect(cnt,o)~=0
                        tip(p,:)=ip(o,:);
                        tbp(p,:)=bp(o,:);
                        output.blind(nblind).sel(k-kmin+1,cnt).beads=[output.blind(nblind).sel(k-kmin+1,cnt).beads,' ',int2str(o)]; %(*)
                        output.blind(nblind).sel(k-kmin+1,cnt).point(p)=o;             %(*)
                        p=p+1;
                    end
                end
                output.blind(nblind).sel(k-kmin+1,cnt).ip=tip;             %(*)
                output.blind(nblind).sel(k-kmin+1,cnt).bp=tbp;             %(*)

                %analysis of the distribution
                output.blind(nblind).sel(k-kmin+1,cnt).stat_used=martin_beads_analysis2(tip,ip2,nblind);       %(*)
                


                 %calculate current transformation
                  tfm=cp2tform(tbp,tip,'linear conformal');
                 output.blind(nblind).sel(k-kmin+1,cnt).tfm=tfm;           %(*)
                  %transform coordinates and estimate
                 bptfm2=tformfwd(tfm,bp2);
                 bptfm=tformfwd(tfm,bp);
                 
                data.closerr(x)=norm(bptfm2(dix)-ip2(dix));
                

     %calculate internal deviations
                ls=0;
%                 for r=1:n
%                    ls=ls+((bptfm(r,1)-ip(r,1))^2+(bptfm(r,2)-ip(r,2))^2)/sqrt(norm(ip(r,:)-ip2(nblind)));
%                 end
expo=1;

                 ls=sum(sum(((bptfm(:,1:2)-ip).^2).^expo,2));%./((sum((ip-repmat(ip2(nblind,:),[length(ip) 1])).^4,2))));



                if ls/n<2000                 
                  m1(k-kmin+1,cnt)=data.closerr(x)*(ls/n)^(0.5);
                end
                m2(k-kmin+1,cnt)=data.closdist(x);
                m3(k-kmin+1,cnt)=data.closerr(x)*ls;
                m4(k-kmin+1,cnt)=ls/n;

               data.prederr(x)=norm(ip2(nblind,:)-bptfm2(nblind)); 
   
               m5(k-kmin+1,cnt)=data.prederr(x);
               
               %counter for data file
               x=x+1;
            end
            


        end

% 
%         
%         for i=1:size((output.blind(nblind).sel),1)
%     for j=1:nchoosek(ntot-1,i+4)
%     
%     m5(i,j)=m1(i,j)*m2(i,j);
%     
%     
%     end
%   
% end
        
        
%         
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
        
        sbp=[output.blind(nblind).sel(rowmin,colmin).bp;bp2(dix,:)];
        sip=[output.blind(nblind).sel(rowmin,colmin).ip;ip2(dix,:)];
        
        tfm=cp2tform(sbp,sip,'linear conformal');
        bptfm2=tformfwd(tfm,bp2);
        
        
        data.optimum(fileidx,nblind)=mat2cell([rowmin,colmin]);
        data.analysis(y).stat=output.blind(nblind).sel(rowmin,colmin).stat_used;
        data.clouddist(y)=norm(output.blind(nblind).sel(rowmin,colmin).stat_used(4:5)-ip2(nblind,:));
        data.cloudell(y)=output.blind(nblind).sel(rowmin,colmin).stat_used(3);
        data.relclouddist(y)=data.clouddist(y)/output.blind(nblind).sel(rowmin,colmin).stat_used(2);
        data.tfmbeads(y)=output.blind(nblind).sel(rowmin,colmin).stat_used(1);
        
        data.imbeads(y)=ntot;
        data.analysis(y).ip=ip;
        data.analysis(y).bp=bp;
        data.analysis(y).blindpos=ip2(nblind,:);
        data.analysis(y).tip=output.blind(nblind).sel(k-kmin+1,cnt).ip;
        data.analysis(y).tbp=output.blind(nblind).sel(k-kmin+1,cnt).bp;
        data.optclosdist(y)=m2(rowmin,colmin);
        data.optcloserr(y)=m3(rowmin,colmin);
        
        if data.optcloserr(y)>200
            
        end
        data.ls(y)=m4(rowmin,colmin);
        data.optprederr(y)=norm(ip2(nblind)-bptfm2(nblind));
        y=y+1;
    end
        
     
    allmean=mean(allprederr);
    allsigma=std(allprederr);
    for all=alli:y-1
        data.allprederr(all)=data.allprederr(all)/(allmean+allsigma);
    end
    end
    end
    
    end
end

%  fclose(file_1);
