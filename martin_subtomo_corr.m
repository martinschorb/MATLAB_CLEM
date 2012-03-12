% function subpos2=martin_subtomo_corr(subtf,outfileroot)

% version MartinSchorb 100713
% 
% usage is martin_subtomo_corr2(slice of interest,subtomogram filename in mrc-format)
% 
% correlates highmag-fluorescence spot onto selected subtomogram.
% 
% uses cross-correlation to find matching tomogram slice.



% 
% st1=tom_mrcread(subtf);
% subtomo=uint8(st1.Value+128);

cc=who;
for jj=1:length(cc)
    clear(cc{jj});
end

[init,db]=martin_dbread(2,1);

if init~=1
    error('MATLAB:martin_subtomo_corr','check database!');
end



[filename, pathname1] = uigetfile({'sliceinfo.mat'},'select highmag correlation','MultiSelect', 'on','/struct/briggs2/schorb/tfm_logs/sliceinfos');

outfileroot=['/home/schorb/m/110110/stc_1/'];


  keep=[who;'keep'];  


for fileidx=1:size(filename,2)
    

filename1=filename{fileidx};



basefind=strfind(filename1,'_all_');
namebase = filename1(1:basefind-1);
cd(pathname1)
cd('../LMHMcoos')
load([namebase,'.lmhmcoos.mat']);

fluo1=strfind(filename1,'FP');
fluotxt=['_',filename1(fluo1(1)-1:fluo1(1)+1)];


subtfmbase=filename1(1:fluo1(2)-3);

namebase=filename1(1:fluo1(1)-3);


dbindex=martin_dbread(namebase,1,db);

if ~isempty(dbindex) & isequal(martin_dbread(dbindex+1,12),{'"inv"'})

subtf=['/struct/briggs2/schorb/tomograms/subtomos/',filename1(1:fluo1(1)-3),'.mrc'];

if exist(subtf)>0
    subtomo=uint8(readMRCfile(subtf));
else
    [stfile, stpath] = uigetfile({'mrc'},['select subtomogram - file: ',subtfmbase(1:end-4)],'/struct/briggs2/schorb/tomograms/subtomos');
    subtomo=uint8(readMRCfile([stpath,stfile]));
end

        load([pathname1,filename1]);

if exist(['/home/schorb/m/tfm_logs/subtomocorr/',subtfmbase,'.subttfm.mat'])>0
            load(['/home/schorb/m/tfm_logs/subtomocorr/',subtfmbase,'.subttfm.mat'])
            disp(['existing stf loaded  - ',namebase])

    if slice>0
            slicerange=slice;
    else
        slicerange = floor(0.1*size(subtomo,3)):ceil(0.9*size(subtomo,3));
    end


else
    cd /home/schorb/m/tfm_logs/subtomocorr/
    ss=dir(['*',namebase,'*.subttfm.mat']);
if ~isempty(ss)
    load(ss(1).name)
    slicerange=slice;
else
    slicerange = floor(0.1*size(subtomo,3)):ceil(0.9*size(subtomo,3));
end
            disp(namebase)
            disp(['slice: ',num2str(slicerange)])



        endpos=strfind(filename1,'.slicei');
        % a=readtext([pathname1,filename1(1:endpos-1),'_hm_transform.log']);

        % if exist(subtf)
        % %     
        % % smfstart=strfind(a{3},'interest:');
        % % 
        % % smf=[a{3}(smfstart+10:end)];
        % % 
        % % if ~exist(smf)
        % % smf=['/struct/briggs/wanda',a{3}(smfstart+18:end)];
        % % if ~exist(smf)
        % % smf=['/struct/briggs/wanda/',a{3}(smfstart+16:end)];
        % % end
        % end

        slicerange = floor(0.1*size(subtomo,3)):ceil(0.9*size(subtomo,3));


end


if exist('smf_sub')
    smf=smf_sub;
end


if smf(8:9)=='./'
    cd('/struct/briggs/wanda/DataLightMicroscopy/091216/corr');
else
    cd('/struct/briggs/wanda/DataLightMicroscopy/091216');
end

sm=imread(smf);

% 
% spotfind=strfind(a{5},'spot:');
% 
% spotlog=a{5}(spotfind+5:spotfind+12);
% 
% spot1=str2num(spotlog);
% 
% spotlog=a{5}(spotfind+14:end);
% 
% spot1(2)=str2num(spotlog);
% 
% 
% if sum(spot1-spotpos)>1
%     disp(['warning: different positions deteceted for: ',filename1]);
% end

% k=msgbox('Select region of interest.','To do...','modal');
% uiwait(k);
% [sub_sm,rect_sm] = imcrop(sm);close(gcf)
% while sum(size(sub_sm)< [size(subtomo,1),size(subtomo,2)])>0
%     k=msgbox('Size of selected region must be bigger than subtomogram area!','Error','modal');
%     uiwait(k);
%     [sub_sm,rect_sm] = imcrop(sm);close(gcf)  
% end

% sm=sm';
% 
% spotpos=fliplr(spotpos);
% sm=imrotate(sm,-90);
% 
% spotpos=[length(sm) 0]+[-1 1].*fliplr(spotpos);

disp([' processing :',namebase,' -- ',fluotxt])
disp([smf,'   --  ',subtf]);

xsearchlim=min([300,floor(spotpos(2)-size(subtomo,2)/2),size(sm,2)-floor(spotpos(2)+size(subtomo,2)/2)])-1;
ysearchlim=min([300,floor(spotpos(1)-size(subtomo,1)/2),size(sm,1)-floor(spotpos(1)+size(subtomo,1)/2)])-1;

searchlim=50;
if xsearchlim<0
   xsearchlim=searchlim;
end 
if ysearchlim<0
    ysearchlim=searchlim;
end

sub_sm=sm(floor(spotpos(2)-size(subtomo,2)/2-xsearchlim):floor(spotpos(2)+size(subtomo,2)/2+xsearchlim),floor(spotpos(1)-size(subtomo,1)/2)-ysearchlim:floor(spotpos(1)+size(subtomo,1)/2+ysearchlim));
rect2_sm1=size(subtomo);
rect_sm=[floor(spotpos(1)-size(subtomo,1)/2-ysearchlim) floor(spotpos(2)-size(subtomo,2)/2)-xsearchlim];
rect_sm(4)=rect2_sm1(2)+xsearchlim;
rect_sm(3)=rect2_sm1(1)+ysearchlim;

disp(['xsl: ',num2str(xsearchlim),'   --   ysl: ',num2str(ysearchlim)]);



% searchlim=170;
% 
% xsearchlim=searchlim;
% ysearchlim=searchlim;
% 
% sub_sm=sm(floor(spotpos(2)-size(subtomo,2)/2-xsearchlim):floor(spotpos(2)+size(subtomo,2)/2+xsearchlim),floor(spotpos(1)-size(subtomo,1)/2)-ysearchlim:floor(spotpos(1)+size(subtomo,1)/2+ysearchlim));
% rect2_sm1=size(subtomo);
% rect_sm=[floor(spotpos(1)-size(subtomo,1)/2-ysearchlim) floor(spotpos(2)-size(subtomo,2)/2)-xsearchlim];
% rect_sm(4)=rect2_sm1(2)+xsearchlim;
% rect_sm(3)=rect2_sm1(1)+ysearchlim;



   
    i=1;
    for subslice=slicerange
         stslice{i}=imrotate(subtomo(:,:,subslice),90);
         slice1(i)=subslice;
         xcorr(i).c=normxcorr2(imadjust(stslice{i}),imadjust(sub_sm));
    %    figure, surf(c), shading flat
% pause;    imshow(stslice{i})
        % offset found by correlation
        [max_c(i), imax] = max(abs(xcorr(i).c(:)));
        [xpeak(i), ypeak(i)] = ind2sub(size(xcorr(i).c),imax(1));
        i=i+1;
    end  
    
    



    
    [totmax, itotmax]=max(max_c);
% end

disp('------------------------------')
disp('optimal correlation coefficient:')
disp(num2str(totmax))

if totmax<0.9

%     searchlim = 100;    
%     clear slice1 stslice 
% 
%     sub_sm=sm(floor(spotpos(2)-size(subtomo,2)/2-searchlim):floor(spotpos(2)+size(subtomo,2)/2+searchlim),floor(spotpos(1)-size(subtomo,1)/2)-searchlim:floor(spotpos(1)+size(subtomo,1)/2+searchlim));
%     rect2_sm1=size(subtomo);
%     rect_sm=[floor(spotpos(1)-size(subtomo,1)/2-searchlim) floor(spotpos(2)-size(subtomo,2)/2)-searchlim];
%     rect_sm(4)=rect2_sm1(2)+2*searchlim;
%     rect_sm(3)=rect2_sm1(1)+2*searchlim;
% 
% 
%     i=1;
%     for subslice=1:size(subtomo,3)
%          stslice{i}=imrotate(subtomo(:,:,subslice),90);
%          slice1(i)=subslice;
%          xcorr(i).c=normxcorr2(stslice{i},sub_sm);
%     %    figure, surf(c), shading flat
%     
%         % offset found by correlation
%         [max_c(i), imax] = max(abs(xcorr(i).c(:)));
%         [xpeak(i), ypeak(i)] = ind2sub(size(xcorr(i).c),imax(1));
%         i=i+1;
%     end  
%     [totmax, itotmax]=max(max_c);
% % end
disp('------------------------------')
disp('optimal correlation coefficient:')
disp(num2str(totmax)) 
    
    
end

% offset of position of subimage
xoffset = ypeak(itotmax)-size(subtomo,1)+rect_sm(1);
yoffset = xpeak(itotmax)-size(subtomo,2)+rect_sm(2);
xbegin = xoffset+1;
xend   = xoffset+ size(subtomo,1);
ybegin = yoffset+1;
yend   = yoffset+size(subtomo,2);

slice=slice1(itotmax);

b=sm;
b(spotpos(2)-5:spotpos(2)+5,spotpos(1)-5:spotpos(1)+5,3)=255;
b(ybegin:yend,xbegin:xend,2)=imadjust(stslice{itotmax});
b(1,1,3)=0;

imshow(b);




subpos=round(spotpos-[xoffset-1 yoffset-1]);

% subpos2=[length(sm) 0]+[-1 1].*fliplr(spotpos);


subpos2=[0 size(subtomo,2)+1]+[1 -1].*subpos;


if min([subpos subpos2])>5 
figure
c=imadjust(stslice{itotmax});
c(subpos(2),subpos(1),2)=255;
c(subpos(2)-5:subpos(2)+5,subpos(1)-5:subpos(1)+5,3)=255;
imshow(c)

figure
d=imrotate(c(:,:,2),270);
f=imrotate(c(:,:,1),270);
f(subpos2(1),subpos2(2),2)=255;
f(subpos2(1)-5:subpos2(1)+5,subpos2(2)-5:subpos2(2)+5,3)=255;
imshow(f)

% k=find(d>3);
% [subpos2(1) subpos2(2)] = ind2sub([size(subtomo,1) size(subtomo,2)],k);


end

subpos2=subpos2+[-1 1]
close all

smf_sub=smf;subtf;
save([outfileroot,namebase,fluotxt,'.subttfm.mat'],'subpos2','slice','smf_sub','subtf');





% 
% 
% switch fluotxt
%     
%     case '_GFP'
%         db{dbindex+1,21} = subpos2(1);
%         db{dbindex+1,22} = subpos2(2);
%     case '_RFP'
%         db{dbindex+1,27} = subpos2(1);
%         db{dbindex+1,28} = subpos2(2);
% end



file_2 = fopen([outfileroot,namebase,fluotxt,'_subtomo_transform.log'],'w');
fprintf(file_2,[outfileroot,namebase,fluotxt,'_subtomo_transform.log      ---   Logfile of Subtomogram-transformation\n\n']);
fprintf(file_2,[' slice of interest: ',smf,' --  subtomogram: ',subtf,'\n\n slice: ',int2str(slice)]);
fprintf(file_2,'   -  coordinates of transformed fluorescence spot:');
fprintf(file_2,'%2.3f %2.3f',subpos2);
fclose(file_2);
end

a=who;
excl=find(ismember(a,[keep;'a']));
a(excl)=[];
for ii=1:length(a)
clear(a{ii});
end

end

% cell2csv('Endocytosis_1.csv',db,';',1);
