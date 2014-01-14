tic

pxs_fm=139;
pxs_em=3.8936;

max_n = 12;


numbad = 2;         % number of potential bad beads to omit.


dist_thr = 50;      % distance threshold in nm while selecting "twin" beads based on distance (should be < accuracy)

d_thr_match = 100;  % distance threshold in nm for possible bead combinations

scale0 = pxs_fm/pxs_em;

scale_thr = 0.005;

num_trafo = 20;



%% ----------------------------------------------------------

sze = size(e,1);
szf = size(f,1);

minsize = min([sze,szf,max_n]);

d={};

[edist edix]=sort(squareform(pdist(e*pxs_em)));


edist=edist';
edix=edix';
 % add twins EM
e_twins = find(abs(diff(edist,1,2)) < dist_thr);
[et_row,et_col]=ind2sub(size(edist)-[0 1],e_twins);

for i=1:length(et_row)
    edix(end+1,:)=edix(et_row(i),:);
    edix(end,et_col(i):et_col(i)+1)=fliplr(edix(end,et_col(i):et_col(i)+1));
end

[fdist fdix] = sort(squareform(pdist(f*pxs_fm)));
fdist = fdist';
fdix = fdix';
%%

%  for i_start = 1:sze
%     e_check = round(edist(1,2:end)/d_thr_match);
%     e_check = [e_check,e_check+1,e_check-1];
%         
%     f_skip = [];
%     
%     f_sel = 1:szf;
%     for i_f_check = 1:szf
%         f_check= round(fdist(i_f_check,2:end)/100);
%         if sum(ismember(e_check,f_check))<4
%             f_skip = [f_skip,i_f_check];
%         end
%     end
% 
%     f_sel(f_skip)=[];
%     
%     
%  end    
%     twins = find(edix(:,1)==i_start);
    
      
%     for i_twins = 1:length(twins)
        
        e_in = e;%(edix(twins(i_twins),:),:);
    

%         for i_f = 1:length(f_sel)
            f_in = f;%(fdix(f_sel,:),:);
            
             
            for nums = 1:numbad+1            
                n(nums) = nums - 1;
                [T{nums},err(:,nums),IX(:,nums),a(:,nums)] = martin_assign_lgd(e,f,minsize-nums);
            end
            
            
            
            
%         end
%     end
% end

toc






toc