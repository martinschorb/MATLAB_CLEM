tic

pxs_fm=139;
pxs_em=3.8936;

imsz = [4096,4096];

n_max = 12;
n_min = 4;

numbad = 2;         % number of potential bad beads to omit.


dist_thr = 50;      % distance threshold in nm while selecting "twin" beads based on distance (should be < accuracy)

d_thr_match = 100;  % distance threshold in nm for possible bead combinations

scale0 = pxs_fm/pxs_em;

scale_thr = 0.03;

num_trafo = 20;



%% ----------------------------------------------------------

sze = size(e,1);
szf = size(f,1);

minsize = min([sze,szf,n_max]);

% d={};
% 
% [edist edix]=sort(squareform(pdist(e*pxs_em)));
% 
% 
% edist=edist';
% edix=edix';
%  % add twins EM
% e_twins = find(abs(diff(edist,1,2)) < dist_thr);
% [et_row,et_col]=ind2sub(size(edist)-[0 1],e_twins);
% 
% for i=1:length(et_row)
%     edix(end+1,:)=edix(et_row(i),:);
%     edix(end,et_col(i):et_col(i)+1)=fliplr(edix(end,et_col(i):et_col(i)+1));
% end
% 
% [fdist fdix] = sort(squareform(pdist(f*pxs_fm)));
% fdist = fdist';
% fdix = fdix';
%%

corrs = [];

% corrs contains the necessary transformation parameters for each matched bead
% position:corrs = [xpos ypos scale rotation_vector_x rv_y x_trans y_trans local prediction error];


    
    e_in = e;
%     e_in(i_check,:)=[];
    
    [T_good,err_good,trans_good] = martin_tfselect(e_in,f,n_min,minsize,scale0,scale_thr);
% 
    if ~isempty(T_good)
 
        for i_check = 1:sze       
        pred_diff =[];   
        for i=1:length(err_good)
            z{i_check}{i}=T_good{i}.b*f*T_good{i}.T+repmat(trans_good{i},[szf,1]);        
            dists{i_check}{i} = sort(sqrt(sum((z{i_check}{i}-repmat(e(i_check,:),[szf,1])).^2,2)));
            pred_diff(i)= dists{i_check}{i}(1);        
        end

        % clean for misassigned beads

        edists = round(sort(sqrt(sum((e-repmat(e(4,:),[sze,1])).^2,2)))/100);

        cleanidx=[];

        for j = 1:length(err_good)
            cdists = round(dists{i_check}{j}/100);
            if sum(ismember(cdists,edists)) < 4
                cleanidx=[cleanidx,j];
            end
        end

        pred_diff(cleanidx)=2*max(pred_diff);
    %     

        [pd,pdi]=sort(pred_diff*pxs_em);



        if pd > 200
           % this bead/feature in EM seems to have no correcponding detected
           % fluorescence signal position
           disp(['mismatch EM bead # ',num2str(i_check),' -> pd = ',num2str(pd),' ... skipping']);
        else


%            av_select = otsu(pd,2);
           sel_idx = find(otsu(pd,4)==1);

           for k = sel_idx
               scales(k) = T_good{pdi(k)}.b;
               rotation(k,:) = T_good{pdi(k)}.T(1:2);
               trans(k,:) = trans_good{pdi(k)};
               z_sel{i_check}(:,:,k) = z{i_check}{pdi(k)};
               
           end

           av_scale{i_check} = mean(scales);
           av_rot{i_check} = mean(rotation);
           av_trans{i_check} = mean(trans);
           
           
           R{i_check} = [av_rot{i_check}(1) -av_rot{i_check}(2);fliplr(av_rot{i_check})];
           
           z_av(:,:,i_check) = av_scale{i_check}* f *R{i_check}/sqrt(det(R{i_check})) + repmat(av_trans{i_check},[szf,1]);
           av_error{i_check} = min(sqrt(sum((z_av(:,:,i_check)-repmat(e(i_check,:),[szf,1])).^2,2)));

           corrs = [corrs; e(i_check,:) av_scale{i_check} av_rot{i_check} av_trans{i_check} av_error{i_check}];

        end
    end
    
end            
      

[X,Y]=meshgrid(linspace(1,imsz(1)),linspace(1,imsz(2)));
% generates a 100x100 mesh covering the em image

gridpos=[1 1;1 imsz(2);imsz(1) 1;imsz;corrs(:,1:2)]; % adds corners

scales = [repmat(scale0,[4,1]);corrs(:,3)];
scale2D = griddata(gridpos(:,1),gridpos(:,2),scales,X,Y,'cubic');

rot_grid = [repmat(mean([corrs(:,4) corrs(:,5)]),[4,1]); corrs(:,4:5)];
rot2D_x = griddata(gridpos(:,1),gridpos(:,2),rot_grid(:,1),X,Y,'cubic');
rot2D_y = griddata(gridpos(:,1),gridpos(:,2),rot_grid(:,2),X,Y,'cubic');

transl_grid = [repmat(mean([corrs(:,6) corrs(:,7)]),[4,1]); corrs(:,6:7)];
transl2D_x = griddata(gridpos(:,1),gridpos(:,2),transl_grid(:,1),X,Y,'cubic');
transl2D_y = griddata(gridpos(:,1),gridpos(:,2),transl_grid(:,2),X,Y,'cubic');

err_grid = [repmat(max(corrs(:,8)),[4,1]); corrs(:,8)];
err2D = griddata(gridpos(:,1),gridpos(:,2),err_grid,X,Y,'cubic');

toc

