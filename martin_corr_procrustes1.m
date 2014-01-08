tic

beadchoiceindex = 1; % need to sort neighboring beads according to distance, otherwise complexity explodes

pxs_fm=64.5;
pxs_em=2.53;


numbad = 3;         % number of potential bad beads to omit.


dist_thr = 10;      % distance threshold in nm while selecting "twin" beads based on distance (should be < accuracy)

scale0 = pxs_fm/pxs_em;

scale_thr = 0.005;

num_trafo = 20;



%% ----------------------------------------------------------


n=min(size(e,1),size(f,1));

data_all=[];


e0 = e - repmat(e(beadchoiceindex,:),[size(e,1) 1]);


[edist, e_idx] = sort(sqrt(sum(e0.^2,2)));

e_twins = find(abs(diff(edist))*pxs_em < dist_thr);

e_ds = e(e_idx,:);

f0 = f - repmat(e(beadchoiceindex,:),[size(e,1) 1]);

[fdist, f_idx] = sort(sqrt(sum(f0.^2,2)));

f_twins = find(abs(diff(fdist))*pxs_fm < dist_thr);

if or(length(e_twins)>1,length(f_twins)>1)
    warning('Possible mismatch in bead pairs! Try to reduce distance threshold.');
end

% sort bead coordinates

f_ds = f(f_idx,:);


%%

for i_bad=0:numbad

    e_sel = martin_combin(size(e,1),n-i_bad);

    % add twins EM
    e_sel1 = e_sel;
    [e_tw e_twx1] = ismember(e_sel,[e_twins e_twins+1]);
    e_twx2 = sum(e_twx1,2)~=3;
    e_tw(e_twx2,:) = [];
    e_twx1(e_twx2,:) = [];
    e_sel1(e_twx2,:) = [];
    e_sel1 = e_sel1+(1-2*e_twx1+2).*e_tw;
    e_sel = [e_sel;e_sel1];


    f_sel=martin_combin(size(f,1),n-i_bad);

    % add twins FM
    f_sel1 = f_sel;
    [f_tw f_twx1] = ismember(f_sel,[f_twins f_twins+1]);
    f_twx2 = sum(f_twx1,2)~=3;
    f_tw(f_twx2,:) = [];
    f_twx1(f_twx2,:) = [];
    f_sel1(f_twx2,:) = [];
    f_sel1 = f_sel1+(1-2*f_twx1+2).*f_tw;
    f_sel = [f_sel;f_sel1];



    e_size{i_bad+1} = size(e_sel,1);
    f_size{i_bad+1} = size(f_sel,1);

    for ix=1:e_size{i_bad+1};
        Xe = e_ds(e_sel(ix,:)',:);
        for jx=1:f_size{i_bad+1}
            Xf = f_ds(f_sel(jx,:)',:);

             [d z tr] = procrustes(Xe,Xf,'reflection',false, 'scaling', true);       % requires the statistics toolbox
                                                % d residual, z transformed
                                                % coordinates, tr details of
                                                % transformation

            dist = sqrt(sum(sum((z-Xe).^2))) / size(Xe,1);           % use geometric distances instead of procrustes because otherwise we cannot compare different numbers of "bad" beads                         

            dmx{i_bad+1}(ix,jx)= dist;            % store the residuals
            scale{i_bad+1}(ix,jx) = tr.b;         % store the scales
            rel_scale{i_bad+1}(ix,jx) = tr.b / scale0;

            T{i_bad+1}{ix,jx} = tr;               % store the transformations
        end  
    end

    [d_sorted{i_bad+1},idx_sorted{i_bad+1}] = sort(dmx{i_bad+1}(:)/(n-i_bad));
    data_all=[data_all;idx_sorted{i_bad+1} repmat(i_bad,size(d_sorted{i_bad+1})) d_sorted{i_bad+1} reshape(abs(rel_scale{i_bad+1}(idx_sorted{i_bad+1})-1),size(d_sorted{i_bad+1}))];


end

scale_score = cellfun(@(x) abs(x-1),rel_scale,'UniformOutput',0);


data_clean=sortrows(data_all,3);

data_clean(data_clean(:,4)>scale_thr,:)=[];

%%

for i_good=1:num_trafo
    
    g_idx = data_clean(i_good,2);
    
    [aa,bb]=ind2sub([e_size{g_idx+1},f_size{g_idx+1}],data_clean(i_good,1));
    
    good_T(i_good) = T{g_idx+1}{aa,bb};
    
    good_z(:,:,i_good) = good_T(i_good).b*f*good_T(i_good).T+repmat(good_T(i_good).c(1,:),[size(f,1) 1]);                 %   b*Y*T+c.
       
end


toc