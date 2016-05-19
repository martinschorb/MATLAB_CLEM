function [T_out,err_out,n_out,f_shift,correction] = martin_assign_lgd(f,e,n_min,n_max)
%% uses local geometric descriptors to assign 2D points in f to those in e
% Input, the two 2-vectors and the number of neighbors to consider
% Output: et -> the transformed vectors proposed to match f
%         tr  -> the transformation that takes e->et (as defined in
%         matlab's procrustes documentation
%         Ix  -> a 2-vector of indices proposed to match
%
verbose = 0;
if n_min<2, error('n must be larger than one');end
%% fluorescence image fiducials

szf = size(f,1);

fd =zeros(szf, szf);    % distance matrix
for ix = 1:szf
    for jx = 1:szf
        fd(ix,jx) = sqrt((f(ix,1)-f(jx,1)).^2 +(f(ix,2)-f(jx,2)).^2);
    end
end


[S B_f] = sort(fd);
% S is the sorted distance matrix, B the indices into rows of f

sze = size(e,1);

ed =zeros(sze, sze);    
for ix = 1:sze
    for jx = 1:sze
        ed(ix,jx) = sqrt((e(ix,1)-e(jx,1)).^2 +(e(ix,2)-e(jx,2)).^2);
    end
end
[S B_e] = sort(ed); % S is the sorted distance matrix, B the indices into rows of e


for n = n_min:n_max
    
    i_n = n-n_min+1;
    
    fn = zeros(szf,n,2);   % list of coordinates of the n neighbors
    for ix = 1:szf
        for in = 2:1+n
%             fn(ix,in-1,:) = [f(B_f(in,ix),1)-f(ix,1) f(B_f(in,ix),2)-f(ix,2)]; % note that we subtract the coordinates of point ix to center the constellation around it
            fn(ix,in-1,:) = f(B_f(in,ix),:);
        end
    end

   
    % generate the lists of immediate neighbors n and store their coordinates in fn (this is the
    % local constellation around every point)
    en = zeros(sze,n,2);   % list of coordinates of the n neighbors
    for ix = 1:sze
        for in = 2:1+n
            en(ix,in-1,:) = [e(B_e(in,ix),1)-e(ix,1) e(B_e(in,ix),2)-e(ix,2)]; % note that we subtract the coordinates of point ix to center the constellation around it
            en(ix,in-1,:) = e(B_e(in,ix),:);        
        end
    end
    %% with fn and en we can now calculate the procrustes correspondence of constellations around every point in f to that in e
    dmx = zeros(szf, sze);  % store the indices of the pair and the dissimilarity measure
    bmx = zeros(szf, sze);  % store the indices of the pair and the dissimilarity measure

    T = {};
    counter = 0;
    for ix = 1:szf,       % loop over the rows of f
        Xf = reshape(fn(ix,:,:), n,2);
        for jx = 1:sze,   % now loop over rows of e
            counter = counter + 1;
            Xe = reshape(en(jx,:,:), n,2);
            %% calculate the dissimilarity using procrustes analysis
            [d z tr] = procrustes(Xf,Xe,'reflection',false, 'scaling', true);       % requires the statistics toolbox
                                                % d residual, z transformed
                                                % coordinates, tr details of
                                                % transformation
            dmx(ix,jx)= d;         % store the residuals
            bmx(ix,jx) = tr.T(1);
            T{ix,jx} = tr;
        end
    end
    % dmx(:) = dmx(:)/max(dmx(:));      % normalize the measure of post-transformation deviation

    %% %% find the unique combinations such that each bead is used only once
    % in combination with its best match: This is a classical "assignment
    % problem" that can be solved using the Hungarian algorithm. 
    % The problem is to match the constellations such that we minimize the total cost.
    [a,c] = munkres(dmx);       % execute the Hungarian algorithm
    % "a" contains the vector of indices matching the rows of dmx (in our case
    % the fluorescence data)
    % Alternatively use: [a,c] = Hungarian(dmx);[a c] = ind2sub(size(a),find(a));disp([c(:) a(:)]);% gives same result as above

    %% We are not necessarily done yet!! 
    % Due to possible differences in the number of fiducials
    % and because we have typically very few (less than 30) fiducials we have
    % missing matches which cannot be dealt with at the level of
    % the matching algorithm due to inherent inaccuracy in the cost matrix

    % To make that point let us look at the mismatch error associated with the assignments
    err = zeros(szf,1);
    for ix = 1:szf,
        if a(ix)
            err(ix) = dmx(ix,a(ix));
        end
    end
    Ix = [[1:szf]' a(:) err(:)];       % output
    if verbose, disp([[1:szf]' a(:) err(:)]);end

    % The strategy is to:
    % [1] identify the "good" transformations (this usually involves the
    % constellations that do NOT include the missing matches -- we do not know
    % which beforehand)
    % [2] transform the full set
    % [3] perform a new match based on simple distances (after transformation)



    [B, IX] = sort(err);
    IX(B==0) = [];
    B(B==0) = [];

    B1 = otsu(B,2);
    IX_good = find(B1-2);
    err_out{i_n} = B(IX_good);

    IX_out{i_n} = IX(IX_good);

    % identify good transformation
    %disp([IX(1) a(IX(1))]);
    %disp(T{IX(1), a(IX(1))}.T);

    for i=1:length(IX_out{i_n})
        tr = T{IX_out{i_n}(i), a(IX_out{i_n}(i))};
%         et=[];
        T_out{i_n}{i}=tr;
        
        et = [];
        for ix = 1:sze
        %     z = tr.b*e(ix,:)*tr.T + f(IX(1),:);
            f_shift{i_n}{i} = f(IX_out{i_n}(i),:);
            z = tr.b*e(ix,:)*tr.T + f(IX_out{i_n}(i),:);
            et = [et;z];
        end

        % % perform a translation to match the point clouds based on the good
        % % transformation information
        % correction = et(a(IX(1)),:)-f(IX(1),:);
        correction{i_n}{i} = et(a(IX_out{i_n}(i)),:) - f(IX_out{i_n}(i),:);
        for ix = 1:size(et,1)
                 %et(ix,:) = et(ix,:)-(et(a(ix),:)-f(ix,:)); %  translate according to good transformation
                et(ix,:) = et(ix,:)-correction{i_n}{i};
        end

         
    end

    a_out(:,i_n)=a';
    
    
    

end

n_out = n_min:n_max;
% 
% % Z = b*Xe*D + c , i.e., z = scale*X2*rotation + translation;
% % % perform this "good" transformation on all em fiducials (also for possible inspection)
% et = [];
% for ix = 1:sze
% %     z = tr.b*e(ix,:)*tr.T + f(IX(1),:);
%     z = tr.b*e(ix,:)*tr.T + f(IX_out{i_n}(i),:);
%     et = [et;z];
% end
% 
% % % perform a translation to match the point clouds based on the good
% % % transformation information
% % correction = et(a(IX(1)),:)-f(IX(1),:);
% correction = et(a(IX_out{i_n}(i)),:) - f(IX_out{i_n}(i),:);
% for ix = 1:size(et,1)
%          %et(ix,:) = et(ix,:)-(et(a(ix),:)-f(ix,:)); %  translate according to good transformation
%         et(ix,:) = et(ix,:)-correction;
% end