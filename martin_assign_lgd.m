function [T_out,B_out,IX_out,a_out] = martin_assign_lgd(f,e,n)
%% uses local geometric descriptors to assign 2D points in f to those in e
% Input, the two 2-vectors and the number of neighbors to consider
% Output: et -> the transformed vectors proposed to match f
%         tr  -> the transformation that takes e->et (as defined in
%         matlab's procrustes documentation
%         Ix  -> a 2-vector of indices proposed to match
%
verbose = 0;
if n<2, error('n must be larger than one');end
%% fluorescence image fiducials
fd =zeros(size(f,1), size(f,1));    % distance matrix
for ix = 1:size(f,1)
    for jx = 1:size(f,1)
        fd(ix,jx) = sqrt((f(ix,1)-f(jx,1)).^2 +(f(ix,2)-f(jx,2)).^2);
    end
end
% generate the lists of immediate neighbors n and store their coordinates in fn (this is the
% local constellation around every point)
fn = zeros(size(f,1),n,2);   % list of coordinates of the n neighbors
[S B] = sort(fd);   % S is the sorted distance matrix, B the indices into rows of f
for ix = 1:size(f,1)
    for in = 2:1+n
        fn(ix,in-1,:) = [f(B(in,ix),1)-f(ix,1) f(B(in,ix),2)-f(ix,2)]; % note that we subtract the coordinates of point ix to center the constellation around it
    end
end

%% do the same for em case
ed =zeros(size(e,1), size(e,1));
for ix = 1:size(e,1)
    for jx = 1:size(e,1)
        ed(ix,jx) = sqrt((e(ix,1)-e(jx,1)).^2 +(e(ix,2)-e(jx,2)).^2);
    end
end
% generate the lists of immediate neighbors n and store their coordinates in fn (this is the
% local constellation around every point)
en = zeros(size(e,1),n,2);   % list of coordinates of the n neighbors
[S B] = sort(ed);   % S is the sorted distance matrix, B the indices into rows of f
for ix = 1:size(e,1)
    for in = 2:1+n
        en(ix,in-1,:) = [e(B(in,ix),1)-e(ix,1) e(B(in,ix),2)-e(ix,2)]; % note that we subtract the coordinates of point ix to center the constellation around it
    end
end
%% with fn and en we can now calculate the procrustes correspondence of constellations around every point in f to that in e
dmx = zeros(size(f,1), size(e,1));  % store the indices of the pair and the dissimilarity measure
bmx = zeros(size(f,1), size(e,1));  % store the indices of the pair and the dissimilarity measure

T = {};
counter = 0;
for ix = 1:size(f,1),       % loop over the rows of f
    Xf = reshape(fn(ix,:,:), n,2);
    for jx = 1:size(e,1),   % now loop over rows of e
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
err = zeros(size(f,1),1);
for ix = 1:size(f,1),
    if a(ix)
        err(ix) = dmx(ix,a(ix));
    end
end
Ix = [[1:size(f,1)]' a(:) err(:)];       % output
if verbose, disp([[1:size(f,1)]' a(:) err(:)]);end

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
B_out = B(IX_good);

IX_out=IX(IX_good);

% identify good transformation
%disp([IX(1) a(IX(1))]);
%disp(T{IX(1), a(IX(1))}.T);

for i=1:length(IX_out)
    T_out{i}=T{IX_out(i), a(IX_out(i))};
end

a_out=a';


% Z = b*Xe*D + c , i.e., z = scale*X2*rotation + translation;
% % perform this "good" transformation on all em fiducials (also for possible inspection)
% et = [];
% for ix = 1:size(e,1)
%     z = tr.b*e(ix,:)*tr.T + f(IX(1),:);
%     et = [et;z];
% end
% 
% % % perform a translation to match the point clouds based on the good
% % % transformation information
% correction = et(a(IX(1)),:)-f(IX(1),:);
% for ix = 1:size(et,1)
%          %et(ix,:) = et(ix,:)-(et(a(ix),:)-f(ix,:)); %  translate according to good transformation
%         et(ix,:) = et(ix,:)-correction;
% end