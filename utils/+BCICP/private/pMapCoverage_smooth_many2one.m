% T1: S2 -> S1; X2 -> X1
% T2: S1 -> S2; X1 -> X2
% deal with the case: many-to-one, smooth it out:
% assume v_i on S2 are all mapped to w on S1
% remap v_i to the neighbors of w
function T21_new = pMapCoverage_smooth_many2one(T21,S2,S1)
X1 = S1.surface.VERT; X2 = S2.surface.VERT;
n1 = S1.nv; n2 = S2.nv;  % nubmer of vertices

freq = get_pMap_range(T21,n1);
% the vertices that are not covered by the pMap T1; on shape S1
% nc_id = reshape(setdiff(1:n1,unique(T1)),[],1);
nc_id = (1:n1)';

% the vertices on S1 that have a lot of preimages - hurts the bijectivity
vid_all = find(freq > 4);
% for each vtx in vid_all, find the vertices on S2 that are mapped to v_id
src_vid_all = cellfun(@(v_id) find(T21 == v_id),num2cell(vid_all),'UniformOutput',false);
displacement_all = cellfun(@(v_id, src_v_id) X1(v_id,:) - mean(X2(src_v_id,:)),...
    num2cell(vid_all),src_vid_all,'UniformOutput',false); % find the displacement vector

% find the new target vertices on S1 for those vertices in src_vid_all
new_target_all = cellfun(@(src_vid,displacement) nc_id(knnsearch(X1(nc_id,:),X2(src_vid,:) + displacement)),...
    src_vid_all,displacement_all,'UniformOutput',false);


source_id = cell2mat(src_vid_all(:)); % vtx on Shape Ss
target_old = T21(source_id); % vtx on shape S1
target_new = cell2mat(new_target_all(:)); % vertex on shape S1

% compare if the new targets are better
% get_pair_dist = @(X,id1,id2) sum((X(id1,:)- X(id2,:)).^2,2);
% err_old = get_pair_dist(X2,T2(target_old),source_id);
% err_new = get_pair_dist(X2,T2(target_new),source_id);

% id = find(err_new <  err_old);
T21_new = T21;
% T1_out(source_id(id)) = target_new(id);
T21_new(source_id) = target_new;
end
