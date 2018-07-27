function T21_new = pMapCoverage6(T21,S2,S1,T12_new)
n1 = S1.nv; n2 = S2.nv;
dist2 = compute_mesh_dist_matrix(S2);
get_mat_entry = @(M,I,J) M(sub2ind(size(M),I,J));

T21_new = T21;

s1_ncid = reshape(setdiff(1:n1,T21_new),[],1); % vertices on S1 that are not covered
s2_ncid = reshape(setdiff(1:n2,T12_new),[],1); % vertices on S2 that are not covered
freq = get_pMap_range(T21_new,n1);
% vertices on Shape S2 that are mapped to a vertex with multiple preimages
src_id = find(ismember(T21_new,find(freq > 1)));
[~,ic] = min(dist2(T12_new(s1_ncid), src_id),[],2);
src1 = src_id(ic);
tar1 = s1_ncid; % new target - uncovered
tar2 = T21_new(src1); % covered

new_pair = [src1,tar1];
new_pair = unique(new_pair,'rows');
if ~isempty(new_pair)
    get_pair_dist = @(X,id1,id2) sum((X(id1,:) - X(id2,:)).^2,2);
    err_new = get_mat_entry(dist2,T12_new(new_pair(:,2)),new_pair(:,1));
    err_old = get_mat_entry(dist2,T12_new(T21_new(new_pair(:,1))),new_pair(:,1));
    
    id = err_new < err_old;
    src_id = unique(new_pair(id,1));
    tmp1 = cellfun(@(x) new_pair(new_pair(:,1) == x,2),num2cell(src_id),'UniformOutput',false);
    tmp2 = cellfun(@(x) err_new(new_pair(:,1)== x),num2cell(src_id),'UniformOutput',false);
    tmp3 = cellfun(@(x) find(x == min(x),1),tmp2);
    tar_id = cellfun(@(x,id) x(id),tmp1,num2cell(tmp3));
    T21_new(src_id) = tar_id;
end
end
