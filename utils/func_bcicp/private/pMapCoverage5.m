% slow!
function T21_new = pMapCoverage5(T21,S2,S1,T12_new)
X1 = S2.surface.VERT; X2 = S2.surface.VERT;
N1 = S1.vtx_neigh; N2 = S2.vtx_neigh;
n1 = S1.nv; n2 = S2.nv;

dist2 = compute_mesh_dist_matrix(S2);
get_mat_entry = @(M,I,J) M(sub2ind(size(M),I,J));


T21_new = T21; 
freq = get_pMap_range(T21_new,n1);

s1_ncid = reshape(setdiff(1:n1,T21_new),[],1); % vertices on S1 that are not covered
s2_ncid = reshape(setdiff(1:n2,T12_new),[],1); % vertices on S2 that are not covered

S2_id = 1:n2;
S2_neigh = cellfun(@(neigh) neigh(freq(T21_new(neigh))>1),N2(S2_id),'UniformOutput',false);
S1_neigh = cellfun(@(neigh) neigh(ismember(neigh,s1_ncid)),N1(T21_new(S2_id)),'UniformOutput',false);

get_length = @(X) cellfun(@(x) length(x),X);
change_id = intersect(find(get_length(S2_neigh)),find(get_length(S1_neigh)));

new_pair = cell2mat(cellfun(@(neigh1,neigh2)[neigh2, neigh1(knnsearch(X1(neigh1,:),X2(neigh2,:)))],...
    S1_neigh(change_id),S2_neigh(change_id),'UniformOutput',false));

new_pair = unique(new_pair,'rows');

err_new = get_mat_entry(dist2,T12_new(new_pair(:,2)),new_pair(:,1));
err_old = get_mat_entry(dist2,T12_new(T21_new(new_pair(:,1))),new_pair(:,1));

id = find(err_new < err_old);
tmp1 = cellfun(@(x) new_pair(new_pair(:,1) == x,2),num2cell(unique(new_pair(id,1))),'UniformOutput',false);
tmp2 = cellfun(@(x) err_new(new_pair(:,1)== x),num2cell(unique(new_pair(id,1))),'UniformOutput',false);
tmp3 = cellfun(@(x) find(x == min(x),1),tmp2);
tar_id = cellfun(@(x,id) x(id),tmp1,num2cell(tmp3));
src_id = unique(new_pair(id,1));

T21_new(src_id) = tar_id;

end

