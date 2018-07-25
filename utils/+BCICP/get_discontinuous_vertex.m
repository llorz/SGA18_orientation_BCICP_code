function [v_id, MAX_EDGE_LENGTH] = get_discontinuous_vertex(T21,S2,S1)
get_mat_entry = @(M,I,J) M(sub2ind(size(M),I,J));
dist1 = compute_mesh_dist_matrix(S1);
dist2 = compute_mesh_dist_matrix(S2);

d2 = get_mat_entry(dist2,S2.Elist(:,1),S2.Elist(:,2)); % edge length on S2
d1 = get_mat_entry(dist1,S1.Elist(:,1),S1.Elist(:,2)); % edge length on S1
MAX_EDGE_LENGTH = max(max(d2),max(d1));

d2_mapped = get_mat_entry(dist1, T21(S2.Elist(:,1)), T21(S2.Elist(:,2)));
% distortion = d2_mapped./d2; % edge distortion
e_id = d2_mapped > MAX_EDGE_LENGTH; % edges with large error
v_id = unique(S2.Elist(e_id,:)); % vertex with large error
v_id = reshape(v_id,[],1);
end