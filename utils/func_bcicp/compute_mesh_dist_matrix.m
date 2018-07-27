function d = compute_mesh_dist_matrix(S)
    if isfield(S,'Gamma') % exists pre-computed geodesic distance matrix
        d = S.Gamma;
    else % compute the Euclidean distance matrix
        X = S.surface.VERT;
        d = dist_xy(X,X);
    end
end