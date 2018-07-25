% T is mapped to a shape S2 with n2 vertices (i.e T: S1 -> S2)
% f(i): how many vertices in S1 are mapped to the i-th vtx on shape S2
function [f,coverage] = get_pMap_range(T,n2)
    f = zeros(n2,1);
    f(unique(T)) = histc(T,unique(T));
    coverage = length(unique(T))/n2; % coverage rate
end