function [C,v,eval] = mat_projection(W)
n1 = size(W,1);
n2 = size(W,2);
[s,v,d] = svd(full(W));
C = s*eye(n1,n2)*d';
v = diag(v);
eval = [range(v),mean(v),var(v)];
end