% x: matrix of n points, each row represents one point
% y: matrix of m points, each row represents one point
% dist(i,j) = norm(x(i,:) - y(j,:)), a n by m matrix
function dist = dist_xy(x,y)
switch nargin
    case 2
        n = size(x,1);
        m = size(y,1);
        dist = diag(x*x')*ones(1,m) - 2*x*y' + ones(n,1)*diag(y*y')';
        dist(dist < 1.0e-12) = 0;   % error comes from the computation?
        dist = sqrt(dist);
    case 1  % dist_xy(x,x)
        dist = squareform(pdist(x));
    otherwise
        warning('Input must be two/one matrix of observations')
end
end