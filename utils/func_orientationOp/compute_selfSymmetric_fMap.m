function [C12] = compute_selfSymmetric_fMap(S1,S2,B1,B2,Ev1,Ev2,fct_src,fct_tar,para)
a = 1; % Descriptors preservation
b = 1; % Commutativity with descriptors
c = 1; % Commutativity with Laplacian 
d = 1; % orientation-reversing/preserving term
numEigsSrc = size(B1,2); numEigsTar = size(B2,2);
if nargin > 8
    if isfield(para,'a'), a = para.a; end
    if isfield(para,'b'), b = para.b; end
    if isfield(para,'c'), c = para.c; end
    if isfield(para,'d'), d = para.d; end
end
%--------------------------------------------------------------------------
% Descriptors
assert(size(fct_src,2)==size(fct_tar,2));
% Normalization
no = sqrt(diag(fct_src'*S1.A*fct_src))';
fct_src = fct_src ./ repmat(no, [S1.nv,1]);
no = sqrt(diag(fct_tar'*S2.A*fct_tar))';
fct_tar = fct_tar ./ repmat(no, [S2.nv,1]);
%--------------------------------------------------------------------------
% Multiplication Operators
numFct = size(fct_src,2);
OpSrc = cell(numFct,1);
OpTar = cell(numFct,1);
for i = 1:numFct
    OpSrc{i} = B1'*S1.A*(repmat(fct_src(:,i), [1,numEigsSrc]).*B1);
    OpTar{i} = B2'*S2.A*(repmat(fct_tar(:,i), [1,numEigsTar]).*B2);
end
Fct_src = B1'*S1.A*fct_src;
Fct_tar = B2'*S2.A*fct_tar;
% Orientation-preserving Operators
compute_all_OrientationOp = @(S,B,fct) ...
    cellfun(@(f) OrientationOp(S,B,f),mat2cell(fct,size(fct,1),ones(size(fct,2),1)),'un',0);
F11_all = compute_all_OrientationOp(S1,B1,fct_src);
F22_all = compute_all_OrientationOp(S2,B2,fct_tar);
%% normalize the operators

max_ev = max(max(Ev1),max(Ev2));
Ev1 = Ev1./max_ev;
Ev2 = Ev2./max_ev;

for i = 1:length(OpSrc)
    scale = max(norm(OpSrc{i}, 'fro'), norm(OpTar{i}, 'fro'));
    OpSrc{i} = OpSrc{i}/scale;
    OpTar{i} = OpTar{i}/scale;
    
    scale = max(norm(F11_all{i}, 'fro'), norm(F22_all{i}, 'fro'));
    F11_all{i} = F11_all{i}/scale;
    F22_all{i} = F22_all{i}/scale;   
end
%% all energy terms and the corresponding gradient
Dlb = (repmat(Ev2, [1,numEigsSrc]) - repmat(Ev1', [numEigsTar,1])).^2;
Dlb = Dlb/norm(Dlb, 'fro');


funOrient_symm = @(C) sum(cellfun(@(F11,F22) 0.5*norm(C*F11 + F22*C, 'fro')^2,F11_all, F22_all));
gradOrient_symm = @(C) sum(cell2mat(cellfun(@(F11,F22) reshape(C*(F11*F11') + F22'*C*F11 + F22*C*F11' + F22'*F22*C,[],1),...
    F11_all, F22_all,'un',0)),2);
% descriptors term
funDesp = @(C) 0.5*norm(C*Fct_src - Fct_tar,'fro')^2;
gradDesp = @(C) reshape((C*Fct_src - Fct_tar)*Fct_src',[],1);
% commutativity with descriptors
funCommDesp = @(C) sum(cell2mat(cellfun(@(X,Y) 0.5*norm(X*C - C*Y,'fro')^2, OpTar', OpSrc', 'UniformOutput', false)), 2);
gradCommDesp = @(C) sum(cell2mat(cellfun(@(X,Y) reshape(X'*(X*C - C*Y) - (X*C - C*Y)*Y',[],1), OpTar', OpSrc', 'UniformOutput', false)), 2);
% commutativity with LB
funCommLB = @(C) sum(sum((C.^2 .* Dlb)/2));
gradCommLB =@(C) reshape(C.*Dlb,[],1);

%% Fmap Computation
%--------------------------------------------------------------------------
constFct = sign(B1(1,1)*B2(1,1))*[sqrt(sum(S2.area)/sum(S1.area)); zeros(numEigsTar-1,1)];
F_lb = zeros(numEigsTar*numEigsSrc, 1); F_lb(1) = constFct(1);

myObj = @(C) a*funDesp(C) + b*funCommDesp(C) + c*funCommLB(C) + d*funOrient_symm(C);
myGrad = @(C) a*gradDesp(C) + b*gradCommDesp(C) + c*gradCommLB(C) + d*gradOrient_symm(C);
funObj = @(F) deal(myObj(reshape(F,numEigsTar, numEigsSrc)),...
    myGrad(reshape(F,numEigsTar, numEigsSrc)));
funProj = @(F) [constFct; F(numEigsTar+1:end)];
% options.maxIter = 1e3;
options.verbose = 1;
fprintf('Optimizing the self-symmetric functional map...');tic;
C12 = reshape(minConf_PQN(funObj, F_lb, funProj, options), [numEigsTar,numEigsSrc]);
t = toc; fprintf('done %.4fs.\n', t);
end
