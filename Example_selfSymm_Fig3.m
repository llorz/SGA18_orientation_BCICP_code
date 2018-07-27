% This script reproduces the Fig.3 in the paper
%"Continuous and Orientation-preserving Correspondence via Functional Maps"
%
% In this example, we want to use the following function:
% compute_fMap_regular_with_orientationOp(S1,S2,B1,B2,Ev1,Ev2,fct1,fct2,'symmetric');
% to compute the self-symmetric map of a mesh, i.e., S1 = S2, B1 = B2 etc.
%%
clc; close all; clear;
mesh_dir = 'data/selfSymm/';
addpath(genpath('utils/'));
%% parameters
k = 10; % the number of basis used
numTimes = 1; % the time-scale to compute WKS descriptors
para.beta = 1; % the weight to control the orientation term (alpha_4)

% options to preprocess
mesh_options = {'IfComputeLB',true,'numEigs',k,... % compute k LB basis
    'IfComputeNormals',true,... % compute vtx normals for orientation term
    'IfComputeGeoDist',false};
%%
test_meshes = {'baby', 'homer', 'gorilla03', 'mesh025', '393', '177', '74'};

for i = 1:length(test_meshes)
    s_name = test_meshes{i};
    S = MESH.preprocess([mesh_dir,s_name],mesh_options{:}); % preprocess the mesh
    B = S.evecs(:,1:k); Ev = S.evals(1:k);
    fct = waveKernelSignature(B, Ev, S.A, numTimes); % compute the WKS descriptors
    [C_symm] = compute_fMap_regular_with_orientationOp(S,S,B,B,Ev,Ev,fct,fct,'symmetric',para); % use orientation-reversing operator
    C_symm = fMAP.icp_refine(B,B,C_symm,5);  % add 5 iteration of ICP
    T_symm = fMAP.fMap2pMap(B,B,C_symm); % convert to pointwise map
    
    % visualize the map
    switch s_name
        case {'baby','homer','gorrilla03'}
            plotOptions = {'CameraPos',[-90,0]};
        case 'mesh025'
            plotOptions = {'OverlayAxis','x','CameraPos',[-90,90]};
        case '393'
            plotOptions = {'OverlayAxis','z','CameraPos',[-90,0]};
        case '177'
            plotOptions = {'OverlayAxis','y','CameraPos',[0,90]};
        case '74'
            plotOptions = {'CameraPos',[-60,15]};    
    end
    figure(i);
    MESH.PLOT.visualize_map_colors(S,S,T_symm,'IfShowCoverage',false,plotOptions{:}); title(s_name); % visualize the map
end
