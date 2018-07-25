function default_param = load_MeshPlot_default_params(S)
default_param.Title = S.name;

% If visualize a function on mesh
default_param.IfShowFunc = false;
default_param.func = [];
default_param.checkFunc = @(T) length(T)==S.nv; % check if the input function has the correct length

% If visualize a set of vertices
default_param.IfShowVtx = false;
default_param.landmarks = [];
default_param.Size_lmk = 50;
% default_param.Col_lmk = [1, 0.7, 0]; % orange
default_param.Col_lmk = [0.4, 0.58, 0.93]; % blue
% If visualize a set of edges
default_param.IfShowEdge = false;
default_param.edgeList = [];
default_param.checkEdgeList = @(edgeList) isnumeric(edgeList) && size(edgeList,2) == 2;
default_param.edgeID = [];
default_param.CameraPos = [0,0];
default_param.LineWidth = 1.2;

% parameters for trimesh
default_param.FaceColor = 'interp';
validFaceColor = {'interp','flat'};
default_param.checkFaceColor = @(x) any(validatestring(x,validFaceColor));

default_param.EdgeColor = 'none';
validEdgeColor = {'none','flat','interp'};
default_param.checkEdgeColor = @(x) any(validatestring(x,validEdgeColor));
default_param.FaceAlpha = 0.6;

end