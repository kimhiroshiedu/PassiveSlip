function AutoCodingGmsh(model)
% This script generates Gmsh script for making trimesh
% Coded by Hiroshi Kimura 2019/11/07

% CodingPAC_NEjp(model)
% CodingPHS_Sagami_NEjp(model)
CodingPAC_SWjp(model)
CodingPHS_Sagami_SWjp(model)
CodingPHS(model)
end

%% Export routine
function [cp,cl,cs,cll] = Export2gmshfile(edge,fid,cp,cl,cs,cll,scale)

% Mesh size
fprintf(fid,'cellsize = %f;\n',scale);

% Points
for np = 1:size(edge,1)
    fprintf(fid,'Point(%i) =  {%f, %f, 0, cellsize};\n',cp, edge(np,:));
    cp = cp + 1;
end
cp = cp - size(edge,1);
fprintf(fid,'\n');

% Lines
for nl = 1:size(edge,1)
    if nl == size(edge,1)
        fprintf(fid,'Line(%i) = {%i, %i};\n',cl,cp,cp+1-size(edge,1));
    else
        fprintf(fid,'Line(%i) = {%i, %i};\n',cl,cp,cp+1             );
    end
    cp = cp + 1;
    cl = cl + 1;
end
cl = cl - size(edge,1);
fprintf(fid,'\n');

% Line loop
fprintf(fid,'Line Loop(%i) = {',cll);
for nll = 1:size(edge,1)
    if nll == size(edge,1)
        fprintf(fid,'%i'  ,cl);
    else
        fprintf(fid,'%i, ',cl);
    end
    cl = cl + 1;
end
fprintf(fid,'};\n');

% Plane surface
fprintf(fid,'Plane Surface(%i) = {%i};\n',cs,cll);
cll = cll + 1;
cs = cs + 1;

fprintf(fid,'\n');
end

%% Generating for PHS from Sagami to N-Ryukyu
function CodingPHS(model)
edges_IMP_NAN = load(['Meshes/model_',model,'/edges_IMP_NAN_rough.txt']);
edges_PHS_NAN = load(['Meshes/model_',model,'/edges_PHS_NAN_rough.txt']);
edges_PHS_ONN  = load(['Meshes/model_',model,'/edges_PHS_ONN_rough.txt']);
edges_PHS_ONC  = load(['Meshes/model_',model,'/edges_PHS_ONC_rough.txt']);
edges_PHS_ONS  = load(['Meshes/model_',model,'/edges_PHS_ONS_rough.txt']);

fid = fopen(['Meshes/model_',model,'/plate_',model,'_phs.geo'],'wt');
fprintf(fid,'// interplate_phs.geo\n\n');
fprintf(fid,'radius = 5.0;\npio2 = Pi/2;\n\n');

cp  =    1;
cl  =  501;
cll = 1001;
cs  = 2001;
% IMP \ NAN
scale = 0.21;
fprintf(fid,'// IMP to NAN\n');
[cp,cl,cs,cll] = Export2gmshfile(edges_IMP_NAN,fid,cp,cl,cs,cll,scale);

% PHS to NAN
scale = 0.21;
fprintf(fid,'// PHS to NAN\n');
[cp,cl,cs,cll] = Export2gmshfile(edges_PHS_NAN,fid,cp,cl,cs,cll,scale);

% PHS to ON-North
scale = 0.4;
fprintf(fid,'// PHS to ON-North\n');
[cp,cl,cs,cll] = Export2gmshfile(edges_PHS_ONN,fid,cp,cl,cs,cll,scale);

% PHS to ON-Center
scale = 1.2;
fprintf(fid,'// PHS to ON-Center\n');
[cp,cl,cs,cll] = Export2gmshfile(edges_PHS_ONC,fid,cp,cl,cs,cll,scale);

% PHS to ON-South
scale = 1.2;
fprintf(fid,'// PHS to ON-South\n');
[~,~,~,~] = Export2gmshfile(edges_PHS_ONS,fid,cp,cl,cs,cll,scale);

fclose(fid);
end

%% Generating for PHS Sagami trough (for SW Japan model)
function CodingPHS_Sagami_SWjp(model)
edges_IMP_OKH = load(['Meshes/model_',model,'/edges_IMP_OKH_rough.txt']);
edges_IOG_OKH = load(['Meshes/model_',model,'/edges_IOG_OKH_rough.txt']);

fid = fopen(['Meshes/model_',model,'/plate_',model,'_phssagami_swjp.geo'],'wt');
fprintf(fid,'// interplate_phssagami.geo\n\n');
fprintf(fid,'radius = 5.0;\npio2 = Pi/2;\n\n');

cp  =    1;
cl  =  501;
cll = 1001;
cs  = 2001;
% IMP \ OKH
scale = 0.5;
fprintf(fid,'// IMP to OKH\n');
[cp,cl,cs,cll] = Export2gmshfile(edges_IMP_OKH,fid,cp,cl,cs,cll,scale);

% Sagami Trough
scale = 0.5;
fprintf(fid,'// Sagami Trough\n');
[~,~,~,~] = Export2gmshfile(edges_IOG_OKH,fid,cp,cl,cs,cll,scale);

fclose(fid);
end

%% Generating for PHS Sagami trough (for NE Japan model)
function CodingPHS_Sagami_NEjp(model)
edges_IOG_THE = load(['Meshes/model_',model,'/edges_IOG_THE.txt']);
edges_IMP_THE = load(['Meshes/model_',model,'/edges_IMP_THE.txt']);
edges_IMP_THW = load(['Meshes/model_',model,'/edges_IMP_THW.txt']);

fid = fopen(['Meshes/model_',model,'/plate_',model,'_phssagami_nejp.geo'],'wt');
fprintf(fid,'// interplate_phssagami.geo\n\n');
fprintf(fid,'radius = 5.0;\npio2 = Pi/2;\n\n');

cp  =    1;
cl  =  501;
cll = 1001;
cs  = 2001;
% Sagami Trough
scale = 0.15;
fprintf(fid,'// Sagami Trough\n');
[cp,cl,cs,cll] = Export2gmshfile(edges_IOG_THE,fid,cp,cl,cs,cll,scale);

% IMP \ THE
scale = 0.15;
fprintf(fid,'// IMP to THE\n');
[cp,cl,cs,cll] = Export2gmshfile(edges_IMP_THE,fid,cp,cl,cs,cll,scale);

% IMP \ THW
scale = 0.15;
fprintf(fid,'// IMP to THW\n');
[~,~,~,~] = Export2gmshfile(edges_IMP_THW,fid,cp,cl,cs,cll,scale);

fclose(fid);
end

%% Generating for PAC subduction zone
function CodingPAC_NEjp(model)
edges_PAC_KUR = load(['Meshes/model_',model,'/edges_PAC_KUR.txt']);
edges_PAC_THE = load(['Meshes/model_',model,'/edges_PAC_THE.txt']);
edges_PAC_IOG = load(['Meshes/model_',model,'/edges_PAC_IOG.txt']);

fid = fopen(['Meshes/model_',model,'/plate_',model,'_pac_nejp.geo'],'wt');
fprintf(fid,'// interplate_pac.geo\n\n');
fprintf(fid,'radius = 5.0;\npio2 = Pi/2;\n\n');

cp  =    1;
cl  =  501;
cll = 1001;
cs  = 2001;

% Kuril trench
scale = 0.15;
fprintf(fid,'// Kuril trench\n');
[cp,cl,cs,cll] = Export2gmshfile(edges_PAC_KUR,fid,cp,cl,cs,cll,scale);

% Japan trench
scale = 0.15;
fprintf(fid,'// Japan trench\n');
[cp,cl,cs,cll] = Export2gmshfile(edges_PAC_THE,fid,cp,cl,cs,cll,scale);

% Izu-Ogasawara trench
scale = 0.15;
fprintf(fid,'// Izu-Ogasawara trench\n');
[~,~,~,~] = Export2gmshfile(edges_PAC_IOG,fid,cp,cl,cs,cll,scale);

fclose(fid);
end

%% Generating for PAC subduction zone (for SW Japan model)
function CodingPAC_SWjp(model)
edges_PAC_OKH = load(['Meshes/model_',model,'/edges_PAC_OKH_rough.txt']);
edges_PAC_IOG = load(['Meshes/model_',model,'/edges_PAC_IOG_rough.txt']);

fid = fopen(['Meshes/model_',model,'/plate_',model,'_pac_swjp.geo'],'wt');
fprintf(fid,'// interplate_pac.geo\n\n');
fprintf(fid,'radius = 5.0;\npio2 = Pi/2;\n\n');

cp  =    1;
cl  =  501;
cll = 1001;
cs  = 2001;

% Japan trench
scale = 2.0;
fprintf(fid,'// Kuril-Japan trench\n');
[cp,cl,cs,cll] = Export2gmshfile(edges_PAC_OKH,fid,cp,cl,cs,cll,scale);

% Izu-Ogasawara trench
scale = 2.0;
fprintf(fid,'// Izu-Ogasawara trench\n');
[~,~,~,~] = Export2gmshfile(edges_PAC_IOG,fid,cp,cl,cs,cll,scale);

fclose(fid);
end

