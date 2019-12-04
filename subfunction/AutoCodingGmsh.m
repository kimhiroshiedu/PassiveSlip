function AutoCodingGmsh(model)
% This script generates Gmsh script for making trimesh
% Coded by Hiroshi Kimura 2019/11/07

CodingPAC(model)
CodingPHS_Sagami(model)
CodingPHS(model)
end

%% Generating for PHS fron Sagami to N-Ryukyu
function CodingPHS(model)
edges_IMP_NEJ = load(['Meshes/model_',model,'/edges_IMP_NEJ.txt']);
edges_IMP_OK  = load(['Meshes/model_',model,'/edges_IMP_OK.txt']);
edges_IMP_SWJ = load(['Meshes/model_',model,'/edges_IMP_SWJ.txt']);
edges_PHS_SWJ = load(['Meshes/model_',model,'/edges_PHS_SWJ.txt']);
edges_PHS_ON  = load(['Meshes/model_',model,'/edges_PHS_ON.txt']);

% keyboard

fid = fopen(['Meshes/model_',model,'/plate_',model,'_phs.geo'],'wt');
fprintf(fid,'// interplate_pac.geo\n\n');
fprintf(fid,'radius = 5.0;\ncellsize = 0.15;\npio2 = Pi/2;\n\n');

% IMP east
fprintf(fid,'// IMP east\n');
countp = 1;
for np = 1:size(edges_IMP_NEJ,1)
    fprintf(fid,'Point(%i) =  {%f, %f, 0, cellsize};\n',countp, edges_IMP_NEJ(np,:));
    countp = countp + 1;
end
countp = countp - size(edges_IMP_NEJ,1);
fprintf(fid,'\n');

countl = 501;
for nl = 1:size(edges_IMP_NEJ,1)
    if nl == size(edges_IMP_NEJ,1)
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1-size(edges_IMP_NEJ,1));
    else
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1                  );
    end
    countp = countp + 1;
    countl = countl + 1;
end
countl = countl - size(edges_IMP_NEJ,1);
fprintf(fid,'\n');

countll = 1001;
fprintf(fid,'Line Loop(%i) = {',countll);
for nll = 1:size(edges_IMP_NEJ,1)
    if nll == size(edges_IMP_NEJ,1)
        fprintf(fid,'%i'  ,countl);
    else
        fprintf(fid,'%i, ',countl);
    end
    countl = countl + 1;
end
fprintf(fid,'};\n');

counts = 2001;
fprintf(fid,'Plane Surface(%i) = {%i};\n\n',counts,countll);
countll = countll + 1;
counts = counts + 1;


% IMP west
fprintf(fid,'// IMP west\n');
% countp = 1;
for np = 1:size(edges_IMP_OK,1)
    fprintf(fid,'Point(%i) =  {%f, %f, 0, cellsize};\n',countp, edges_IMP_OK(np,:));
    countp = countp + 1;
end
countp = countp - size(edges_IMP_OK,1);
fprintf(fid,'\n');

% countl = 501;
for nl = 1:size(edges_IMP_OK,1)
    if nl == size(edges_IMP_OK,1)
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1-size(edges_IMP_OK,1));
    else
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1                  );
    end
    countp = countp + 1;
    countl = countl + 1;
end
countl = countl - size(edges_IMP_OK,1);
fprintf(fid,'\n');

% countll = 1001;
fprintf(fid,'Line Loop(%i) = {',countll);
for nll = 1:size(edges_IMP_OK,1)
    if nll == size(edges_IMP_OK,1)
        fprintf(fid,'%i'  ,countl);
    else
        fprintf(fid,'%i, ',countl);
    end
    countl = countl + 1;
end
fprintf(fid,'};\n');

% counts = 2001;
fprintf(fid,'Plane Surface(%i) = {%i};\n',counts,countll);
countll = countll + 1;
counts = counts + 1;


% IMP to FA
fprintf(fid,'// IMP to FA\n');
% countp = 1;
for np = 1:size(edges_IMP_SWJ,1)
    fprintf(fid,'Point(%i) =  {%f, %f, 0, cellsize};\n',countp, edges_IMP_SWJ(np,:));
    countp = countp + 1;
end
countp = countp - size(edges_IMP_SWJ,1);
fprintf(fid,'\n');

% countl = 501;
for nl = 1:size(edges_IMP_SWJ,1)
    if nl == size(edges_IMP_SWJ,1)
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1-size(edges_IMP_SWJ,1));
    else
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1                  );
    end
    countp = countp + 1;
    countl = countl + 1;
end
countl = countl - size(edges_IMP_SWJ,1);
fprintf(fid,'\n');

% countll = 1001;
fprintf(fid,'Line Loop(%i) = {',countll);
for nll = 1:size(edges_IMP_SWJ,1)
    if nll == size(edges_IMP_SWJ,1)
        fprintf(fid,'%i'  ,countl);
    else
        fprintf(fid,'%i, ',countl);
    end
    countl = countl + 1;
end
fprintf(fid,'};\n');

% counts = 2001;
fprintf(fid,'Plane Surface(%i) = {%i};\n',counts,countll);
countll = countll + 1;
counts = counts + 1;


% PHS to FA
fprintf(fid,'// PHS to FA\n');
% countp = 1;
for np = 1:size(edges_PHS_SWJ,1)
    fprintf(fid,'Point(%i) =  {%f, %f, 0, cellsize};\n',countp, edges_PHS_SWJ(np,:));
    countp = countp + 1;
end
countp = countp - size(edges_PHS_SWJ,1);
fprintf(fid,'\n');

% countl = 501;
for nl = 1:size(edges_PHS_SWJ,1)
    if nl == size(edges_PHS_SWJ,1)
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1-size(edges_PHS_SWJ,1));
    else
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1                  );
    end
    countp = countp + 1;
    countl = countl + 1;
end
countl = countl - size(edges_PHS_SWJ,1);
fprintf(fid,'\n');

% countll = 1001;
fprintf(fid,'Line Loop(%i) = {',countll);
for nll = 1:size(edges_PHS_SWJ,1)
    if nll == size(edges_PHS_SWJ,1)
        fprintf(fid,'%i'  ,countl);
    else
        fprintf(fid,'%i, ',countl);
    end
    countl = countl + 1;
end
fprintf(fid,'};\n');

% counts = 2001;
fprintf(fid,'Plane Surface(%i) = {%i};\n',counts,countll);
countll = countll + 1;
counts = counts + 1;


% PHS to N-ON
fprintf(fid,'// PHS to N-ON\n');
% countp = 1;
for np = 1:size(edges_PHS_ON,1)
    fprintf(fid,'Point(%i) =  {%f, %f, 0, cellsize};\n',countp, edges_PHS_ON(np,:));
    countp = countp + 1;
end
countp = countp - size(edges_PHS_ON,1);
fprintf(fid,'\n');

% countl = 501;
for nl = 1:size(edges_PHS_ON,1)
    if nl == size(edges_PHS_ON,1)
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1-size(edges_PHS_ON,1));
    else
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1                  );
    end
    countp = countp + 1;
    countl = countl + 1;
end
countl = countl - size(edges_PHS_ON,1);
fprintf(fid,'\n');

% countll = 1001;
fprintf(fid,'Line Loop(%i) = {',countll);
for nll = 1:size(edges_PHS_ON,1)
    if nll == size(edges_PHS_ON,1)
        fprintf(fid,'%i'  ,countl);
    else
        fprintf(fid,'%i, ',countl);
    end
    countl = countl + 1;
end
fprintf(fid,'};\n');

% counts = 2001;
fprintf(fid,'Plane Surface(%i) = {%i};\n',counts,countll);
countll = countll + 1;
counts = counts + 1;


fclose(fid);
end

%% Generating for PHS Sagami trough
function CodingPHS_Sagami(model)
edges_IOG_NEJ = load(['Meshes/model_',model,'/edges_IOG_NEJ.txt']);

fid = fopen(['Meshes/model_',model,'/plate_',model,'_phssagami.geo'],'wt');
fprintf(fid,'// interplate_phssagami.geo\n\n');
fprintf(fid,'radius = 5.0;\ncellsize = 0.15;\npio2 = Pi/2;\n\n');

% Sagami Trough
fprintf(fid,'// Sagami Trough\n');
countp = 1;
for np = 1:size(edges_IOG_NEJ,1)
    fprintf(fid,'Point(%i) =  {%f, %f, 0, cellsize};\n',countp, edges_IOG_NEJ(np,:));
    countp = countp + 1;
end
countp = countp - size(edges_IOG_NEJ,1);
fprintf(fid,'\n');

countl = 501;
for nl = 1:size(edges_IOG_NEJ,1)
    if nl == size(edges_IOG_NEJ,1)
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1-size(edges_IOG_NEJ,1));
    else
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1                  );
    end
    countp = countp + 1;
    countl = countl + 1;
end
countl = countl - size(edges_IOG_NEJ,1);
fprintf(fid,'\n');

countll = 1001;
fprintf(fid,'Line Loop(%i) = {',countll);
for nll = 1:size(edges_IOG_NEJ,1)
    if nll == size(edges_IOG_NEJ,1)
        fprintf(fid,'%i'  ,countl);
    else
        fprintf(fid,'%i, ',countl);
    end
    countl = countl + 1;
end
fprintf(fid,'};\n');

counts = 2001;
fprintf(fid,'Plane Surface(%i) = {%i};\n\n',counts,countll);
countll = countll + 1;
counts = counts + 1;

fclose(fid);
end

%% Generating for PAC subduction zone
function CodingPAC(model)
edges_PAC_KUR = load(['Meshes/model_',model,'/edges_PAC_KUR.txt']);
edges_PAC_NEJ = load(['Meshes/model_',model,'/edges_PAC_NEJ.txt']);
edges_PAC_IOG = load(['Meshes/model_',model,'/edges_PAC_IOG.txt']);
% keyboard

fid = fopen(['Meshes/model_',model,'/plate_',model,'_pac.geo'],'wt');
fprintf(fid,'// interplate_pac.geo\n\n');
fprintf(fid,'radius = 5.0;\ncellsize = 0.15;\npio2 = Pi/2;\n\n');

% Kuril trench
fprintf(fid,'// Kuril trench\n');
countp = 1;
for np = 1:size(edges_PAC_KUR,1)
    fprintf(fid,'Point(%i) =  {%f, %f, 0, cellsize};\n',countp, edges_PAC_KUR(np,:));
    countp = countp + 1;
end
countp = countp - size(edges_PAC_KUR,1);
fprintf(fid,'\n');

countl = 501;
for nl = 1:size(edges_PAC_KUR,1)
    if nl == size(edges_PAC_KUR,1)
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1-size(edges_PAC_KUR,1));
    else
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1                  );
    end
    countp = countp + 1;
    countl = countl + 1;
end
countl = countl - size(edges_PAC_KUR,1);
fprintf(fid,'\n');

countll = 1001;
fprintf(fid,'Line Loop(%i) = {',countll);
for nll = 1:size(edges_PAC_KUR,1)
    if nll == size(edges_PAC_KUR,1)
        fprintf(fid,'%i'  ,countl);
    else
        fprintf(fid,'%i, ',countl);
    end
    countl = countl + 1;
end
fprintf(fid,'};\n');

counts = 2001;
fprintf(fid,'Plane Surface(%i) = {%i};\n\n',counts,countll);
countll = countll + 1;
counts = counts + 1;


% Japan trench
fprintf(fid,'// Japan trench\n');
% countp = 1;
for np = 1:size(edges_PAC_NEJ,1)
    fprintf(fid,'Point(%i) =  {%f, %f, 0, cellsize};\n',countp, edges_PAC_NEJ(np,:));
    countp = countp + 1;
end
countp = countp - size(edges_PAC_NEJ,1);
fprintf(fid,'\n');

% countl = 501;
for nl = 1:size(edges_PAC_NEJ,1)
    if nl == size(edges_PAC_NEJ,1)
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1-size(edges_PAC_NEJ,1));
    else
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1                  );
    end
    countp = countp + 1;
    countl = countl + 1;
end
countl = countl - size(edges_PAC_NEJ,1);
fprintf(fid,'\n');

% countll = 1001;
fprintf(fid,'Line Loop(%i) = {',countll);
for nll = 1:size(edges_PAC_NEJ,1)
    if nll == size(edges_PAC_NEJ,1)
        fprintf(fid,'%i'  ,countl);
    else
        fprintf(fid,'%i, ',countl);
    end
    countl = countl + 1;
end
fprintf(fid,'};\n');

% counts = 2001;
fprintf(fid,'Plane Surface(%i) = {%i};\n',counts,countll);
countll = countll + 1;
counts = counts + 1;


% Izu-Ogasawara trench
fprintf(fid,'// Izu-Ogasawara trench\n');
% countp = 1;
for np = 1:size(edges_PAC_IOG,1)
    fprintf(fid,'Point(%i) =  {%f, %f, 0, cellsize};\n',countp, edges_PAC_IOG(np,:));
    countp = countp + 1;
end
countp = countp - size(edges_PAC_IOG,1);
fprintf(fid,'\n');

% countl = 501;
for nl = 1:size(edges_PAC_IOG,1)
    if nl == size(edges_PAC_IOG,1)
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1-size(edges_PAC_IOG,1));
    else
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1                  );
    end
    countp = countp + 1;
    countl = countl + 1;
end
countl = countl - size(edges_PAC_IOG,1);
fprintf(fid,'\n');

% countll = 1001;
fprintf(fid,'Line Loop(%i) = {',countll);
for nll = 1:size(edges_PAC_IOG,1)
    if nll == size(edges_PAC_IOG,1)
        fprintf(fid,'%i'  ,countl);
    else
        fprintf(fid,'%i, ',countl);
    end
    countl = countl + 1;
end
fprintf(fid,'};\n');

% counts = 2001;
fprintf(fid,'Plane Surface(%i) = {%i};\n',counts,countll);
countll = countll + 1;
counts = counts + 1;

fclose(fid);
end

