function AutoCodingGmsh
% This script generates Gmsh script for making trimesh
% Coded by Hiroshi Kimura 2019/11/07
CodingPAC
CodingPHS_Sagami
end

%% Generating for PHS Sagami trough
function CodingPHS_Sagami
load('MODEL_JP/BLOCK_Int_ne_japan/edges_6_7.txt')

fid = fopen('MODEL_JP/BLOCK_Int_ne_japan/interplate_phssagami.geo','wt');
fprintf(fid,'// interplate_phssagami.geo\n\n');
fprintf(fid,'radius = 5.0;\ncellsize = 0.15;\npio2 = Pi/2;\n\n');

% Sagami Trough
fprintf(fid,'// Sagami Trough\n');
countp = 1;
for np = 1:size(edges_6_7,1)
    fprintf(fid,'Point(%i) =  {%f, %f, 0, cellsize};\n',countp, edges_6_7(np,:));
    countp = countp + 1;
end
countp = countp - size(edges_6_7,1);
fprintf(fid,'\n');

countl = 501;
for nl = 1:size(edges_6_7,1)
    if nl == size(edges_6_7,1)
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1-size(edges_6_7,1));
    else
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1                  );
    end
    countp = countp + 1;
    countl = countl + 1;
end
countl = countl - size(edges_6_7,1);
fprintf(fid,'\n');

countll = 1001;
fprintf(fid,'Line Loop(%i) = {',countll);
for nll = 1:size(edges_6_7,1)
    if nll == size(edges_6_7,1)
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
function CodingPAC
load('MODEL_JP/BLOCK_Int_ne_japan/edges_5_8.txt')
load('MODEL_JP/BLOCK_Int_ne_japan/edges_6_8.txt')
load('MODEL_JP/BLOCK_Int_ne_japan/edges_7_8.txt')
% keyboard

fid = fopen('MODEL_JP/BLOCK_Int_ne_japan/interplate_pac.geo','wt');
fprintf(fid,'// interplate_pac.geo\n\n');
fprintf(fid,'radius = 5.0;\ncellsize = 0.15;\npio2 = Pi/2;\n\n');

% Kuril trench
fprintf(fid,'// Kuril trench\n');
countp = 1;
for np = 1:size(edges_5_8,1)
    fprintf(fid,'Point(%i) =  {%f, %f, 0, cellsize};\n',countp, edges_5_8(np,:));
    countp = countp + 1;
end
countp = countp - size(edges_5_8,1);
fprintf(fid,'\n');

countl = 501;
for nl = 1:size(edges_5_8,1)
    if nl == size(edges_5_8,1)
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1-size(edges_5_8,1));
    else
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1                  );
    end
    countp = countp + 1;
    countl = countl + 1;
end
countl = countl - size(edges_5_8,1);
fprintf(fid,'\n');

countll = 1001;
fprintf(fid,'Line Loop(%i) = {',countll);
for nll = 1:size(edges_5_8,1)
    if nll == size(edges_5_8,1)
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
for np = 1:size(edges_6_8,1)
    fprintf(fid,'Point(%i) =  {%f, %f, 0, cellsize};\n',countp, edges_6_8(np,:));
    countp = countp + 1;
end
countp = countp - size(edges_6_8,1);
fprintf(fid,'\n');

% countl = 501;
for nl = 1:size(edges_6_8,1)
    if nl == size(edges_6_8,1)
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1-size(edges_6_8,1));
    else
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1                  );
    end
    countp = countp + 1;
    countl = countl + 1;
end
countl = countl - size(edges_6_8,1);
fprintf(fid,'\n');

% countll = 1001;
fprintf(fid,'Line Loop(%i) = {',countll);
for nll = 1:size(edges_6_8,1)
    if nll == size(edges_6_8,1)
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
for np = 1:size(edges_7_8,1)
    fprintf(fid,'Point(%i) =  {%f, %f, 0, cellsize};\n',countp, edges_7_8(np,:));
    countp = countp + 1;
end
countp = countp - size(edges_7_8,1);
fprintf(fid,'\n');

% countl = 501;
for nl = 1:size(edges_7_8,1)
    if nl == size(edges_7_8,1)
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1-size(edges_7_8,1));
    else
        fprintf(fid,'Line(%i) = {%i, %i};\n',countl,countp,countp+1                  );
    end
    countp = countp + 1;
    countl = countl + 1;
end
countl = countl - size(edges_7_8,1);
fprintf(fid,'\n');

% countll = 1001;
fprintf(fid,'Line Loop(%i) = {',countll);
for nll = 1:size(edges_7_8,1)
    if nll == size(edges_7_8,1)
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

