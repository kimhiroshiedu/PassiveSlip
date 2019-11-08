function TriMeshGenerate
% This script generate tri_*.txt file of GMT format
% Coded by Hiroshi Kimura 2019/11/08
GeneratePACmesh_iwa;
GeneratePHSmesh_iwa;
GeneratePHSsagamimesh_iwa;
end

%% Iwasaki model
% PHS Sagami
function GeneratePHSsagamimesh_iwa
plate_iwasaki_phssagami;
if ispc
    phs = load('/MasterResearch/plate_data/plate_iwasaki/PHS_Plate/phs_regional/phs_2015_5a_r_2017.xyz');
else
    phs = load('~/MasterResearch/plate_data/plate_iwasaki/PHS_Plate/phs_regional/phs_2015_5a_r_2017.xyz');
end
F = scatteredInterpolant(phs(:,1),phs(:,2),phs(:,3));
msh.POS(:,3) = F(msh.POS(:,1),msh.POS(:,2));

file = 'Meshes/model_iwasaki/tri_phssagami.txt';
savetri(file,msh);

end

% PHS Izu to Ryukyu
function GeneratePHSmesh_iwa
plate_iwasaki_phs;
if ispc
    phs = load('/MasterResearch/plate_data/plate_iwasaki/PHS_Plate/phs_regional/phs_2015_5a_r_2017.xyz');
else
    phs = load('~/MasterResearch/plate_data/plate_iwasaki/PHS_Plate/phs_regional/phs_2015_5a_r_2017.xyz');
end
F = scatteredInterpolant(phs(:,1),phs(:,2),phs(:,3));
msh.POS(:,3) = F(msh.POS(:,1),msh.POS(:,2));

file = 'Meshes/model_iwasaki/tri_phs.txt';
savetri(file,msh);

end

% PAC
function GeneratePACmesh_iwa
plate_iwasaki_pac;
if ispc
    pac = load('/MasterResearch/plate_data/plate_iwasaki/PAC_Plate/pac_regional/pac_2017_4a.xyz');
else
    pac = load('~/MasterResearch/plate_data/plate_iwasaki/PAC_Plate/pac_regional/pac_2017_4a.xyz');
end
F = scatteredInterpolant(pac(:,1),pac(:,2),pac(:,3));
msh.POS(:,3) = F(msh.POS(:,1),msh.POS(:,2));

file = 'Meshes/model_iwasaki/tri_pac.txt';
savetri(file,msh);

end

%% Save part
function savetri(file,msh)
fid = fopen(file,'wt');
for ntri = 1:size(msh.TRIANGLES,1)
    fprintf(fid,'%f %f %f\n',msh.POS(msh.TRIANGLES(ntri,1),1),msh.POS(msh.TRIANGLES(ntri,1),2),msh.POS(msh.TRIANGLES(ntri,1),3));
    fprintf(fid,'%f %f %f\n',msh.POS(msh.TRIANGLES(ntri,2),1),msh.POS(msh.TRIANGLES(ntri,2),2),msh.POS(msh.TRIANGLES(ntri,2),3));
    fprintf(fid,'%f %f %f\n',msh.POS(msh.TRIANGLES(ntri,3),1),msh.POS(msh.TRIANGLES(ntri,3),2),msh.POS(msh.TRIANGLES(ntri,3),3));
    fprintf(fid,'%f %f %f\n',msh.POS(msh.TRIANGLES(ntri,1),1),msh.POS(msh.TRIANGLES(ntri,1),2),msh.POS(msh.TRIANGLES(ntri,1),3));
    fprintf(fid,'>\n');
end
fclose(fid);

end
