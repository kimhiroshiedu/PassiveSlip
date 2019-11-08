function TriMeshGenerate
% This script generate tri_*.txt file of GMT format
% Coded by Hiroshi Kimura 2019/11/08
GeneratePACmesh;
GeneratePHSmesh;
GeneratePHSsagamimesh;
end

%% PHS Sagami
function GeneratePHSsagamimesh
interplate_phssagami
phs = load('~/MasterResearch/plate_data/plate_iwasaki/PHS_Plate/phs_regional/phs_2015_5a_r_2017.xyz');
F = scatteredInterpolant(phs(:,1),phs(:,2),phs(:,3));
msh.POS(:,3) = F(msh.POS(:,1),msh.POS(:,2));

fid = fopen('tri_phssagami.txt','wt');
for ntri = 1:size(msh.TRIANGLES,1)
    fprintf(fid,'%f %f %f\n',msh.POS(msh.TRIANGLES(ntri,1),1),msh.POS(msh.TRIANGLES(ntri,1),2),msh.POS(msh.TRIANGLES(ntri,1),3));
    fprintf(fid,'%f %f %f\n',msh.POS(msh.TRIANGLES(ntri,2),1),msh.POS(msh.TRIANGLES(ntri,2),2),msh.POS(msh.TRIANGLES(ntri,2),3));
    fprintf(fid,'%f %f %f\n',msh.POS(msh.TRIANGLES(ntri,3),1),msh.POS(msh.TRIANGLES(ntri,3),2),msh.POS(msh.TRIANGLES(ntri,3),3));
    fprintf(fid,'%f %f %f\n',msh.POS(msh.TRIANGLES(ntri,1),1),msh.POS(msh.TRIANGLES(ntri,1),2),msh.POS(msh.TRIANGLES(ntri,1),3));
    fprintf(fid,'>\n');
end
fclose(fid);

end

%% PHS Izu to Ryukyu
function GeneratePHSmesh
interplate_phs
phs = load('~/MasterResearch/plate_data/plate_iwasaki/PHS_Plate/phs_regional/phs_2015_5a_r_2017.xyz');
F = scatteredInterpolant(phs(:,1),phs(:,2),phs(:,3));
msh.POS(:,3) = F(msh.POS(:,1),msh.POS(:,2));

fid = fopen('tri_phs.txt','wt');
for ntri = 1:size(msh.TRIANGLES,1)
    fprintf(fid,'%f %f %f\n',msh.POS(msh.TRIANGLES(ntri,1),1),msh.POS(msh.TRIANGLES(ntri,1),2),msh.POS(msh.TRIANGLES(ntri,1),3));
    fprintf(fid,'%f %f %f\n',msh.POS(msh.TRIANGLES(ntri,2),1),msh.POS(msh.TRIANGLES(ntri,2),2),msh.POS(msh.TRIANGLES(ntri,2),3));
    fprintf(fid,'%f %f %f\n',msh.POS(msh.TRIANGLES(ntri,3),1),msh.POS(msh.TRIANGLES(ntri,3),2),msh.POS(msh.TRIANGLES(ntri,3),3));
    fprintf(fid,'%f %f %f\n',msh.POS(msh.TRIANGLES(ntri,1),1),msh.POS(msh.TRIANGLES(ntri,1),2),msh.POS(msh.TRIANGLES(ntri,1),3));
    fprintf(fid,'>\n');
end
fclose(fid);

end

%% PAC
function GeneratePACmesh
interplate_pac
pac = load('~/MasterResearch/plate_data/plate_iwasaki/PAC_Plate/pac_regional/pac_2017_4a.xyz');
F = scatteredInterpolant(pac(:,1),pac(:,2),pac(:,3));
msh.POS(:,3) = F(msh.POS(:,1),msh.POS(:,2));

fid = fopen('tri_pac.txt','wt');
for ntri = 1:size(msh.TRIANGLES,1)
    fprintf(fid,'%f %f %f\n',msh.POS(msh.TRIANGLES(ntri,1),1),msh.POS(msh.TRIANGLES(ntri,1),2),msh.POS(msh.TRIANGLES(ntri,1),3));
    fprintf(fid,'%f %f %f\n',msh.POS(msh.TRIANGLES(ntri,2),1),msh.POS(msh.TRIANGLES(ntri,2),2),msh.POS(msh.TRIANGLES(ntri,2),3));
    fprintf(fid,'%f %f %f\n',msh.POS(msh.TRIANGLES(ntri,3),1),msh.POS(msh.TRIANGLES(ntri,3),2),msh.POS(msh.TRIANGLES(ntri,3),3));
    fprintf(fid,'%f %f %f\n',msh.POS(msh.TRIANGLES(ntri,1),1),msh.POS(msh.TRIANGLES(ntri,1),2),msh.POS(msh.TRIANGLES(ntri,1),3));
    fprintf(fid,'>\n');
end
fclose(fid);

end
