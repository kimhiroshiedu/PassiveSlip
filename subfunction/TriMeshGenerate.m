function TriMeshGenerate(model)
% This script generate tri_*.txt file of GMT format
% Coded by Hiroshi Kimura 2019/11/08
figure(10); clf(10)
GeneratePACmesh(model);
GeneratePHSmesh(model);
GeneratePHSsagamimesh(model);
hold off
end

%% PHS Sagami
function GeneratePHSsagamimesh(model)
switch model
    case 'iwasaki'
        plate_iwasaki_phssagami;
        if ispc
            phs = load('/MasterResearch/plate_data/plate_iwasaki/PHS_Plate/phs_regional/phs_2015_5a_r_2017.xyz');
        else
            phs = load('~/MasterResearch/plate_data/plate_iwasaki/PHS_Plate/phs_regional/phs_2015_5a_r_2017.xyz');
        end
    case 'hirose'
        plate_hirose_phssagami;
        if ispc
            phs = load('/MasterResearch/plate_data/plate_hirose/combined_contour_swjp/phs_hirose_slab1_combine.xyz');
        else
            phs = load('~/MasterResearch/plate_data/plate_hirose/combined_contour_swjp/phs_hirose_slab1_combine.xyz');
        end
    otherwise
        fprintf('Not applicable.\n'); return
end
F = scatteredInterpolant(phs(:,1),phs(:,2),phs(:,3));
msh.POS(:,3) = F(msh.POS(:,1),msh.POS(:,2));

file = ['Meshes/model_',model,'/tri_phssagami.txt'];
savetri(file,msh);

figure(10)
trisurf(msh.TRIANGLES(:,1:3),msh.POS(:,1),msh.POS(:,2),msh.POS(:,3))
hold on

end

%% PHS Izu to Ryukyu
function GeneratePHSmesh(model)
switch model
    case 'iwasaki'
        plate_iwasaki_phs;
        if ispc
            phs = load('/MasterResearch/plate_data/plate_iwasaki/PHS_Plate/phs_regional/phs_2015_5a_r_2017.xyz');
        else
            phs = load('~/MasterResearch/plate_data/plate_iwasaki/PHS_Plate/phs_regional/phs_2015_5a_r_2017.xyz');
        end
    case 'hirose'
        plate_hirose_phs;
        if ispc
            phs = load('/MasterResearch/plate_data/plate_hirose/combined_contour_swjp/phs_hirose_slab1_combine.xyz');
        else
            phs = load('~/MasterResearch/plate_data/plate_hirose/combined_contour_swjp/phs_hirose_slab1_combine.xyz');
        end
    otherwise
        fprintf('Not applicable.\n'); return
end

F = scatteredInterpolant(phs(:,1),phs(:,2),phs(:,3));
msh.POS(:,3) = F(msh.POS(:,1),msh.POS(:,2));

file = ['Meshes/model_',model,'/tri_phs.txt'];
savetri(file,msh);

figure(10)
trisurf(msh.TRIANGLES(:,1:3),msh.POS(:,1),msh.POS(:,2),msh.POS(:,3))
hold on

end

%% PAC
function GeneratePACmesh(model)
switch model
    case 'iwasaki'
        plate_iwasaki_pac;
        if ispc
            pac = load('/MasterResearch/plate_data/plate_iwasaki/PAC_Plate/pac_regional/pac_2017_4a.xyz');
        else
            pac = load('~/MasterResearch/plate_data/plate_iwasaki/PAC_Plate/pac_regional/pac_2017_4a.xyz');
        end
    case 'hirose'
        plate_hirose_pac;
        if ispc
            pac = load('/MasterResearch/plate_data/plate_hirose/PAC/pac_hirose_slab1_combine.xyz');
        else
            pac = load('~/MasterResearch/plate_data/plate_hirose/PAC/pac_hirose_slab1_combine.xyz');
        end
    otherwise
        fprintf('Not applicable.\n'); return
end

F = scatteredInterpolant(pac(:,1),pac(:,2),pac(:,3));
msh.POS(:,3) = F(msh.POS(:,1),msh.POS(:,2));

file = ['Meshes/model_',model,'/tri_pac.txt'];
savetri(file,msh);

figure(10)
trisurf(msh.TRIANGLES(:,1:3),msh.POS(:,1),msh.POS(:,2),msh.POS(:,3))
hold on

end

%% Save part
function savetri(file,msh)
fid = fopen(file,'wt');
for ntri = 1:size(msh.TRIANGLES,1)
    fprintf(fid,'%f %f %f\n',msh.POS(msh.TRIANGLES(ntri,1),1),msh.POS(msh.TRIANGLES(ntri,1),2),msh.POS(msh.TRIANGLES(ntri,1),3));
    fprintf(fid,'%f %f %f\n',msh.POS(msh.TRIANGLES(ntri,2),1),msh.POS(msh.TRIANGLES(ntri,2),2),msh.POS(msh.TRIANGLES(ntri,2),3));
    fprintf(fid,'%f %f %f\n',msh.POS(msh.TRIANGLES(ntri,3),1),msh.POS(msh.TRIANGLES(ntri,3),2),msh.POS(msh.TRIANGLES(ntri,3),3));
    fprintf(fid,'%f %f %f\n',msh.POS(msh.TRIANGLES(ntri,1),1),msh.POS(msh.TRIANGLES(ntri,1),2),msh.POS(msh.TRIANGLES(ntri,1),3));
    if ntri ~= size(msh.TRIANGLES,1)
        fprintf(fid,'>\n');
    end
end
fclose(fid);

end
