red  = [           0: 1/32:1     ones(1,32)]';
green= [           0: 1/32:1 1-1/32:-1/32:0]';
blue = [ones(1,32) 1:-1/32:0               ]';
rwb = [red green blue];

% load('G:MasterResearch/inversion/PassiveSlip/Result_local/Test_01/tcha.mat')
%% correlation of zu, zd, pole, e
sample = [tcha.smpasp(:,:,1); tcha.smppol(:,:,1); tcha.smpine(:,:,1)];
ro = corr(sample(:,round(0.2*size(tcha.smpasp,2)):end)');
figure(1); clf(1)
imagesc(ro)
ax1 = gca;
ax1.PlotBoxAspectRatio = [1,1,1];
colormap(ax1,rwb);
colorbar
id_ud = size(tcha.smpasp,1);
xl_ud = xline(id_ud+0.5,'-y');
yl_ud = yline(id_ud+0.5,'-y');
id_pole = id_ud + size(tcha.smppol,1);
xl_pole = xline(id_pole+0.5,'-y');
yl_pole = yline(id_pole+0.5,'-y');
savefolder = 'Result_local/Test_01/figure/pdf_graphs';
saveas(1,fullfile(savefolder,'T_01_correlation_all'))
print(fullfile(savefolder,'T_01_correlation_all'),'-dpdf','-painters')
print(fullfile(savefolder,'T_01_correlation_all'),'-dpng')

%% correlation of dz, pole, e
dz = tcha.smpasp(1:size(tcha.smpasp,1)/2,:,1) - tcha.smpasp(size(tcha.smpasp,1)/2+1:end,:,1);
sample2 = [dz; tcha.smppol(:,:,1); tcha.smpine(:,:,1)];
ro2 = corr(sample2(:,round(0.2*size(tcha.smpasp,2)):end)');
figure(2); clf(2)
imagesc(ro2)
ax2 = gca;
ax2.PlotBoxAspectRatio = [1,1,1];
colormap(ax2,rwb);
colorbar
id_ud2 = size(dz,1);
xl_ud2 = xline(id_ud2+0.5,'-y');
yl_ud2 = yline(id_ud2+0.5,'-y');
id_pole2 = id_ud2 + size(tcha.smppol,1);
xl_pole2 = xline(id_pole2+0.5,'-y');
yl_pole2 = yline(id_pole2+0.5,'-y');
savefolder = 'Result_local/Test_01/figure/pdf_graphs';
saveas(2,fullfile(savefolder,'T_01_correlation_all_dz'))
print(fullfile(savefolder,'T_01_correlation_all_dz'),'-dpdf','-painters')
print(fullfile(savefolder,'T_01_correlation_all_dz'),'-dpng')

%% correlation of dz, pole, e (cut matrix)
% select id
cutid_ud = 9:27;
cutid_pole = [2,7:13];
cutid_E = [2,7:12];
% make id matrix
Cutid_ud = cutid_ud;
Cutid_pole = [];
for i=1:size(cutid_pole,2)
  Cutid_pole = [Cutid_pole,3*cutid_pole(i)-2,3*cutid_pole(i)-1,3*cutid_pole(i)];
end
Cutid_E = [];
for i=1:size(cutid_E,2)
  Cutid_E = [Cutid_E,3*cutid_E(i)-2,3*cutid_E(i)-1,3*cutid_E(i)];
end
sample_cut = [dz(Cutid_ud,:); tcha.smppol(Cutid_pole,:,1); tcha.smpine(Cutid_E,:,1)];
ro_cut = corr(sample_cut(:,round(0.2*size(tcha.smpasp,2)):end)');
figure(3); clf(3)
imagesc(ro_cut)
ax3 = gca;
ax3.PlotBoxAspectRatio = [1,1,1];
colormap(ax3,rwb);
colorbar
id_ud_cut = size(dz(Cutid_ud,:),1);
xl_ud_cut = xline(id_ud_cut+0.5,'-y');
yl_ud_cut = yline(id_ud_cut+0.5,'-y');
id_pole_cut = id_ud_cut + size(tcha.smppol(Cutid_pole,:,1),1);
xl_pole_cut = xline(id_pole_cut+0.5,'-y');
yl_pole_cut = yline(id_pole_cut+0.5,'-y');
savefolder = 'Result_local/Test_01/figure/pdf_graphs';
saveas(3,fullfile(savefolder,'T_01_correlation_cut_dz'))
print(fullfile(savefolder,'T_01_correlation_cut_dz'),'-dpdf','-painters')
print(fullfile(savefolder,'T_01_correlation_cut_dz'),'-dpng')
