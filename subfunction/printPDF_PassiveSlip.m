function printPDF_PassiveSlip(tcha,burn,id_zd,binw)

burn1 = size(tcha.smpasp,2) * burn / 100 + 1;
nline = size(tcha.smpasp,1) / 2;
id_zu = id_zd + nline;

figure(1); clf(1)
figure(1);histogram(tcha.smpasp(id_zd,burn1:end),'Normalization','pdf','BinWidth',binw,'FaceColor','b','LineStyle','none','FaceAlpha',0.5);hold on
figure(1);histogram(tcha.smpasp(id_zu,burn1:end),'Normalization','pdf','BinWidth',binw,'FaceColor','r','LineStyle','none','FaceAlpha',0.5);hold on
ax = gca;
ax.PlotBoxAspectRatio = [1,1,1];
ax.FontSize = 20;
ax.XLabel.String = 'Depth (km)';
ax.YLabel.String = 'PDF';

dep = 0:0.5:80;

pd_zd = fitdist(tcha.smpasp(id_zd,burn1:end)','Normal');
p_zd = pdf(pd_zd,dep);
plot(dep,p_zd,'Color','b','LineWidth',2)

pd_zu = fitdist(tcha.smpasp(id_zu,burn1:end)','Normal');
p_zu = pdf(pd_zu,dep);
plot(dep,p_zu,'Color','r','LineWidth',2)

saveas(1,['udline_pdf_',num2str(id_zd)])
print(['udline_pdf_',num2str(id_zd)],'-dpdf','-painters')

close all; clear
end