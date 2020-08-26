function printPDF_PassiveSlip(folder,blk,G,tcha,ud,T)
% tcha, blk and grn.mat must be loaded before using the function.

% Given parameters
sint   =  1;
binint =  2;
nud    = size(tcha.smpasp,1) / 2;
W      = blk(1).aline_zd(ud) - blk(1).aline_zu(ud);
sid    = round(size(tcha.smpasp,2)*tcha.burnin/100);

% Calc PDF of z_u and z_d
edges_d = blk(1).aline_zu(ud)    :binint:blk(1).aline_zd(ud)+2*W+binint;
edges_u = blk(1).aline_zu(ud)-2*W:binint:blk(1).aline_zd(ud)    +binint;
pdf_d = histcounts(tcha.smpasp(    ud,sid:sint:end,T),edges_d,'Normalization','pdf');
pdf_u = histcounts(tcha.smpasp(nud+ud,sid:sint:end,T),edges_u,'Normalization','pdf');
d_d = (edges_d(1:end-1) + edges_d(2:end)) ./ 2;
d_u = (edges_u(1:end-1) + edges_u(2:end)) ./ 2;
npdf_d = pdf_d ./ max([pdf_d,pdf_u]);
npdf_u = pdf_u ./ max([pdf_d,pdf_u]);

% Calc Confidence Intervals
sort_d = sort(tcha.smpasp(    ud,sid:sint:end,T));
sort_u = sort(tcha.smpasp(nud+ud,sid:sint:end,T));
ci95_d = [sort_d(round(size(sort_d,2)*0.025)),sort_d(round(size(sort_d,2)*0.975))];
ci95_u = [sort_u(round(size(sort_u,2)*0.025)),sort_u(round(size(sort_u,2)*0.975))];

% Calc PDF of Pc
ma = zeros(size(tcha.smpasp,1),1);
ma(ud) = 1;
idtri = G(1).zd * ma;
idtri = idtri == 1;
dep = G(1).zc(idtri);
pc  = tcha.aveaid(idtri,:,T);
x = blk(1).aline_zu(ud):1:blk(1).aline_zd(ud);
y = interp1(dep,pc,x,'linear','extrap');

%% Plot whole depth range
figure(230); clf(230)
subplot(1,2,1)
stairs(d_d,npdf_d,'-b','LineWidth',1); hold on
stairs(d_u,npdf_u,'-r','LineWidth',1);
xl_mean_zd = xline(tcha.aveasp(    ud,1,T),'-b');
xl_mean_zu = xline(tcha.aveasp(nud+ud,1,T),'-r');
xl_ci95_zd = xline(ci95_d(1),'--b'); xline(ci95_d(2),'--b');
xl_ci95_zu = xline(ci95_u(1),'--r'); xline(ci95_u(2),'--r');
plot(x,y,'Color','g','LineWidth',1)
xl_bot = xline(blk(1).aline_zd(ud),'-k');
xl_top = xline(blk(1).aline_zu(ud),'-k');
legend([xl_mean_zu, xl_mean_zd, xl_ci95_zu, xl_ci95_zd],...
    {['Mean: ',num2str(tcha.aveasp(nud+ud,1,T),'%-4.0f'),'km'],...
     ['Mean: ',num2str(tcha.aveasp(    ud,1,T),'%-4.0f'),'km'],...
     ['95% CI: ',num2str(diff(ci95_u),'%-4.0f'),'km'],...
     ['95% CI: ',num2str(diff(ci95_d),'%-4.0f'),'km']});
view(90,90)
hold off

xl_bot.LineWidth = 1; xl_bot.Alpha = 1;
xl_top.LineWidth = 1; xl_top.Alpha = 1;
ax1 = gca;
ax1.PlotBoxAspectRatio = [2,1,1];
ax1.XLim = [blk(1).aline_zu(ud)-2*W,blk(1).aline_zd(ud)+2*W];
ax1.YLim = [0,1.1];
ax1.XColor = 'k'; ax1.YColor = 'k';
ax1.FontSize = 11;
ax1.XLabel.String = 'Depth (km)'; ax1.XLabel.Color = 'k';
ax1.YLabel.String = 'PDF'       ; ax1.YLabel.Color = 'k';


%% Plot zoomed depth range
figure(230)
subplot(1,2,2)
stairs(d_d,npdf_d,'-b','LineWidth',1); hold on
stairs(d_u,npdf_u,'-r','LineWidth',1);
xl_d = xline(tcha.aveasp(    ud,1,T),'-b',[num2str(tcha.aveasp(    ud,1,T),'%-4.0f\n'),'km']);
xl_u = xline(tcha.aveasp(nud+ud,1,T),'-r',[num2str(tcha.aveasp(nud+ud,1,T),'%-4.0f\n'),'km']);
plot(x,y,'Color','g','LineWidth',1)
view(90,90)
hold off

xl_d.LabelHorizontalAlignment =  'right'; xl_d.LabelVerticalAlignment   = 'middle';
xl_u.LabelHorizontalAlignment =   'left'; xl_u.LabelVerticalAlignment   = 'middle';
ax2 = gca;
ax2.PlotBoxAspectRatio = [2,1,1];
ax2.XLim = [blk(1).aline_zu(ud),blk(1).aline_zd(ud)];
ax2.YLim = [0,1.1];
ax2.XColor = 'k'; ax2.YColor = 'k';
ax2.FontSize = 11;
ax2.XLabel.String = 'Depth (km)'; ax2.XLabel.Color = 'k';
ax2.YLabel.String = 'PDF'       ; ax2.YLabel.Color = 'k';

% return
%% Save figures
savefolder = fullfile(pwd,folder,'figure','pdf_graphs');
if exist(savefolder,'dir') ~= 7; mkdir(savefolder); end
saveas(230,fullfile(savefolder,['T',num2str(T,'%02i'),'_pdf_ud_',num2str(ud)]))
print(fullfile(savefolder,['T',num2str(T,'%02i'),'_pdf_ud_',num2str(ud)]),'-dpdf','-painters')
print(fullfile(savefolder,['T',num2str(T,'%02i'),'_pdf_ud_',num2str(ud)]),'-dpng')
end

function printPDF_PassiveSlip_old(tcha,burn,id_zd,binw)

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

dep = 0:0.05:80;

pd_zd = fitdist(tcha.smpasp(id_zd,burn1:end)','Normal');
p_zd = pdf(pd_zd,dep);
plot(dep,p_zd,'Color','b','LineWidth',2)

pd_zu = fitdist(tcha.smpasp(id_zu,burn1:end)','Normal');
p_zu = pdf(pd_zu,dep);
plot(dep,p_zu,'Color','r','LineWidth',2)

saveas(1,['udline_pdf_',num2str(id_zd)])
print(['udline_pdf_',num2str(id_zd)],'-dpdf','-painters')

end