function printPDF_PassiveSlip(folder,blk,G,tcha,udlist,T)
% tcha, blk and grn.mat must be loaded before using the function.

% Given parameters
sint   =  1;
binint =  2;
nud    = size(tcha.smpasp,1) / 2;
sid    = round(size(tcha.smpasp,2)*tcha.burnin/100);
prm.sint   = sint;
prm.binint = binint;
prm.nud    = nud;
prm.sid    = sid;

% Data prepare
maxpdf = 0;
for ud = udlist
  W = blk(1).aline_zd(ud) - blk(1).aline_zu(ud);
  udl(ud).W = W;
  udl(ud).edges_d = blk(1).aline_zu(ud)    :binint:blk(1).aline_zd(ud)+2*W+binint;
  udl(ud).edges_u = blk(1).aline_zu(ud)-2*W:binint:blk(1).aline_zd(ud)    +binint;
  z_d = tcha.smpasp(    ud,sid:sint:end,T); pick_z_d = z_d;
  z_u = tcha.smpasp(nud+ud,sid:sint:end,T); pick_z_u = z_u;
  udl(ud).pdf_d = histcounts(z_d,udl(ud).edges_d,'Normalization','pdf');
  udl(ud).pdf_u = histcounts(z_u,udl(ud).edges_u,'Normalization','pdf');
  id_cross = z_d < z_u;
  pick_z_d(id_cross) = 99999;
  pick_z_u(id_cross) = 99999;
  udl(ud).pdf_dcross = histcounts(pick_z_d,udl(ud).edges_d,'Normalization','pdf');
  udl(ud).pdf_ucross = histcounts(pick_z_u,udl(ud).edges_u,'Normalization','pdf');
  if sum(~id_cross) == 0
    udl(ud).pdf_dcross(:) = 0;
    udl(ud).pdf_ucross(:) = 0;
  end
  udl(ud).mean_pick_z_d = mean(pick_z_d);
  udl(ud).mean_pick_z_u = mean(pick_z_u);
  maxpdf_ud = max([udl(ud).pdf_d,udl(ud).pdf_u]);
  if maxpdf_ud > maxpdf; maxpdf = maxpdf_ud; end
end

% Make, save figures
for ud = udlist
  printpdf_ud(folder,blk,G,tcha,prm,udl,maxpdf,ud,T)
end

savelog(folder,udlist)
end

function savelog(folder,udlist)
savefolder = fullfile(pwd,folder,'figure','pdf_graphs');
if exist(savefolder,'dir') ~= 7; mkdir(savefolder); end
fid = fopen(fullfile(savefolder,'makepdflog.txt'),'a+');
dt = datetime('now','TimeZone','Asia/Tokyo');
DateString = datestr(dt,'yyyy/mm/dd HH:MM:SS');
fprintf(fid,'%s\n ud_No.:',DateString);
for ud = udlist
  fprintf(fid,' %i',ud);
end
fprintf(fid,'\n\n');
fclose(fid);
end

function printpdf_ud(folder,blk,G,tcha,prm,udl,maxpdf,ud,T)
sint   = prm.sint;
binint = prm.binint;
nud    = prm.nud;
sid    = prm.sid;
W      = udl(ud).W;
% Calc PDF of z_u and z_d
edges_d = udl(ud).edges_d;
edges_u = udl(ud).edges_u;
pdf_d = udl(ud).pdf_d;
pdf_u = udl(ud).pdf_u;
pdf_dcross = udl(ud).pdf_dcross;
pdf_ucross = udl(ud).pdf_ucross;
d_d = (edges_d(1:end-1) + edges_d(2:end)) ./ 2;
d_u = (edges_u(1:end-1) + edges_u(2:end)) ./ 2;
npdf_d = pdf_d ./ max([pdf_d,pdf_u]);
npdf_u = pdf_u ./ max([pdf_d,pdf_u]);

% Calc mean and 95% CI
zd_mean = tcha.aveasp(    ud,1,T);
zu_mean = tcha.aveasp(nud+ud,1,T);
sort_d = sort(tcha.smpasp(    ud,sid:sint:end,T));
sort_u = sort(tcha.smpasp(nud+ud,sid:sint:end,T));
ci95_d = [sort_d(round(size(sort_d,2)*0.025)),sort_d(round(size(sort_d,2)*0.975))];
ci95_u = [sort_u(round(size(sort_u,2)*0.025)),sort_u(round(size(sort_u,2)*0.975))];

% Calc PDF of Pc
ma = zeros(size(tcha.smpasp,1),1);
ma(ud) = 1;
idtri = G(1).zd * ma;
idtri = idtri == 1;
tri_dep = G(1).zc(idtri);
tri_pl  = tcha.aveaid(idtri,:,T);
pl_x = blk(1).aline_zu(ud):binint:blk(1).aline_zd(ud);
pl_y = interp1(tri_dep,tri_pl,pl_x,'linear');
ids = find(~isnan(pl_y));
pl_y(       1:ids(1)) = pl_y(ids( 1 ));
pl_y(ids(end):   end) = pl_y(ids(end));

%% Plot whole depth range
fig = figure(230); clf(230)
fig.Position = [488 200 560 520];
subplot(2,1,2)
yyaxis right
hd  = histogram('BinEdges',edges_d,'BinCounts',pdf_d,'DisplayStyle','stairs','EdgeColor','b','EdgeAlpha',1,'LineWidth',1); hold on
hu  = histogram('BinEdges',edges_u,'BinCounts',pdf_u,'DisplayStyle','stairs','EdgeColor','r','EdgeAlpha',1,'LineWidth',1);
hdc = histogram('BinEdges',edges_d,'BinCounts',pdf_dcross,'DisplayStyle','bar','EdgeColor','none','FaceColor','b','FaceAlpha',0.3);
huc = histogram('BinEdges',edges_u,'BinCounts',pdf_ucross,'DisplayStyle','bar','EdgeColor','none','FaceColor','r','FaceAlpha',0.3);
xl_mean_zd = xline(zd_mean,'-b',[num2str(zd_mean,'%-4.0f'),'(',num2str(diff(ci95_d),'%-3.0f'),') km'],'FontName','Helvetica');
xl_mean_zu = xline(zu_mean,'-r',[num2str(zu_mean,'%-4.0f'),'(',num2str(diff(ci95_u),'%-3.0f'),') km'],'FontName','Helvetica');
xl_ci95_zd = xline(ci95_d(1),':b','Linewidth',1); xline(ci95_d(2),':b','Linewidth',1);
xl_ci95_zu = xline(ci95_u(1),':r','Linewidth',1); xline(ci95_u(2),':r','Linewidth',1);
yyaxis left
plot(pl_x,pl_y,'Color','g','LineWidth',1)
xl_bot = xline(blk(1).aline_zd(ud),'-k','Linewidth',1,'Alpha',1);
xl_top = xline(blk(1).aline_zu(ud),'-k','Linewidth',1,'Alpha',1);
ax1 = gca;
ax1.XLim = [blk(1).aline_zu(ud)-2*W,blk(1).aline_zd(ud)+2*W];
ax1.XLabel.String = 'Depth (km)';
ax1.XLabel.Rotation = 180;
ax1.XLabel.VerticalAlignment = "bottom";
ax1 = graph_settings(ax1,maxpdf);
% legend([xl_mean_zu, xl_mean_zd, xl_ci95_zu, xl_ci95_zd],...
%     {['Mean: ',num2str(tcha.aveasp(nud+ud,1,T),'%-4.0f'),'km'],...
%      ['Mean: ',num2str(tcha.aveasp(    ud,1,T),'%-4.0f'),'km'],...
%      ['95% CI: ',num2str(diff(ci95_u),'%-4.0f'),'km'],...
%      ['95% CI: ',num2str(diff(ci95_d),'%-4.0f'),'km']});
hold off

%% Plot zoomed depth range
figure(230);
subplot(2,1,1)
yyaxis right
hd  = histogram('BinEdges',edges_d,'BinCounts',pdf_d,'DisplayStyle','stairs','EdgeColor','b','EdgeAlpha',1,'LineWidth',1); hold on
hu  = histogram('BinEdges',edges_u,'BinCounts',pdf_u,'DisplayStyle','stairs','EdgeColor','r','EdgeAlpha',1,'LineWidth',1);
hdc = histogram('BinEdges',edges_d,'BinCounts',pdf_dcross,'DisplayStyle','bar','EdgeColor','none','FaceColor','b','FaceAlpha',0.3);
huc = histogram('BinEdges',edges_u,'BinCounts',pdf_ucross,'DisplayStyle','bar','EdgeColor','none','FaceColor','r','FaceAlpha',0.3);
xl_mean_zd = xline(zd_mean,'-b',[num2str(zd_mean,'%-4.0f'),'(',num2str(diff(ci95_d),'%-3.0f'),') km'],'FontName','Helvetica');
xl_mean_zu = xline(zu_mean,'-r',[num2str(zu_mean,'%-4.0f'),'(',num2str(diff(ci95_u),'%-3.0f'),') km'],'FontName','Helvetica');
xl_ci95_zd = xline(ci95_d(1),':b','Linewidth',1); xline(ci95_d(2),':b','Linewidth',1);
xl_ci95_zu = xline(ci95_u(1),':r','Linewidth',1); xline(ci95_u(2),':r','Linewidth',1);
yyaxis left
plot(pl_x,pl_y,'Color','g','LineWidth',1)

ax2 = gca;
ax2.XLim = [blk(1).aline_zu(ud),blk(1).aline_zd(ud)];
ax2 = graph_settings(ax2,maxpdf);
hold off

%% Save figures
savefolder = fullfile(pwd,folder,'figure','pdf_graphs');
if exist(savefolder,'dir') ~= 7; mkdir(savefolder); end
saveas(230,fullfile(savefolder,['T',num2str(T,'%02i'),'_pdf_ud_',num2str(ud)]))
print(fullfile(savefolder,['T',num2str(T,'%02i'),'_pdf_ud_',num2str(ud)]),'-dpdf','-painters')
print(fullfile(savefolder,['T',num2str(T,'%02i'),'_pdf_ud_',num2str(ud)]),'-dpng')

end

function [ax] = graph_settings(ax,maxpdf)
% axis
ax.FontName = 'Helvetica';
ax.PlotBoxAspectRatio = [2,1,1];
ax.FontSize = 11;
% x-axis
ax.XColor = 'k'; 
ax.XTickLabelRotation = 90;
% y-axis
yyaxis left
ax.YLim = [0,1.1];
ax.YColor = 'k';
ax.YTickLabelRotation = 90;
yyaxis right
ax.YLim = [0,1.1*maxpdf];
ax.YColor = 'k';
ax.YLabel.String = 'Probability'; ax.YLabel.Color = 'k';
ax.YTickLabelRotation = 90;
end

%% old version
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