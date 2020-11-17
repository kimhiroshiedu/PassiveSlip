function printPDF_PassiveSlip(folder,blk,G,tcha,udlist,T)
% tcha, blk and grn.mat must be loaded before using the function.

% Given parameters
sint   =  1;
binint =  2;
nud    =  size(tcha.smpasp,1) / 2;
sid    = (size(tcha.smpasp,2) / tcha.nit) * (tcha.burnin + 1);
prm.sint   = sint;
prm.binint = binint;
prm.nud    = nud;
prm.sid    = sid;

% Data prepare
n = 0;
maxpdf = 0;
for ud = udlist
  n = n + 1;
  flud = floor(ud);
  rat  = ud - flud;
  if rat == 0
    Zd = blk(1).aline_zd(ud);
    Zu = blk(1).aline_zu(ud);
    Xd = blk(1).aline_lond(ud);
    Xu = blk(1).aline_lonu(ud);
    Yd = blk(1).aline_latd(ud);
    Yu = blk(1).aline_latu(ud);
    z_d = tcha.smpasp(    ud,sid:sint:end,T); pick_z_d = z_d;
    z_u = tcha.smpasp(nud+ud,sid:sint:end,T); pick_z_u = z_u;    
  else
    Zd = blk(1).aline_zd(flud) + (blk(1).aline_zd(flud+1) - blk(1).aline_zd(flud)) * rat;
    Zu = blk(1).aline_zu(flud) + (blk(1).aline_zu(flud+1) - blk(1).aline_zu(flud)) * rat;
    Xd = blk(1).aline_lond(flud) + (blk(1).aline_lond(flud+1) - blk(1).aline_lond(flud)) * rat;
    Xu = blk(1).aline_lonu(flud) + (blk(1).aline_lonu(flud+1) - blk(1).aline_lonu(flud)) * rat;
    Yd = blk(1).aline_latd(flud) + (blk(1).aline_latd(flud+1) - blk(1).aline_latd(flud)) * rat;
    Yu = blk(1).aline_latu(flud) + (blk(1).aline_latu(flud+1) - blk(1).aline_latu(flud)) * rat;
    z_d = tcha.smpasp(    flud,sid:sint:end,T) + (tcha.smpasp(    flud+1,sid:sint:end,T)-tcha.smpasp(    flud,sid:sint:end,T)) .* rat; pick_z_d = z_d;
    z_u = tcha.smpasp(nud+flud,sid:sint:end,T) + (tcha.smpasp(nud+flud+1,sid:sint:end,T)-tcha.smpasp(nud+flud,sid:sint:end,T)) .* rat; pick_z_u = z_u;
  end
  W = Zd - Zu;
  udl(n).W = W;
  udl(n).Zd = Zd;
  udl(n).Zu = Zu;
  udl(n).Xd = Xd;
  udl(n).Xu = Xu;
  udl(n).Yd = Yd;
  udl(n).Yu = Yu;
  udl(n).edges_d = Zu    :binint:Zd+2*W+binint;
  udl(n).edges_u = Zu-2*W:binint:Zd    +binint;
  udl(n).pdf_d = histcounts(z_d,udl(n).edges_d,'Normalization','pdf');
  udl(n).pdf_u = histcounts(z_u,udl(n).edges_u,'Normalization','pdf');
  id_cross = z_d < z_u;
  pick_z_d(id_cross) = 99999;
  pick_z_u(id_cross) = 99999;
  udl(n).pdf_dcross = histcounts(pick_z_d,udl(n).edges_d,'Normalization','pdf');
  udl(n).pdf_ucross = histcounts(pick_z_u,udl(n).edges_u,'Normalization','pdf');
  if sum(~id_cross) == 0
    udl(n).pdf_dcross(:) = 0;
    udl(n).pdf_ucross(:) = 0;
  end
  udl(n).mean_pick_z_d = mean(pick_z_d);
  udl(n).mean_pick_z_u = mean(pick_z_u);
  maxpdf_ud = max([udl(n).pdf_d,udl(n).pdf_u]);
  if maxpdf_ud > maxpdf; maxpdf = maxpdf_ud; end
end

% Make, save figures
n = 0;
for ud = udlist
  n = n + 1;
  printpdf_ud(folder,blk,G,tcha,prm,udl,maxpdf,ud,n,T)
end

savelog(folder,udlist,udl)

end

function savelog(folder,udlist,udl)
savefolder = fullfile(pwd,folder,'figure','pdf_graphs');
if exist(savefolder,'dir') ~= 7; mkdir(savefolder); end
% Save log
fidlog = fopen(fullfile(savefolder,'makepdflog.txt'),'a+');
dt = datetime('now','TimeZone','Asia/Tokyo');
DateString = datestr(dt,'yyyy/mm/dd HH:MM:SS');
fprintf(fidlog,'%s\n ud_No.:',DateString);
for ud = udlist
  fprintf(fidlog,' %6.3f',ud);
end
fprintf(fidlog,'\n\n');
fclose(fidlog);
% Save interped udline
fidudl = fopen(fullfile(savefolder,'udline_interp.txt'),'wt');
fprintf(fidudl,'#    N      lon_d      lat_d      lon_u      lat_u\n');
n = 0;
for ud = udlist
  n = n + 1;
  fprintf(fidudl,'%6.3f %10.4f %10.4f %10.4f %10.4f\n',[ud,udl(n).Xd,udl(n).Yd,udl(n).Xu,udl(n).Yu]);
end
fclose(fidudl);
end

function printpdf_ud(folder,blk,G,tcha,prm,udl,maxpdf,ud,n,T)
flud = floor(ud);
rat  = ud - flud;
sint   = prm.sint;
binint = prm.binint;
nud    = prm.nud;
sid    = prm.sid;
W      = udl(n).W;
% Calc PDF of z_u and z_d
edges_d = udl(n).edges_d;
edges_u = udl(n).edges_u;
pdf_d = udl(n).pdf_d;
pdf_u = udl(n).pdf_u;
pdf_dcross = udl(n).pdf_dcross;
pdf_ucross = udl(n).pdf_ucross;
d_d = (edges_d(1:end-1) + edges_d(2:end)) ./ 2;
d_u = (edges_u(1:end-1) + edges_u(2:end)) ./ 2;
npdf_d = pdf_d ./ max([pdf_d,pdf_u]);
npdf_u = pdf_u ./ max([pdf_d,pdf_u]);

% Calc mean and 95% CI
if rat == 0
  zd_mean = tcha.aveasp(    ud,1,T);
  zu_mean = tcha.aveasp(nud+ud,1,T);
  sort_d = sort(tcha.smpasp(    ud,sid:sint:end,T));
  sort_u = sort(tcha.smpasp(nud+ud,sid:sint:end,T));
else
  z_d = tcha.smpasp(    flud,sid:sint:end,T) + (tcha.smpasp(    flud+1,sid:sint:end,T)-tcha.smpasp(    flud,sid:sint:end,T)) .* rat;
  z_u = tcha.smpasp(nud+flud,sid:sint:end,T) + (tcha.smpasp(nud+flud+1,sid:sint:end,T)-tcha.smpasp(nud+flud,sid:sint:end,T)) .* rat;
  zd_mean = mean(z_d);
  zu_mean = mean(z_u);
  sort_d = sort(z_d);
  sort_u = sort(z_u);
end
ci95_d = [sort_d(round(size(sort_d,2)*0.025)),sort_d(round(size(sort_d,2)*0.975))];
ci95_u = [sort_u(round(size(sort_u,2)*0.025)),sort_u(round(size(sort_u,2)*0.975))];

% Calc PDF of Pc
ma = zeros(size(tcha.smpasp,1),1);
ma(flud) = 1;
idtri = G(1).zd * ma;
if rat == 0
  idtri = idtri == 1;
else
  ma2 = zeros(size(tcha.smpasp,1),1);
  ma2(flud+1) = 1;
  idtri1 = round(idtri,3);
  idtri2 = round(G(1).zd * ma2,3);
  uqidtri1 = unique(idtri1);
  uqidtri2 = unique(idtri2);
  [M1,I1] = min(abs(1 - uqidtri1 - rat)); % minimum distance from left interped line
  [M2,I2] = min(abs(    uqidtri2 - rat)); % minimum distance from right interped line
  idtri = idtri1 == uqidtri1(I1) & idtri2 == uqidtri2(I2);
end
tri_dep = G(1).zc(idtri);
tri_pl  = tcha.aveaid(idtri,:,T);
pl_x = udl(n).Zu:binint:udl(n).Zd;
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
xl_bot = xline(udl(n).Zd,'-k','Linewidth',1,'Alpha',1);
xl_top = xline(udl(n).Zu,'-k','Linewidth',1,'Alpha',1);
ax1 = gca;
ax1.XLim = [udl(n).Zu-2*W,udl(n).Zd+2*W];
ax1.XLabel.String = 'Depth (km)';
ax1.XLabel.Rotation = 180;
ax1.XLabel.VerticalAlignment = "bottom";
ax1.YLabel.String = ['Line: ',num2str(ud)];
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
ax2.XLim = [udl(n).Zu,udl(n).Zd];
ax2 = graph_settings(ax2,maxpdf);
hold off

%% Save figures
savefolder = fullfile(pwd,folder,'figure','pdf_graphs');
if exist(savefolder,'dir') ~= 7; mkdir(savefolder); end
saveas(230,fullfile(savefolder,['T',num2str(T,'%02i'),'_pdf_ud_',num2str(floor(ud)),'-',num2str(round(ud-floor(ud),3)*1000,'%03i')]))
print(fullfile(savefolder,['T',num2str(T,'%02i'),'_pdf_ud_',num2str(floor(ud)),'-',num2str(round(ud-floor(ud),3)*1000,'%03i')]),'-dpdf','-painters')
print(fullfile(savefolder,['T',num2str(T,'%02i'),'_pdf_ud_',num2str(floor(ud)),'-',num2str(round(ud-floor(ud),3)*1000,'%03i')]),'-dpng')

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