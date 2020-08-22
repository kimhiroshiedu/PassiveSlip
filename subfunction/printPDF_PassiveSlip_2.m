function printPDF_PassiveSlip_2(blk,G,tcha,ud,T)
% tcha, blk and grn.mat must be loaded before using the function.

nud    = size(tcha.smpasp,1) / 2;
sint   = 10;
binint =  5;
W       = blk(1).aline_zd(ud) - blk(1).aline_zu(ud);
edges_d = blk(1).aline_zu(ud)    :binint:blk(1).aline_zd(ud)+2*W+binint;
edges_u = blk(1).aline_zu(ud)-2*W:binint:blk(1).aline_zd(ud)    +binint;

pdf_d = histcounts(tcha.smpasp(    ud,1:sint:end,T),edges_d,'Normalization','pdf');
pdf_u = histcounts(tcha.smpasp(nud+ud,1:sint:end,T),edges_u,'Normalization','pdf');
d_d = (edges_d(1:end-1) + edges_d(2:end)) ./ 2;
d_u = (edges_u(1:end-1) + edges_u(2:end)) ./ 2;
figure(230); clf(230)
npdf_d = pdf_d ./ max([pdf_d,pdf_u]);
npdf_u = pdf_u ./ max([pdf_d,pdf_u]);
stairs(d_d,npdf_d,'-b','LineWidth',1); hold on
stairs(d_u,npdf_u,'-r','LineWidth',1);
xline(tcha.aveasp(nud+ud,1,T),'--r');
xline(tcha.aveasp(    ud,1,T),'--b');
xlbot = xline(blk(1).aline_zd(ud),'-k');
xltop = xline(blk(1).aline_zu(ud),'-k');
xlbot.LineWidth = 1;
xltop.LineWidth = 1;
xlbot.Alpha = 1;
xltop.Alpha = 1;

ma = zeros(size(tcha.smpasp,1),1);
ma(ud) = 1;
idtri = G(1).zd * ma;
idtri = idtri == 1;
dep = G(1).zc(idtri);
pc  = tcha.aveaid(idtri,:,T);
x = blk(1).aline_zu(ud):1:blk(1).aline_zd(ud);
y = interp1(dep,pc,x,'linear','extrap');
plot(x,y,'Color','g','LineWidth',1)
view(90,90)

hold off

ax = gca;
ax.PlotBoxAspectRatio = [2,1,1];
ax.XLim = [blk(1).aline_zu(ud)-2*W,blk(1).aline_zd(ud)+2*W];
ax.YLim = [0,1.1];
ax.XColor = 'k';
ax.YColor = 'k';
ax.FontSize = 11;
ax.XLabel.String = 'Depth (km)';
ax.YLabel.String = 'PDF';
ax.XLabel.Color = 'k';
ax.YLabel.Color = 'k';

return
% saveas(1,['pdf_ud_',num2str(ud)])
% print(['pdf_ud_',num2str(ud)],'-dpdf','-painters')
end