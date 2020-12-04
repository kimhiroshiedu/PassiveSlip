%% Pole location relative to Nankai block
% PHS
phs.x = tcha.smppol(37,sid:end,1)-tcha.smppol(28,sid:end,1);
phs.y = tcha.smppol(38,sid:end,1)-tcha.smppol(29,sid:end,1);
phs.z = tcha.smppol(39,sid:end,1)-tcha.smppol(30,sid:end,1);
phs.lat = rad2deg(atan2(phs.z,sqrt(phs.x.^2+phs.y.^2)));
phs.lon = rad2deg(atan2(phs.y,phs.x));
phs.ang = rad2deg(sqrt(phs.x.^2+phs.y.^2+phs.z.^2).*1e6);

% IMP
imp.x = tcha.smppol(31,sid:end,1)-tcha.smppol(28,sid:end,1);
imp.y = tcha.smppol(32,sid:end,1)-tcha.smppol(29,sid:end,1);
imp.z = tcha.smppol(33,sid:end,1)-tcha.smppol(30,sid:end,1);
imp.lat = rad2deg(atan2(imp.z,sqrt(imp.x.^2+imp.y.^2)));
imp.lon = rad2deg(atan2(imp.y,imp.x));
imp.ang = rad2deg(sqrt(imp.x.^2+imp.y.^2+imp.z.^2).*1e6);

%% Augular histogram
figure(91); clf(91)
subplot(1,2,1)
h01 = histogram(phs.ang,'BinWidth',0.025,'Normalization','pdf','FaceAlpha',1,'FaceColor',[128,128,128]./255);
xl1 = xline(mean(phs.ang));
ax1 = gca;
ax1.PlotBoxAspectRatio = [1,0.6,1];
ax1.FontSize = 6;
ax1.FontName = 'Helvetica';
ax1.TickDir = "out";
ax1.XLabel.String = "(deg/myr)";
ax1.YLabel.String = "Frequency";
ax1.Title.String = "PHS-NAN";
xlim1 = ax1.XLim;
ylim1 = ax1.YLim;

subplot(1,2,2)
h02 = histogram(imp.ang,'BinWidth',0.025,'Normalization','pdf','FaceAlpha',1,'FaceColor',[128,128,128]./255);
xl2 = xline(mean(imp.ang));
ax2 = gca;
ax2.PlotBoxAspectRatio = [1,0.6,1];
ax2.FontSize = 6;
ax2.FontName = 'Helvetica';
ax2.TickDir = "out";
ax2.XLabel.String = "(deg/myr)";
ax2.YLabel.String = "Frequency";
ax2.Title.String = "IMP-NAN";
xlim2 = ax2.XLim;
ylim2 = ax2.YLim;

xlim = [min(xlim1(1),xlim2(1)),max(xlim1(2),xlim2(2))];
ylim = [min(ylim1(1),ylim2(1)),max(ylim1(2),ylim2(2))];
ax1.XLim = xlim;
ax1.YLim = ylim;
ax2.XLim = xlim;
ax2.YLim = ylim;

saveas(91,'angular_phs_nan')
print('angular_phs_nan','-dpdf','-painters')

%% 2D histogram of pole location
% PHS
figure(101); clf(101)
h1 = histogram2(phs.lon,phs.lat,'DisplayStyle',"tile",'EdgeColor',"none",'BinWidth',0.5);
xbin = (h1.XBinEdges(1:end-1)+h1.XBinEdges(2:end))./2;
ybin = (h1.YBinEdges(1:end-1)+h1.YBinEdges(2:end))./2;
[Xbin,Ybin] = meshgrid(xbin,ybin);
Xbin = Xbin';
Ybin = Ybin';
bw = h1.BinWidth;
Value = h1.Values(:);
meanx = mean(phs.x);
meany = mean(phs.y);
meanz = mean(phs.z);
meanlat = rad2deg(atan2(meanz,sqrt(meanx.^2+meany.^2)));
meanlon = rad2deg(atan2(meany,meanx));

fid = fopen('Euler_PAC.txt','wt');
fprintf(fid,'## lon  lat\n# %f %f\n',meanlon,meanlat);
for n = 1:size(Value,1)
    fprintf(fid,'> -Z %f\n',100*Value(n)./max(Value));
    fprintf(fid,'%f %f\n',[Xbin(n)-bw(1)/2,Ybin(n)-bw(2)/2]);
    fprintf(fid,'%f %f\n',[Xbin(n)+bw(1)/2,Ybin(n)-bw(2)/2]);
    fprintf(fid,'%f %f\n',[Xbin(n)+bw(1)/2,Ybin(n)+bw(2)/2]);
    fprintf(fid,'%f %f\n',[Xbin(n)-bw(1)/2,Ybin(n)+bw(2)/2]);
    fprintf(fid,'%f %f\n',[Xbin(n)-bw(1)/2,Ybin(n)-bw(2)/2]);
end
fclose(fid);

% IMP
figure(201); clf(102)
h2 = histogram2(imp.lon,imp.lat,'DisplayStyle',"tile",'EdgeColor',"none",'BinWidth',0.5);
xbin = (h2.XBinEdges(1:end-1)+h2.XBinEdges(2:end))./2;
ybin = (h2.YBinEdges(1:end-1)+h2.YBinEdges(2:end))./2;
[Xbin,Ybin] = meshgrid(xbin,ybin);
Xbin = Xbin';
Ybin = Ybin';
bw = h2.BinWidth;
Value = h2.Values(:);
meanx = mean(imp.x);
meany = mean(imp.y);
meanz = mean(imp.z);
meanlat = rad2deg(atan2(meanz,sqrt(meanx.^2+meany.^2)));
meanlon = rad2deg(atan2(meany,meanx));

fid = fopen('Euler_IMP.txt','wt');
fprintf(fid,'## lon  lat\n# %f %f\n',meanlon,meanlat);
for n = 1:size(Value,1)
    fprintf(fid,'> -Z %f\n',100*Value(n)./max(Value));
    fprintf(fid,'%f %f\n',[Xbin(n)-bw(1)/2,Ybin(n)-bw(2)/2]);
    fprintf(fid,'%f %f\n',[Xbin(n)+bw(1)/2,Ybin(n)-bw(2)/2]);
    fprintf(fid,'%f %f\n',[Xbin(n)+bw(1)/2,Ybin(n)+bw(2)/2]);
    fprintf(fid,'%f %f\n',[Xbin(n)-bw(1)/2,Ybin(n)+bw(2)/2]);
    fprintf(fid,'%f %f\n',[Xbin(n)-bw(1)/2,Ybin(n)-bw(2)/2]);
end
fclose(fid);

%% Relative motion at a point
% PHS
phs.gps = deg2rad([134,32.5]);
phs.gpsxy = [cos(phs.gps(2)) .* cos(phs.gps(1));...
             cos(phs.gps(2)) .* sin(phs.gps(1));...
             sin(phs.gps(2))];
phs.gpsxy = phs.gpsxy .* 6371e6;
gg = [           0   phs.gpsxy(3) -phs.gpsxy(2);...
      -phs.gpsxy(3)            0   phs.gpsxy(1);...
       phs.gpsxy(2) -phs.gpsxy(1)            0];
phs.v = gg * [phs.x;phs.y;phs.z];
vx = phs.v(1,:);
vy = phs.v(2,:);
vz = phs.v(3,:);
cn = cos(phs.gps(1));
sn = sin(phs.gps(1));
ct = cos(phs.gps(2));
st = sin(phs.gps(2));
vew=   -sn*vx   +cn*vy;
vns=-st*cn*vx-st*sn*vy+ct*vz;
vud= ct*cn*vx+ct*sn*vy+st*vz;
phs.venu = [vew;vns;vud];

figure(1); clf(1)
subplot(2,2,2)  % 2D histogram
h12 = histogram2(phs.venu(1,:),phs.venu(2,:),'BinWidth',0.5,'DisplayStyle',"tile",'EdgeColor',"none");
ax = gca;
ax.TickDir = "out";
ax.GridLineStyle = "none";
ax.GridLineStyle = ":";
ax.FontName = 'Helvetica';
ax.PlotBoxAspectRatio = [1,1,1];
ax.FontSize = 6;
xlim = ax.XLim;
ylim = ax.YLim;
colorbar
subplot(2,2,1)  % 1D histogram
h11 = histogram(phs.venu(2,:),'BinWidth',0.1,'FaceColor',[128,128,128]./255,'FaceAlpha',1);
ax = gca;
ax.PlotBoxAspectRatio = [1,0.6,1];
ax.FontName = 'Helvetica';
ax.FontSize = 6;
ax.TickDir = "out";
ax.XLim = ylim;
subplot(2,2,4)  % 1D histogram
h22 = histogram(phs.venu(1,:),'BinWidth',0.1,'FaceColor',[128,128,128]./255,'FaceAlpha',1);
ax = gca;
ax.PlotBoxAspectRatio = [1,0.6,1];
ax.FontName = 'Helvetica';
ax.FontSize = 6;
ax.TickDir = "out";
ax.XLim = xlim;
saveas(1,'vrel_phs-nan')
print('vrel_phs-nan','-dpdf','-painters')

% IMP
imp.gps = deg2rad([138,34]);
imp.gpsxy = [cos(imp.gps(2)) .* cos(imp.gps(1));...
             cos(imp.gps(2)) .* sin(imp.gps(1));...
             sin(imp.gps(2))];
imp.gpsxy = imp.gpsxy .* 6371e6;
gg = [           0   imp.gpsxy(3) -imp.gpsxy(2);...
      -imp.gpsxy(3)            0   imp.gpsxy(1);...
       imp.gpsxy(2) -imp.gpsxy(1)            0];
imp.v = gg * [imp.x;imp.y;imp.z];
vx = imp.v(1,:);
vy = imp.v(2,:);
vz = imp.v(3,:);
cn = cos(imp.gps(1));
sn = sin(imp.gps(1));
ct = cos(imp.gps(2));
st = sin(imp.gps(2));
vew=   -sn*vx   +cn*vy;
vns=-st*cn*vx-st*sn*vy+ct*vz;
vud= ct*cn*vx+ct*sn*vy+st*vz;
imp.venu = [vew;vns;vud];

figure(1); clf(1)
subplot(2,2,2)  % 2D histogram
h12 = histogram2(imp.venu(1,:),imp.venu(2,:),'BinWidth',0.5,'DisplayStyle',"tile",'EdgeColor',"none");
ax = gca;
ax.TickDir = "out";
ax.GridLineStyle = "none";
ax.GridLineStyle = ":";
ax.FontName = 'Helvetica';
ax.PlotBoxAspectRatio = [1,1,1];
ax.FontSize = 6;
xlim = ax.XLim;
ylim = ax.YLim;
colorbar
subplot(2,2,1)  % 1D histogram
h11 = histogram(imp.venu(2,:),'BinWidth',0.1,'FaceColor',[128,128,128]./255,'FaceAlpha',1);
ax = gca;
ax.PlotBoxAspectRatio = [1,0.6,1];
ax.FontName = 'Helvetica';
ax.FontSize = 6;
ax.TickDir = "out";
ax.XLim = ylim;
subplot(2,2,4)  % 1D histogram
h22 = histogram(imp.venu(1,:),'BinWidth',0.1,'FaceColor',[128,128,128]./255,'FaceAlpha',1);
ax = gca;
ax.PlotBoxAspectRatio = [1,0.6,1];
ax.FontName = 'Helvetica';
ax.FontSize = 6;
ax.TickDir = "out";
ax.XLim = xlim;
saveas(1,'vrel_imp-nan')
print('vrel_imp-nan','-dpdf','-painters')

