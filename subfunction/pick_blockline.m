load('D:\MasterResearch\GMT\data\faults_jb1.txt')
load('D:\MasterResearch\inversion\PassiveSlip\Result_win\Test_09\blk.mat')

figure(10); clf(10)
% Coast line
latlim   = [ 10  50];
lonlim   = [120 150];
filename = gunzip('gshhs_l.b.gz', tempdir);
japan    = gshhs(filename{1},latlim,lonlim);
geoshow([japan.Lat], [japan.Lon])
ax = gca;
ax.XLim = [138 152]; % NE Japan
ax.YLim = [ 32  48]; % NE Japan
% ax.XLim = [130 142]; % SW Japan
% ax.YLim = [ 30  37]; % SW Japan

% JB1 block line
alon1 = faults_jb1(:,2);
alat1 = faults_jb1(:,3);
alon2 = faults_jb1(:,4);
alat2 = faults_jb1(:,5);
alon = [alon1; alon2];
alat = [alat1; alat2];
for nf = 1:size(faults_jb1,1)
  lon = faults_jb1(nf,[2,4]);
  lat = faults_jb1(nf,[3,5]);
  plot(lon,lat, 'k'); hold on
  plot(lon,lat,'xk'); hold on
end

% Kimura block line
blon = [];
blat = [];
for nb = 1:size(blk,2)
  blon = [blon; blk(nb).lon];
  blat = [blat; blk(nb).lat];
  plot(blk(nb).lon,blk(nb).lat, 'r'); hold on
  plot(blk(nb).lon,blk(nb).lat,'xr'); hold on
end

% Pick points
jlon = [];
jlat = [];
klon = [];
klat = [];
count = 0;
while 1
  [x,y,id1] = ginput(1);
  if id1 == 1
    count = count + 1;
    [dp,Ip] = min((alon-x).^2+(alat-y).^2);
    xc = alon(Ip);
    yc = alat(Ip);
    jlon = [jlon; [count, xc]];
    jlat = [jlat; [count, yc]];
    plot(xc,yc,'ob')
  elseif id1 == 3
    count = count + 1;
    [dp,Ip] = min((blon-x).^2+(blat-y).^2);
    xc = blon(Ip);
    yc = blat(Ip);
    klon = [klon; [count, xc]];
    klat = [klat; [count, yc]];
    plot(xc,yc,'sb')
  else
    break
  end
end

% Combine
newlon = [jlon; klon];
newlat = [jlat; klat];
[~,id] = sort(newlon(:,1));
newlon = newlon(id,:);
newlat = newlat(id,:);

% Combine2
Newlon = [];
Newlat = [];
Newlon = [Newlon; newlon];
Newlat = [Newlat; newlat];

% Final plot
plot(Newlon, Newlat, 'g')
