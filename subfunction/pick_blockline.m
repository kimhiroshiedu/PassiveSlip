% load('D:\MasterResearch\GMT\data\faults_jb1.txt')
% load('D:\MasterResearch\inversion\PassiveSlip\Result_win\Test_09\blk.mat')
load('~/MasterResearch/GMT/data/faults_jb1.txt')
load('~/MasterResearch/inversion/PassiveSlip/Result_red/Test_11/blk.mat')
load('~/MasterResearch/inversion/MCMC_Estslip/Result_red/Test_58/BLK.mat')

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
hold on

% % 1. JB1 block line
% alon1 = faults_jb1(:,2);
% alat1 = faults_jb1(:,3);
% alon2 = faults_jb1(:,4);
% alat2 = faults_jb1(:,5);
% alon = [alon1; alon2];
% alat = [alat1; alat2];
% for nf = 1:size(faults_jb1,1)
%   lon = faults_jb1(nf,[2,4]);
%   lat = faults_jb1(nf,[3,5]);
%   plot(lon,lat, 'k'); hold on
%   plot(lon,lat,'xk'); hold on
% end

% 1. Kimura block line (SW Japan)
load('~/MasterResearch/inversion/MCMC_Estslip/Result_red/Test_58/BLK.mat')
alon = [];
alat = [];
for nb = 1:size(blk,2)
alon = [alon; BLK(nb).LON];
alat = [alat; BLK(nb).LAT];
plot(BLK(nb).LON,BLK(nb).LAT, 'm'); hold on
plot(BLK(nb).LON,BLK(nb).LAT,'xm'); hold on
end

% 2. Kimura block line (NE Japan)
blon = [];
blat = [];
for nb = 1:size(blk,2)
  blon = [blon; blk(nb).lon];
  blat = [blat; blk(nb).lat];
  plot(blk(nb).lon,blk(nb).lat, 'c'); hold on
  plot(blk(nb).lon,blk(nb).lat,'xc'); hold on
end

% Pick points
jlon = [];
jlat = [];
klon = [];
klat = [];
flon = [];
flat = [];
count = 0;
while 1
  [x,y,id1] = ginput(1);
  if id1 == 1       % Mouse "L" click
    count = count + 1;
    [dp,Ip] = min((alon-x).^2+(alat-y).^2);
    xc = alon(Ip);
    yc = alat(Ip);
    jlon = [jlon; [count, xc]];
    jlat = [jlat; [count, yc]];
    plot(xc,yc,'ob')
  elseif id1 == 3   % Mouse "R" click
    count = count + 1;
    [dp,Ip] = min((blon-x).^2+(blat-y).^2);
    xc = blon(Ip);
    yc = blat(Ip);
    klon = [klon; [count, xc]];
    klat = [klat; [count, yc]];
    plot(xc,yc,'sb')
  elseif id1 == 32  % "space-key" type
    count = count + 1;
    xc = x;
    yc = y;
    flon = [flon; [count, xc]];
    flat = [flat; [count, yc]];
    plot(xc,yc,'sb')
  else
    break
  end
end

% Combine
newlon = [jlon; klon; flon];
newlat = [jlat; klat; flat];
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


function temp_block
load('MODEL_JP/BLOCK_ne_japan/02_Niigata.txt')
load('MODEL_JP/BLOCK_ne_japan/03_NE_Honshu_backarc.txt')
load('MODEL_JP/BLOCK_ne_japan/04_Okhotsk.txt')
load('MODEL_JP/BLOCK_ne_japan/05_Kuril_Forearc.txt')
load('MODEL_JP/BLOCK_ne_japan/06_NE_Honshu_forearc.txt')
load('MODEL_JP/BLOCK_ne_japan/01_Amur.txt')

plot(X05_Kuril_Forearc(:,1),X05_Kuril_Forearc(:,2), 'r')
plot(X05_Kuril_Forearc(:,1),X05_Kuril_Forearc(:,2),'xr')
plot(X06_NE_Honshu_forearc(:,1),X06_NE_Honshu_forearc(:,2), 'r')
plot(X06_NE_Honshu_forearc(:,1),X06_NE_Honshu_forearc(:,2),'xr')
plot(X03_NE_Honshu_backarc(:,1),X03_NE_Honshu_backarc(:,2), 'r')
plot(X03_NE_Honshu_backarc(:,1),X03_NE_Honshu_backarc(:,2),'xr')
plot(X02_Niigata(:,1),X02_Niigata(:,2), 'r')
plot(X02_Niigata(:,1),X02_Niigata(:,2),'xr')
plot(X04_Okhotsk(:,1),X04_Okhotsk(:,2), 'r')
plot(X04_Okhotsk(:,1),X04_Okhotsk(:,2),'xr')
plot(X01_Amur(:,1),X01_Amur(:,2), 'r')
plot(X01_Amur(:,1),X01_Amur(:,2),'xr')

blon = [blon; X05_Kuril_Forearc(:,1)];
blat = [blat; X05_Kuril_Forearc(:,2)];
blon = [blon; X06_NE_Honshu_forearc(:,1)];
blat = [blat; X06_NE_Honshu_forearc(:,2)];
blon = [blon; X03_NE_Honshu_backarc(:,1)];
blat = [blat; X03_NE_Honshu_backarc(:,2)];
blon = [blon; X02_Niigata(:,1)];
blat = [blat; X02_Niigata(:,2)];
blon = [blon; X04_Okhotsk(:,1)];
blat = [blat; X04_Okhotsk(:,2)];
blon = [blon; X01_Amur(:,1)];
blat = [blat; X01_Amur(:,2)];

end