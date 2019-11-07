% Trimesh
prm.dirblock = 'MODEL_JP/BLOCK_Int_ne_japan';
ext = 'triB_*.txt';
file = dir([prm.dirblock,'/',ext]);
[nbound,~] = size(file);
tri(1).nbound = nbound;
figure(10); clf(10)
for nb = 1:tri(1).nbound
  fid = fopen(fullfile(prm.dirblock,file(nb).name),'r');
  nf  = 0;
  tri(nb).blon = zeros(1,3);
  tri(nb).blat = zeros(1,3);
  tri(nb).bdep = zeros(1,3);
  while 1
    nf    = nf+1;
    loc_f = fscanf(fid,'%f %f %f \n', [3 3]);
    tline = fgetl(fid); if ~ischar(tline); break; end
    tri(nb).blon(nf,:) = loc_f(1,:);
    tri(nb).blat(nf,:) = loc_f(2,:);
    tri(nb).bdep(nf,:) = loc_f(3,:);
    tline = fgetl(fid); if ~ischar(tline); break; end
  end
  fclose(fid);
  patch(tri(nb).blon',tri(nb).blat',zeros(size(tri(nb).blon))'); hold on
  colormap('white')
end

% pre_tri_f = 'MODEL_JP/BLOCK_Int_ne_japan/triB_5_8.txt';
% fid = fopen(pre_tri_f,'r');
% nf   = 0;
% blon = zeros(1,3);
% blat = zeros(1,3);
% bdep = zeros(1,3);
% figure(10); clf(10)
% while 1
%   nf    = nf+1;
%   loc_f = fscanf(fid,'%f %f %f \n', [3 3]);
%   tline = fgetl(fid); if ~ischar(tline); break; end
%   blon(nf,:) = loc_f(1,:);
%   blat(nf,:) = loc_f(2,:);
%   bdep(nf,:) = loc_f(3,:);
%   tline = fgetl(fid); if ~ischar(tline); break; end
% end
% fclose(fid);
% tri.lon = blon;
% tri.lat = blat;
% tri.dep = bdep;
% patch(tri.lon',tri.lat',zeros(size(tri.lon))'); hold on
% colormap('white')

% Coast line
figure(10)
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

% Block line
prm.dirblock = 'MODEL_JP/BLOCK_ne_japan';
ext = '*.txt';
file = dir([prm.dirblock,'/',ext]);
[nblock,~] = size(file);
blk(1).nblock = nblock;
blon = [];
blat = [];
figure(10)
for nb = 1:blk(1).nblock
  tmp = load(fullfile(prm.dirblock,file(nb).name));
  blk(nb).name = file(nb).name;
  blk(nb).lon  = tmp(:,1);
  blk(nb).lat  = tmp(:,2);
  blon = [blon; blk(nb).lon];
  blat = [blat; blk(nb).lat];
  plot(blk(nb).lon,blk(nb).lat, 'r'); hold on
  plot(blk(nb).lon,blk(nb).lat,'xr'); hold on
end

% Pick points
alon = [];
alat = [];
slon = [];
slat = [];
count = 0;
while 1
  [x,y,id1] = ginput(1);
  if id1 == 1
    count = count + 1;
    xc = x;
    yc = y;
    alon = [alon; [count, xc]];
    alat = [alat; [count, yc]];
    plot(xc,yc,'om')
    text(xc,yc,[num2str(count)])
  elseif id1 == 3
    count = count + 1;
    [dp,Ip] = min((blon-x).^2+(blat-y).^2);
    xc = blon(Ip);
    yc = blat(Ip);
    slon = [slon; [count, xc]];
    slat = [slat; [count, yc]];
    plot(xc,yc,'sm')
    text(xc,yc,[num2str(count)])
  else
    break
  end
end

% Combine
newlon = [alon; slon];
newlat = [alat; slat];
[~,id] = sort(newlon(:,1));
newlon = newlon(id,:);
newlat = newlat(id,:);

% Combine2
Newlon = [];
Newlat = [];
Newlon = [Newlon; newlon];
Newlat = [Newlat; newlat];

% Final plot
figure(10)
plot(Newlon, Newlat, 'm')
