lon = [139,148];
lat = [46,33];

img = imread('D:/Users/hkimu/Downloads/asperity.png');
BW = im2bw(img);
Ired = img(:,:,1) == 255 & img(:,:,2) == 0 & img(:,:,3) == 0;
[x, y] = find(Ired);
dX = size(img,2);
dY = size(img,1);
dx = diff(lon);
dy = diff(lat);
ax = dX/dx;
ay = dY/dy;
ox = lon(1);
oy = lat(1);
X = y/ax + ox;
Y = x/ay + oy;
figure(10); clf(10)
image(lon,lat,img); hold on
plot(X,Y,'k.')
ax = gca;
ax.YDir = "normal";

% Pick points
alon = [];
alat = [];
slon = [];
slat = [];
count = 0;
while 1
  [x,y,id1] = ginput(1);
  if id1 == 3       % R-click
    count = count + 1;
    xc = x;
    yc = y;
    alon = [alon; [count, xc]];
    alat = [alat; [count, yc]];
    plot(xc,yc,'ob')
    %     text(xc,yc,[num2str(count)])
  elseif id1 == 1   % L-click
    count = count + 1;
    [dp,Ip] = min((X-x).^2+(Y-y).^2);
    xc = X(Ip);
    yc = Y(Ip);
    slon = [slon; [count, xc]];
    slat = [slat; [count, yc]];
    plot(xc,yc,'sb')
    %     text(xc,yc,[num2str(count)])
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

% Add start point
plon = [newlon; newlon(1,:)];
plat = [newlat; newlat(1,:)];

% Test plot
plot(plon(:,2),plat(:,2),'k')

% Save
fid = fopen('contour_1.txt','wt');
fprintf(fid,'%8.3f %7.3f\n',[plon(:,2),plat(:,2)]');
fclose(fid);