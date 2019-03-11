fid  = fopen('MODEL_OGnew/BLOCK_Int_sw_japan/triB_8_11.txt');
nf   = 0;
blon = zeros(1,3);
blat = zeros(1,3);
bdep = zeros(1,3);
while 1
  nf    = nf+1;
  loc_f = fscanf(fid,'%f %f %f \n', [3 3]);
  tline = fgetl(fid); if ~ischar(tline); break; end
  blon(nf,:) = loc_f(1,:);  % lon
  blat(nf,:) = loc_f(2,:);  % lat
  bdep(nf,:) = loc_f(3,:);  % hight
  tline = fgetl(fid); if ~ischar(tline); break; end
end
fclose(fid);
% 
latlim   = [ 10  40];
lonlim   = [120 150];
filename = gunzip('gshhs_i.b.gz', tempdir);
japan    = gshhs(filename{1},latlim,lonlim);
% 
figure(10); clf(10)
% patch(blon',blat',bdep')
patch(blon',blat','white')
hold on
geoshow([japan.Lat], [japan.Lon])
ax = gca;
ax.XLim = [130 142];
ax.YLim = [ 30  37];
% 
fid = fopen('PHS_plate/plate_nankai.txt');
np    = 0;
n     = 0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline) ; break; end
    if tline(1) ~= '>'
        n   = n+1;
        tmp = strsplit(strtrim(tline));
        contour.line(np+1).lon(n) = str2double(cellstr(tmp(1)));
        contour.line(np+1).lat(n) = str2double(cellstr(tmp(2)));
        contour.line(np+1).dep(n) = str2double(cellstr(tmp(3)));
    else
        np = np+1;
        n  = 0;
        continue;
    end
end
fclose(fid);
% 
hold on
for nl=1:size(contour.line,2)
    plot(contour.line(nl).lon,contour.line(nl).lat,'m')
    hold on
end
% 
lon=[];
lat=[];
while 1
  [x,y,id]=ginput(1);
  if id~=1; break; end
  lon = [lon,x];
  lat = [lat,y];
  hold on
  plot(x,y,'ro')
end
% 
lon=[lon, lon(1)];
lat=[lat, lat(1)];
hold on
plot(lon,lat,'r','LineWidth',5)
