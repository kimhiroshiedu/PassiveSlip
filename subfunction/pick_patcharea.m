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
patch(blon',blat',bdep')
hold on
geoshow([japan.Lat], [japan.Lon])
ax = gca;
ax.XLim = [130 142];
ax.YLim = [ 30  37];
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
