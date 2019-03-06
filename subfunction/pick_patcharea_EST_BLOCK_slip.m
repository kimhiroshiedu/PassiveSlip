fprintf('Now loading %s ...',[DIR,'/BLK.mat'])
load([DIR,'/BLK.mat']);fprintf('load\n')
fprintf('Now loading %s ...',[DIR,'/TCHA.mat'])
load([DIR,'/TCHA.mat']);fprintf('load\n')
% 
LO_Mc=0;
UP_Mc=1;
red=[0:1/32:1 ones(1,32)]';
green=[0:1/32:1 1-1/32:-1/32:0]';
blue=[ones(1,32) 1:-1/32:0]';
rwb=[red green blue];
rw =rwb(33:end,:);
if LO_Mc==-1
    cmap=rwb;
else
    cmap=rw;
end
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
NN=1;
for NB1=1:BLK(1).NBlock
    for NB2=NB1+1:BLK(1).NBlock
        NF=size(BLK(1).BOUND(NB1,NB2).blon,1);
        if NF~=0
            cp=TCHA.AVEFLT(NN:NN+NF-1,:);
            cp(cp<0.975)=0;
            patch(BLK(1).BOUND(NB1,NB2).blon',BLK(1).BOUND(NB1,NB2).blat',BLK(1).BOUND(NB1,NB2).bdep',cp);
            NN=NN+NF;
            hold on
        end
    end
end
ax=gca;
ax.CLim=[LO_Mc UP_Mc];
colormap(cmap)
colorbar
% 
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
plot(lon,lat,'g','LineWidth',5)
% 
figure(10); clf(10)
