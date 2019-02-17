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