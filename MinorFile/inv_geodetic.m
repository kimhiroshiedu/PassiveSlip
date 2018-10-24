function inv_geodetic
%---------------------------------
Base_dir='/Users/take/Dropbox/yellow/kikai_jima';
%Base_dir='/home/take/Dropbox/yellow/kikai_jima';
LIST=fullfile(Base_dir,'data/list_all');
data_dir=fullfile(Base_dir,'data/pos');
fcomb=fullfile(Base_dir,'data/combination');
corr_F3=fullfile(Base_dir,'data/corrdat_F3.txt');
sub_f=fullfile(Base_dir,'data/nankai.dep');
bound_f=fullfile(Base_dir,'data/ryu.bo');
n_mesh=200;
%---------------------------------
%[GPS]=READ_OBS_LIST_F3(LIST,data_dir,corr_F3);
%save('save_data_all','GPS');
load('save_data_all');
%---------------------------------
[S]=mk_green_tri(GPS,sub_f,bound_f,n_mesh);
save('save_S','S');
[COMB,NCOMB]=SLOPE_DIS_CHANG_GPS(GPS,fcomb);
%---------------------------------
%{
st_num=datenum('1996/05/01','yyyy/mm/dd');
en_num=datenum('2013/01/01','yyyy/mm/dd');
figure
for N=1:NCOMB
  index=COMB(N).DATE>st_num & COMB(N).DATE<en_num;
  subplot(NCOMB./2,2,N);
  plot(COMB(N).DATE(index),100.*detrend(COMB(N).SDC(index)-mean(COMB(N).SDC(index))),'.')
  datetick('x','yyyy/mm','keeplimits')
  xlim([st_num en_num])
  xlabel('Date')
  ylabel('[cm]')
  title(strcat('Slope Distance Change ',COMB(N).trgID,'--',COMB(N).refID));
end
%}
%---------------------------------
FAULTS(1).lat=28.0;
FAULTS(1).lon=130.0;
FAULTS(1).str=45.0;
FAULTS(1).dep=20.0;
FAULTS(1).dip=20.0;
FAULTS(1).al=100.0;
FAULTS(1).aw=100.0;
FAULTS(1).rake=90.0;

[GREEN]=GREEN_SLOPE_DIS_CHANGE(COMB,GPS,FAULTS)

%---------------------------------
end
%====================================================
function [s]=mk_green_tri(GPS,sub_f,bound_f,n_mesh)
[s]=init_interface_tri(sub_f,bound_f,n_mesh.*10);
[s]=down_tri(s,GPS,n_mesh);
end
%====================================================
function [S]=down_tri(S,GPS,n_mesh)
alat0=(mean([GPS.lat])+mean([S.lat]))./2;
alon0=(mean([GPS.lon])+mean([S.lon]))./2;
[gx,gy]=PLTXY([GPS.lat],[GPS.lon],alat0,alon0);
[sx,sy]=PLTXY(S.lat,S.lon,alat0,alon0);
gz=[GPS.hight]./1000;
sz=[S.dep];
Ua=zeros(length([S.tri]),1);

tic
for n=1:length(Ua)
  [U]=CalcTriDisps(gx',gy',gz',sx(S.tri(n,:)),sy(S.tri(n,:)),sz(S.tri(n,:)),0.25,0,1,0);%
%  LENG(1)=sqrt((X(1)-X(2)).^2+(Y(1)-Y(2)).^2+(Z(1)-Z(2)).^2);
%  LENG(2)=sqrt((X(3)-X(2)).^2+(Y(3)-Y(2)).^2+(Z(3)-Z(2)).^2);
%  LENG(3)=sqrt((X(3)-X(4)).^2+(Y(3)-Y(4)).^2+(Z(3)-Z(4)).^2);
  Ua(n)=sum(sqrt(U.x.^2+U.y.^2+U.z.^2));
end
toc
Ntri=length([S.tri]);
while Ntri > n_mesh
  r_index=ones(length([S.lat]),1);
  [~,index]=min(Ua);
  min_tri=[S.tri(index,:)];
  f_tri=zeros(3,1);
  for n=1:3
    f_tri(n)=sum(find(S.tri,min_tri(n)));
  end
  [~,index]=max(f_tri);
  r_index(min_tri(index))=0;
  r_index=logical(r_index);
  S.lat=S.lat(r_index);
  S.lon=S.lon(r_index);
  S.dep=S.dep(r_index);
  [sx,sy]=PLTXY(S.lat,S.lon,alat0,alon0);
  sz=[S.dep];
%------------
  tri=delaunay(S.lon,S.lat);
  ntri=length(tri);
  nn=0;
  Stri=[];
  for n=1:ntri
    glon=mean(S.lon(tri(n,:)));  
    glat=mean(S.lat(tri(n,:)));
    ID=inpolygon(glon,glat,S.bound(:,1),S.bound(:,2));
    if ID==1
      nn=nn+1;
      Stri(nn,:)=tri(n,:);
    end  
  end
%------------
  Ua_tmp=Ua;
  Ua=zeros(nn,1);
  for n=1:nn
    if Stri(n,1)~=S.tri(n,1) || Stri(n,2)~=S.tri(n,2) || Stri(n,3)~=S.tri(n,3)
     [U]=CalcTriDisps(gx',gy',gz',sx(Stri(n,:)),sy(Stri(n,:)),sz(Stri(n,:)),0.25,0,1,0);
     Spx=sx(Stri(n,:));
     Spy=sy(Stri(n,:));
     Spz=sz(Stri(n,:));
     SS=(Spx(1)-Spx(2)).^2+(Spx(2)-Spx(3)).^2+(Spx(3)-Spx(1)).^2;
     
     Ua(n)=sum(sqrt(U.x.^2+U.y.^2+U.z.^2));
    else
     Ua(n)=Ua_tmp(n);
    end
  end
  S.tri=Stri;
  Ntri=length([S.tri]);
  Fid=figure;
  plot(S.bound(:,1),S.bound(:,2),'r');
  hold on;
  plot([GPS(:).lon],[GPS(:).lat],'.g');
  triplot(S.tri,S.lon,S.lat);
  title(['Number of triangels= ',num2str(Ntri)]);
  print(Fid,'-depsc ',['./figs/Mesh',num2str(Ntri)]);
  close(Fid)
end
Fid=figure;
plot(S.bound(:,1),S.bound(:,2),'r')
hold on
plot([GPS(:).lon],[GPS(:).lat],'.g')
triplot(S.tri,S.lon,S.lat);
title(['Number of triangels= ',num2str(Ntri)])
print(Fid,'-depsc ',['Mesh',num2str(Ntri)]);
pause(.1)
gz=[GPS.hight]./1000;
sz=[S.dep];
Ua=zeros(length([S.tri]),1);
for n=1:length(Ua)
  [U]=CalcTriDisps(gx',gy',gz',sx(S.tri(n,:)),sy(S.tri(n,:)),sz(S.tri(n,:)),0.25,0,1,0);
  S.Green_x(:,n)=U.x;
  S.Green_y(:,n)=U.y;
  S.Green_z(:,n)=U.z;
end
end
%====================================================
function [GREEN]=GREEN_SLOPE_DIS_CHANGE(COMB,GPS,FAULTS)
%
NFAU=length(FAULTS);
index_t=[COMB.trgNID];
index_r=[COMB.refNID];
GREEN=zeros(length(COMB),NFAU);
POIS=0.25;
for N=1:NFAU
  [Us,Ud]=DISLOC_RECT([GPS.lat]',[GPS.lon]',FAULTS(N).lat,...
          FAULTS(N).lon,FAULTS(N).dep,FAULTS(N).str,...
          FAULTS(N).dip,FAULTS(N).al,FAULTS(N).aw,POIS);
  GREEN(:,N)=cosd(FAULTS(N).rake).^2.*(sqrt((Us(index_t,1)-Us(index_r,1)).^2 + ...
             (Us(index_t,2)-Us(index_r,2)).^2+(Us(index_t,3)-Us(index_r,3)).^2)) + ...
             sind(FAULTS(N).rake).^2.*(sqrt((Ud(index_t,1)-Ud(index_r,1)).^2 + ...
             (Ud(index_t,2)-Ud(index_r,2)).^2+(Ud(index_t,3)-Ud(index_r,3)).^2));
end
end
%====================================================
function [COMB,NC]=SLOPE_DIS_CHANG_GPS(GPS,fcomb)
%---------------------------------
% slope distance change of combination list
%
% INPUT
% GPS       : GPS DATA (F3 format)
% fcomb     : combination list file
% 
% OUTPUT
% comb.trg      : Target  
% comb.ref      : Refrence  
% comb.???ID    : GPS ID 
% comb.???NID   : Array number for GPS ID 
% comb.???Index : Index for GPS 
% comb.DATE     : DATE for slope distance change
% comb.SDC      : slope distance change
% 
% Code by T.ITO 2012/06/13 ver0.1
%
Fid=fopen(fcomb,'r');
combin=textscan(Fid,'%s%s','delimiter',' ');
fclose(Fid);
NC=length(combin{1});
NGPS=length(GPS);
trgID=combin{1};
refID=combin{2};
for N=1:NC
  COMB(N).trgID=trgID(N);
  COMB(N).refID=refID(N);
  for ND=1:NGPS
    if strcmp(COMB(N).trgID,GPS(ND).ID); COMB(N).trgNID=ND; end
    if strcmp(COMB(N).refID,GPS(ND).ID); COMB(N).refNID=ND; end
  end
  [COMB(N).DATE,COMB(N).trgIndex,COMB(N).refIndex]=intersect(GPS(COMB(N).trgNID).DATE,GPS(COMB(N).refNID).DATE);
   COMB(N).SDC=sqrt((GPS(COMB(N).trgNID).XYZ(COMB(N).trgIndex,1)-GPS(COMB(N).refNID).XYZ(COMB(N).refIndex,1)).^2+...
                    (GPS(COMB(N).trgNID).XYZ(COMB(N).trgIndex,2)-GPS(COMB(N).refNID).XYZ(COMB(N).refIndex,2)).^2+...
                    (GPS(COMB(N).trgNID).XYZ(COMB(N).trgIndex,3)-GPS(COMB(N).refNID).XYZ(COMB(N).refIndex,3)).^2);
end
end
%====================================================
function [GPS]=READ_OBS_LIST_F3(LIST,data_dir,corr_F3)
%---------------------------------
% READ GEONET F3 Format from POS file list
%
% INPUT 
% LIST : File name of observation list
% data_dir : data dir
%
% Example
% LIST='/Users/take/Documents/MATLAB/data/filelist';
% data_dir='/Users/take/Documents/MATLAB/GEONET/';
% corr_F3='/Users/take/Documents/MATLAB/data/corrdat_F3.txt
%
% OUTPUT
% GPS       : GPS DATA (F3 format)
% GPS.XYZ   : XYZ
% GPS.DATE  : Observation date
% GPS.COUNT : Number of GPS data for a site
% GPS.ID    : ID Number
% GPS.lat   : latitude of GPS site
% GPS.lon   : longitude of GPS site
% GPS.hight : hight of GPS site
% see GEONET F3 format 
%
% You should update F3 of GEONET using wget
% http://192.168.15.35/~rsvdstaff/share/data/coordinates_F3/
% Code by T.ITO 2011/05/13 ver0.1
% Modified by T.ITO 2012/06/19 ver0.21
%
Fid_LIST=fopen(LIST,'r');
NGPS=0;
pre_ID=[];
while 1
  fname=fgetl(Fid_LIST);
  if ~ischar(fname); break; end
  [T_DATA]=READ_F3_FORMAT(fullfile(data_dir,fname));
  if strcmp(T_DATA.ID,pre_ID);
    GPS(NGPS).XYZ(COUNT+1:COUNT+T_DATA.COUNT,:)=T_DATA.XYZ;  
    GPS(NGPS).DATE(COUNT+1:COUNT+T_DATA.COUNT)=T_DATA.DATE;
    if GPS(NGPS).EPOCH_ST>T_DATA.EPOCH_ST; GPS(NGPS).EPOCH_ST=T_DATA.EPOCH_ST; end;
    if GPS(NGPS).EPOCH_ED<T_DATA.EPOCH_ED; GPS(NGPS).EPOCH_ED=T_DATA.EPOCH_ED; end;
    COUNT=COUNT+T_DATA.COUNT;
    GPS(NGPS).COUNT=COUNT;
  else
    NGPS=NGPS+1;
    GPS(NGPS)=T_DATA;
    COUNT=T_DATA.COUNT;
    pre_ID=T_DATA.ID;
  end;
end
fclose(Fid_LIST);
[GPS]=Corrdat_F3(GPS,corr_F3);
for N=1:NGPS
  tmp=xyz2blh(mean(GPS(N).XYZ));
  GPS(N).lat=tmp(1);
  GPS(N).lon=tmp(2);
  GPS(N).hight=tmp(3);
  [tmp,index]=sort(GPS(N).DATE);
  GPS(N).DATE=tmp;
  GPS(N).XYZ(:,1)=GPS(N).XYZ(index,1);
  GPS(N).XYZ(:,2)=GPS(N).XYZ(index,2);
  GPS(N).XYZ(:,3)=GPS(N).XYZ(index,3);
  GPS(N).NEU=single(xyz2neu([GPS(N).lat,GPS(N).lon],GPS(N).XYZ));  
end
end
%====================================================
function [DATA]=READ_F3_FORMAT(fname)
%
% READ F3 FORMAT OF GEONET
% Code by T.ITO 2011/05/03 ver0.1
% Modified by T.ITO 2012/06/13 ver0.2
%---------------------------------
%
fprintf('READING FILE %s \n',fname)
Fid=fopen(fname,'r');
%
DATA.ID=[]; DATA.J_NAME=[]; DATA.E_NAME=[];
DATA.DATE=[]; DATA.XYZ=[];
DATA.EPOCH_ST=[]; DATA.EPOCH_ED=[]; DATA.COUNT=[];
DATA.VERSION=[]; DATA.SOFT_NAME=[];
DATA.EPHEMERIS=[]; DATA.COORDINATE=[];
DATA.SOLUTION_ID=[];
%
while 1
tline=fgetl(Fid);
if isempty(tline);
elseif strcmp(tline(2:3),'ID'); DATA.ID=tline(15:end);
elseif strcmp(tline(1:5),'-DATA'); break;
elseif strcmp(tline(1:5),'+DATA');
  tline = fgetl(Fid); tline = fgetl(Fid);
  N=0;
  DATA.DATE=zeros(DATA.COUNT,1);
  DATA.XYZ=zeros(DATA.COUNT,3);
  while 1
    tline = fgetl(Fid);
    if strcmp(tline(1:4),'*---'); break; end
    N=N+1;
    DATA.DATE(N)=datenum(tline(2:19),'yyyy mm dd HH:MM:SS');
    DATA.XYZ(N,1)=str2double(tline(22:38));
    DATA.XYZ(N,2)=str2double(tline(40:56));
    DATA.XYZ(N,3)=str2double(tline(58:74));
  end
elseif strcmp(tline(2:6),'RINEX'); DATA.RINEX=tline(15:end);
elseif strcmp(tline(2:6),'EPOCH');
  DATA.EPOCH_ST=datenum(tline(21:42),'yyyy/mm/dd HH:MM:SS');
  DATA.EPOCH_ED=datenum(tline(46:64),'yyyy/mm/dd HH:MM:SS');
  DATA.COUNT=str2double(tline(73:end));
elseif strcmp(tline(2:7),'J_NAME'); DATA.J_NAME=tline(15:end);
elseif strcmp(tline(2:7),'E_NAME'); DATA.E_NAME=tline(15:end);
elseif strcmp(tline(2:8),'VERSION'); DATA.VERSION=tline(15:end);
elseif strcmp(tline(1:9),'+SITE/INF');
elseif strcmp(tline(1:9),'-SITE/INF');
elseif strcmp(tline(2:10),'SOFT_NAME'); DATA.SOFT_NAME=tline(15:end);
elseif strcmp(tline(2:10),'EPHEMERIS'); DATA.EPHEMERIS=tline(15:end);
elseif strcmp(tline(2:10),'ELLIPSOID'); DATA.ELLIPSOID=tline(15:end);
elseif strcmp(tline(2:11),'COORDINATE'); DATA.COORDINATE=tline(15:end);
elseif strcmp(tline(1:11),'+SOLVER/INF');
elseif strcmp(tline(1:11),'-SOLVER/INF');
elseif strcmp(tline(2:12),'SOLUTION_ID'); DATA.SOLUTION_ID=tline(15:end);
elseif strcmp(tline(1:12),'%*----+--+--');
else
  break
end
end
fclose(Fid);
end
%====================================================
function [GPS]=Corrdat_F3(GPS,corr_F3)
%
% You should be download corr_F3 
% http://mekira.gsi.go.jp/JAPANESE/corrf3o.dat
%
% Code by T.ITO 2012/06/13 ver0.1
%
Fid=fopen(corr_F3,'r');
cHead=0;
while 1
  tline=fgetl(Fid);
  if strcmp(tline(1),'#')
    cHead=cHead+1;
  else
    break;
  end
end
frewind(Fid);
for N=1:cHead;
  tline=fgetl(Fid);
end
Ngps=length(GPS);
corrdat=textscan(Fid,'%s%s%f%f%f%f%f%f%s','delimiter',',');
fclose(Fid);
Ncorr=length(corrdat{1,1});
for Nc=1:Ncorr
  for Ng=1:Ngps
    if strcmp(cell2mat(corrdat{1,1}(Nc)),GPS(Ng).ID);
      Gdate=datenum(corrdat{1,2}(Nc),'YYYY/mm/dd');
      index=find(GPS(Ng).DATE > Gdate);
      GPS(Ng).XYZ(index,1)=GPS(Ng).XYZ(index,1)+corrdat{1,3}(Nc);
      GPS(Ng).XYZ(index,2)=GPS(Ng).XYZ(index,2)+corrdat{1,4}(Nc);
      GPS(Ng).XYZ(index,3)=GPS(Ng).XYZ(index,3)+corrdat{1,5}(Nc);
      break
    end
  end
end
end
%====================================================
function dneu=xyz2neu(bl,dxyz)
%
% TRANSFORMATION FROM (DX,DY,DZ) => (DN,DE,DU)
% CODE BY T.ITO 2006/12/13 ver0.1
%
deg2rad=pi/180;
bl=bl(1:2).*deg2rad;
dneu=zeros(size(dxyz));
R=[-sin(bl(1)).*cos(bl(2)) -sin(bl(2)).*sin(bl(2)) cos(bl(1)) ; ...
   -sin(bl(2))              cos(bl(2))             0.0        ; ...
    cos(bl(1)).*cos(bl(2))  cos(bl(1)).*sin(bl(2)) sin(bl(1))];
for j=1:length(dxyz(:,1))
    dneu(j,:)=R*dxyz(j,:)';
end
end
%===================================================
function blh=xyz2blh(xyz)
%
% TRANSFORMATION FROM (X,Y,Z) => (LAT,LON,H)
% GRS80 = Geodetic Reference System 1980 
% Code by T.ITO 2006/12/13 ver0.1
% Bug fixed by T.ITO 2012/06/13 ver0.2
%
a=6378137.0;
f=1/298.257222101;
rad2deg=180/pi;
e2=f*(2-f);
r=sqrt(xyz(:,1).^2+xyz(:,2).^2);
lat0=atan2(xyz(:,3),r.*(1.0-e2));
h0=r.*cos(lat0);
diff=1;
while eps < diff
    dn=a./sqrt(1.0-e2.*sin(lat0).^2);
    h=r.*cos(lat0)+xyz(:,3).*sin(lat0)-(a.*a)./dn;   
    lat=atan2(xyz(:,3),r.*(1.0-e2.*dn./(dn+h)));
    diff=max(abs(lat-lat0)+abs(h-h0));
    lat0=lat;
    h0=h;
end
blh(:,1)=lat0.*rad2deg;
blh(:,2)=atan2(xyz(:,2),xyz(:,1)).*rad2deg;
blh(:,3)=h;
end
%===================================================
function Geodetic_fault_slip_MCMC
% Geodetic fault slip inversion using Monte Calro 
% For computaion using GPU 
% Need CalcTriDisps.m 
%
%
clear all
%
home_dir = pwd;
%
%====================================================
% INPUT FILES
%====================================================
cos_o_f=fullfile(home_dir,'obs_cos_jpl_sort.pos');
fault_f=fullfile(home_dir,'subfault4new.llz');
%====================================================
% READ FILES
%====================================================
[YYY,EEE,OST,OBS_NAME]=read_cos(cos_o_f);
[FLOC,DA,STR,DIP,NV,AL,AW,FLOC_GMT]=read_faults(fault_f);

pltmap([30,40],[130 140],FLOC(:,2),FLOC(:,1),'o','r');
pltmap([30,40],[130 140],OST(:,2),OST(:,1),'','b');

%====================================================
% MAKE GREEN FUNCTION
%====================================================
[GG]=MAKE_GREEN_FUNC_TRI(OST,FLOC_GMT);
GREEN_RECT;
end
%====================================================
function [GG]=MAKE_GREEN_FUNC_TRI(OST,FLOC,pr)
%
for N_O=1:length(FLOC)
  [sx,sy]=PLTXY(FLOC,FLOC,OST(NN,1),OST(NN,2));
  tic
  [U]=CalcTriDisps(sx,sy,sz,x,y,z,pr,ss,ts,ds);
  toc
end

% Arguments
%  sx : x-coordinates of observation points
%  sy : y-coordinates of observation points
%  sz : z-coordinates of observation points
%  x  : x-coordinates of triangle vertices.
%  y  : y-coordinates of triangle vertices.
%  z  : z-coordinates of triangle vertices.
%  pr : Poisson's ratio
%  ss : strike slip displacement
%  ts : tensile slip displacement
%  ds : dip slip displacement


end
%====================================================
function [FLOC,DA,STR,DIP,NV,AL,AW]=est_fault_p(loc_f)
[X,Y]=PLTXY(loc_f(:,2),loc_f(:,1),loc_f(1,2),loc_f(1,1));
[DA]=AREA(X,Y,loc_f(:,3));
[AL,AW]=f_leng(X,Y,loc_f(:,3),DA);
FLOC=mean(loc_f);
[STR,DIP,NV]=EST_STRDIP(X,Y,loc_f(:,3));
end
%====================================================
function [AL,AW]=f_leng(X,Y,Z,DA)
ALt(1)=sqrt((X(1)-X(2)).^2+(Y(1)-Y(2)).^2+(Z(1)-Z(2)).^2);
ALt(2)=sqrt((X(3)-X(4)).^2+(Y(3)-Y(4)).^2+(Z(3)-Z(4)).^2);
AL=(ALt(1)+ALt(2))./2;
AW=DA./AL;
AL=AL./2;
AW=AW./2;
end
%====================================================
function [DA]=AREA(X,Y,Z)
%==========
% CALC. AREA IN THREE DIMENSION USING HERON'S FOMULA
% CODE BY T.ITO (2006/3/4)
% Modified by T.ITO (2011/04/14)
%==========
LENG(1)=sqrt((X(1)-X(2)).^2+(Y(1)-Y(2)).^2+(Z(1)-Z(2)).^2);
LENG(2)=sqrt((X(3)-X(2)).^2+(Y(3)-Y(2)).^2+(Z(3)-Z(2)).^2);
LENG(3)=sqrt((X(3)-X(4)).^2+(Y(3)-Y(4)).^2+(Z(3)-Z(4)).^2);
LENG(4)=sqrt((X(1)-X(4)).^2+(Y(1)-Y(4)).^2+(Z(1)-Z(4)).^2);
LENG(5)=sqrt((X(2)-X(4)).^2+(Y(2)-Y(4)).^2+(Z(2)-Z(4)).^2);
S1=(LENG(1)+LENG(2)+LENG(5))./2;
S2=(LENG(3)+LENG(4)+LENG(5))./2;
DA1=sqrt(S1*(S1-LENG(1))*(S1-LENG(2))*(S1-LENG(5)));
DA2=sqrt(S2*(S2-LENG(3))*(S2-LENG(4))*(S2-LENG(5)));
DA=real(DA1+DA2);
end
%====================================================
function [STR,DIP,NV]=EST_STRDIP(X,Y,Z)
%==========
% CALC. STR AND DIP ON FAULT
% CODE BY T.ITO (2006/03/04)
% Modified by T.ITO (2007/01/12)
% Modified by T.ITO (2011/04/14)
% DEPTH IS MINUS
%==========
XV=[X(1)-X(2),Y(1)-Y(2),Z(1)-Z(2)];
YV=[X(3)-X(2),Y(3)-Y(2),Z(3)-Z(2)];
NV1=cross(XV,YV);
DP=[NV1(1),NV1(2),0];
ST=cross(DP,NV1);
DP=cross(NV1,ST);
STR1=atan2(ST(2),ST(1));
DIP1=atan2(DP(3),sqrt(DP(1).^2+DP(2).^2));
%
XV=[X(3)-X(4),Y(3)-Y(4),Z(3)-Z(4)];
YV=[X(1)-X(4),Y(1)-Y(4),Z(1)-Z(4)];
NV2=cross(XV,YV);
DP=[NV2(1),NV2(2),0];
ST=cross(DP,NV2);
DP=cross(NV2,ST);
STR2=atan2(ST(2),ST(1));
DIP2=atan2(DP(3),sqrt(DP(1).^2+DP(2).^2));
%
STR=(180./pi).*(STR1+STR2)./2; 
DIP=(180./pi).*(DIP1+DIP2)./2;
NV=(NV1+NV2)./2;
end
%====================================================
function [X,Y]=PLTXY(ALAT,ALON,ALAT0,ALON0)
%-------------------
%  PLTXY TRANSFORMS (ALAT,ALONG) TO (X,Y)
%  WHEN ICORD.NE.0  PLTXY MAKES NO CHANGE IN 
%  TRANSFORMATION BETWEEN (X,Y) AND (ALAT,ALONG).
%-------------------
A=6.378160e3;
E2=6.6944541e-3;
E12=6.7395719e-3;
D=5.72958e1;
RD=1.0/D;
RLAT = RD.*ALAT;
SLAT = sin(RLAT);
CLAT = cos(RLAT);
V2   = 1.0 + E12.*CLAT.^2;
AL   = ALON-ALON0;
PH1  = ALAT + (V2.*AL.^2.*SLAT.*CLAT)./(2.0*D);
RPH1 = PH1.*RD;
RPH2 = (PH1 + ALAT0).*0.5.*RD;
R    = A.*(1.0-E2)./sqrt((1.0-E2.*sin(RPH2).^2).^3);
AN   = A./sqrt(1.0-E2.*sin(RPH1).^2);
C1   = D./R;
C2   = D./AN;
Y    = (PH1-ALAT0)./C1;
X    = (AL.*CLAT)./C2+(AL.^3.*CLAT.*cos(2.0.*RLAT))./(6.0.*C2.*D.^2);
end
%====================================================
function [Mslip,Mlamda,coff_slip]=GPU_CALC_MCMC(GG,YYY,EEE,DA,Pr_Mw,PP,i_str,i_dip,i_lamda,home_dir,Nitr)
%
g=gpuDevice(1);
GPUMem=g.FreeMemory./4;
%
seed=1234;
rng(seed);              % Reset the CPU random number generator.
parallel.gpu.rng(seed); % Reset the GPU random number generator.
%
MaxRtime=50;
walk_dis_slip=0.5;
walk_dis_lamda=0.1;
%
Chain_length=128;
%
%%%%%%%%%%%%
% GPU Initialize 
%%%%%%%%%%%%
[NN,MM]=size(GG);
Max_itr=MaxRtime*Chain_length;
S_length=min([fix((GPUMem-(4*PP*(NN+2*MM)+MM*(NN+MM)))/(PP*MM+1)),Max_itr]);
%
G_slip_old=parallel.gpu.GPUArray.zeros(MM,PP,'single');
G_slip_old(1:2:MM,:)=repmat(i_str,1,PP);
G_slip_old(2:2:MM,:)=repmat(i_dip,1,PP);
G_GG=gpuArray(single(GG));
%
G_Mu=Mu.*parallel.gpu.GPUArray.ones(1,PP,'single');
%
YYM=repmat(YYY',1,PP);        YYM=single(YYM); G_YYM=gpuArray(YYM);
EEM=repmat(abs(1./EEE)',1,PP);EEM=single(EEM); G_EEM=gpuArray(EEM);
DAM=repmat(DA,1,PP).*1E6;     DAM=single(DAM); G_DAM=gpuArray(DAM);
%
G_C91=parallel.gpu.GPUArray.ones(1,PP,'single').*9.1;
G_C15=parallel.gpu.GPUArray.ones(1,PP,'single').*1.5;
%
G_chain_slip =parallel.gpu.GPUArray.zeros(S_length*PP,MM,'single');
G_chain_lamda=parallel.gpu.GPUArray.zeros(S_length*PP,1,'single');
G_lamda_old=parallel.gpu.GPUArray.zeros(PP,1,'single');
%
G_Re2w_old=parallel.gpu.GPUArray.inf(1,PP,'single');
G_dPr2_old=parallel.gpu.GPUArray.inf(1,PP,'single');
tRe2w=YYY*diag(abs(1./EEE))*YYY';
%
fprintf('Chain_length=%4d PP=%4d Save_length=%5d M=%3d \n',Chain_length,PP,S_length,MM);
%
%%%%%%%%%%%%
% GPU CALC
%%%%%%%%%%%%
Rtime=0;
iT=0;
%
while Rtime < MaxRtime
  tic;
  Nacc=0;
  Rtime=Rtime+1;
  for Chain=1:Chain_length
    G_slip=G_slip_old+walk_dis_slip.*(randn(MM,PP,'single')-0.5);
    G_Re=G_GG*G_slip-G_YYM;
    G_Re2w=sum(G_EEM.*G_Re.*G_Re,1);
%%%%%%%%%%%%
% Calculate Mw
%%%%%%%%%%%%
    G_Mw=(log10(G_Mu.*sum(G_DAM.*sqrt(G_slip(1:2:MM,:).^2+G_slip(2:2:MM,:).^2)))-G_C91)./G_C15;
%%%%%%%%%%%%
    if lamda_opt==1
      G_lamda=G_lamda_old+walk_dis_lamda.*(randn(1,PP,'single')-0.5);
      G_dPr2=(G_Mw-Pr_Mw).^2;
      G_tst=exp(-0.5*(G_Re2w-G_Re2w_old)+exp(G_lamda).*(G_dPr2-G_dPr2_old));
    else
      G_tst=exp(-0.5*(G_Re2w-G_Re2w_old));
      G_tst(tMw>G_Mw)=0;
    end      
%-----------
    index_M=G_tst>randn(PP,1);
    G_slip_old(:,index_M)=G_slip(:,index_M);
    G_Re2w_old(index_M)=G_Re2w(index_M);
%-----------
    iT=rem(iT,S_length)+1;
    G_chain_slip(:,(iT-1)*PP+1:iT*PP)=G_slip;
%-----------
    if lamda_opt==1
      G_lamda_old(index_M)=G_lamda(index_M);
      G_dPr2_old(index_M)=G_dPr2(index_M);
      G_chain_lamda((iT-1)*PP+1:iT*PP,1)=G_lamda(:);
    end
%-----------
    Nacc=Nacc+sum(index_M);
  end
  G_slip_old=repmat(mean(G_slip,1),PP);
  if lamda_opt==1
    G_lamda_old=repmat(mean(G_lamda),PP);
  end
  Re2w=gather(G_Re2w);
  maxRes=max(1-(Re2w./tRe2w));
  minRes=min(1-(Re2w./tRe2w));
  Racc=Nacc./(Chain_length.*PP);
  fprintf('T=%3d maxRes=%5.3f minRes=%5.3f Accept=%5.2f Time=%4.1fsec\n',Rtime,maxRes,minRes,100*Racc,toc)
end
%
chain_slip=int16(gather(G_chain_slip));
Mslip=mean(chain_slip,1);
Ofile=strcat('ITR_D1s_chain_',num2str(Nitr,'%02d'));
if lamda_opt==1
  chain_lamda=gather(G_chain_lamda);
  Mlamda=mean(chain_lamda);
  save_f=fullfile(home_dir,'dat_save_lamda',Ofile);
  save(save_f,'chain_slip','chain_lamda','-v7.3');
else
  save_f=fullfile(home_dir,'dat_save',Ofile);
  save(save_f,'chain_slip','-v7.3');
end
%%%%%%%%%%%%
% cofficient of corrlation 
%%%%%%%%%%%%
%
clear G_lamda G_lamda_old G_dPr2 G_dPr2_old G_tst G_chain_lamda
clear G_slip G_slip_old G_GG G_Mw G_Mu G_YYM G_EEM G_DAM G_C91 G_C15 G_Mw
%
G_chain_slip=G_chain_slip(1:5:S_length*PP,:);
Gxslip=bsxfun(@minus,G_chain_slip,sum(G_chain_slip,1)./S_length*PP);
%
Gcoff_slip = (Gxslip' * Gxslip) / (S_length*PP-1);
Gd_slip=sqrt(diag(Gcoff_slip));
Gcoff_slip=Gcoff_slip./(Gd_slip*Gd_slip');
coff_slip=gather(Gcoff_slip);
%
end
%====================================================
function [Us,Ud]=DISLOC_RECT(OLAT,OLON,ALAT0,ALON0,DEP,PHAI,DIP,AL,AW,POIS)
%
%	   E,N    : Coordinates of observation points (relative to fault centroid)
%	   DEP    : Depth of the fault centroid (DEP > 0)
%	   PHAI   : Strike-angle from North (in degrees)
%	   DIP    : Dip-angle (in degrees, must be scalar)
%	   AL     : Fault length in strike direction (LEN > 0)
%	   AW     : Fault width in dip direction (WIDTH > 0)
%	   ALP    : (LAMBDA+MYU)/(LAMBDA+2*MYU) 0.5
%      ALP    :  1.0/(2.0*(1.0-POIS)) ;
%
RAD=pi./180;
ALP=1.0/(2.0*(1.0-POIS));
nsite=length(OLAT);
%
%%PHAI= 90-67.1978097765380+180;
%PHAI= 0;
%%DIP=    9.1319772299803;
%%AL=     2.*7.2863138103036;
%%AW=     2.*7.6445662751718;
%%DEP=   17.55;
%%ALAT0= 38.18045;
%%ALON0=142.631;
%
SD=sin(DIP.*RAD);
CD=cos(DIP.*RAD);
%
THETA=(PHAI-90).*RAD;
ST=sin(THETA);
CT=cos(THETA);
%
CD(abs(CD)<=1.0e-3)=0;
CD(CD==0 && SD > 0)=1;
CD(CD==0 && SD < 0)=-1;
%
[Xobs,Yobs]=PLTXY(OLAT,OLON,ALAT0,ALON0);
%
X=(Xobs.*CT-Yobs.*ST);
Y=(Xobs.*ST+Yobs.*CT);
AL2=AL./2;
%
P=Y.*CD+DEP.*SD*ones(nsite,1);
Q=Y.*SD-DEP.*CD*ones(nsite,1);
%
% XI,ET,Q : FAULT COORDINATE
%           ORIGIN AT CENTER OF LOWER SIDE OF FAULT
%
Us=zeros(nsite,3);
Ud=zeros(nsite,3);
%
Uxs=zeros(nsite,1); Uys=zeros(nsite,1); Uzs=zeros(nsite,1);
Uxd=zeros(nsite,1); Uyd=zeros(nsite,1); Uzd=zeros(nsite,1);
%
[uxs uys uzs uxd uyd uzd]=STATIC(ALP,X+AL2,P,Q,SD,CD);
Uxs=Uxs+uxs; Uys=Uys+uys; Uzs=Uzs+uzs;
Uxd=Uxd+uxd; Uyd=Uyd+uyd; Uzd=Uzd+uzd;
%
[uxs uys uzs uxd uyd uzd]=STATIC(ALP,X-AL2,P,Q,SD,CD);
Uxs=Uxs-uxs; Uys=Uys-uys; Uzs=Uzs-uzs;
Uxd=Uxd-uxd; Uyd=Uyd-uyd; Uzd=Uzd-uzd;
%
[uxs uys uzs uxd uyd uzd]=STATIC(ALP,X+AL2,P-AW,Q,SD,CD);
Uxs=Uxs-uxs; Uys=Uys-uys; Uzs=Uzs-uzs;
Uxd=Uxd-uxd; Uyd=Uyd-uyd; Uzd=Uzd-uzd;
%
[uxs uys uzs uxd uyd uzd]=STATIC(ALP,X-AL2,P-AW,Q,SD,CD);
Uxs=Uxs+uxs; Uys=Uys+uys; Uzs=Uzs+uzs;
Uxd=Uxd+uxd; Uyd=Uyd+uyd; Uzd=Uzd+uzd;
% 
Us(:,1)= Uxs(:).*CT+Uys(:).*ST; %E
Us(:,2)=-Uxs(:).*ST+Uys(:).*CT; %N
Us(:,3)= Uzs(:);                %D
%
Ud(:,1)= Uxd(:).*CT+Uyd(:).*ST; %E
Ud(:,2)=-Uxd(:).*ST+Uyd(:).*CT; %N
Ud(:,3)= Uzd(:);                %D
end
%====================================================
function [Uxs Uys Uzs Uxd Uyd Uzd]=STATIC(ALP,XI,ET,Qc,SD,CD)
% Subroutine STATIC
% coded by Okada(1985)
%
XI2=XI.^2;
ET2=ET.^2;
Q2=Qc.^2;
R2=XI2+ET2+Q2;
R=sqrt(R2);
Dc=ET.*SD-Qc.*CD;
Yc=ET.*CD+Qc.*SD;
RET=R+ET;
%
RET(RET<0)=0;
%
nn=length(RET);
RE= zeros(nn,1);
DLE=zeros(nn,1);
TT= zeros(nn,1);
%
index=find(RET~=0);
RE(index)=1./RET(index);
DLE(index)=log(RET(index));
DLE(~index)=-log(R(~index)-ET(~index));
%
index=find(Qc~=0);
TT(index)=atan(XI(index).*ET(index)./Qc(index)./R(index));
%
RD=R+Dc;
%
if CD ~= 0; % FOR INCLINED FAULT
  TD=SD./CD;
  Xc=sqrt(XI2+Q2);
  A5=zeros(nn,1);
  index=find(XI~=0);
  A5(index)=2./CD.*atan((ET(index).*(Xc(index)+Qc(index).*CD)+Xc(index).*(R(index)+Xc(index)).*SD)./(XI(index).*(R(index)+Xc(index)).*CD));
  A4= 1./CD.*(log(RD)-SD.*DLE);
  A3= 1.*(Yc./RD./CD-DLE)+TD.*A4;
  A1=-1./CD.*XI./RD-TD.*A5;
elseif CD == 0; % FOR VERTICAL FAULT
  RD2=RD.*RD;
  A1=-1./2.*XI.*Qc./RD2;
  A3= 1./2.*(ET./RD+Yc.*Qc./RD2-DLE);
  A4=-1.*Qc./RD;
  A5=-1.*XI.*SD./RD;
end
% COMMON
A2=-1.*DLE-A3;
UN=1/pi/2;
%
% STRIKE-SLIP CONTRIBUTION
REQ=RE./R.*Qc;
Uxs=-UN.*(REQ.*XI+TT+ALP.*A1.*SD);
Uys=-UN.*(REQ.*Yc+Qc*CD.*RE+ALP.*A2.*SD);
Uzs=-UN.*(REQ.*Dc+Qc*SD.*RE+ALP.*A4.*SD);
%
% DIP-SLIP CONTRIBUTION
SDCD=SD.*CD;
RRX=1./(R.*(R+XI));
Uxd=-UN.*(Qc./R-ALP.*A3.*SDCD);
Uyd=-UN.*(Yc.*Qc.*RRX+CD.*TT-ALP.*A1.*SDCD);
Uzd=-UN.*(Dc.*Qc.*RRX+SD.*TT-ALP.*A5.*SDCD);
end
%====================================================
function [s]=init_interface_tri(sub_f,bound_f,int_mesh)
%====================================================
Fid=fopen(sub_f);
dep_sub=textscan(Fid,'%f%f%f');
fclose(Fid);
dep_sub=cell2mat(dep_sub);
%====================================================
Fid=fopen(bound_f);
bound=textscan(Fid,'%f%f%f');
fclose(Fid);
bound=cell2mat(bound);
%====================================================
F=TriScatteredInterp(dep_sub(:,1),dep_sub(:,2),dep_sub(:,3),'natural');
min_lon=min(bound(:,1)); max_lon=max(bound(:,1));
min_lat=min(bound(:,2)); max_lat=max(bound(:,2));
figure
plot(bound(:,1),bound(:,2),'r')
hold on
n=0;
while n<int_mesh
  slat=(max_lat-min_lat).*rand(1)+min_lat;
  slon=(max_lon-min_lon).*rand(1)+min_lon;
  ID=inpolygon(slon,slat,bound(:,1),bound(:,2));
  if ID==1
    n=n+1;
    s.lat(n)=slat;
    s.lon(n)=slon;
    s.dep(n)=F(slon,slat);
    if rem(n,round(int_mesh/10))==1;
      plot3(s.lon,s.lat,s.dep,'.')
      pause(.1)
    end
  end
end
plot3(s.lon,s.lat,s.dep,'.')
%====================================================
tri = delaunay(s.lon,s.lat);
%====================================================
ntri=length(tri);
nn=0;
for n=1:ntri
  glon=mean(s.lon(tri(n,:)));  
  glat=mean(s.lat(tri(n,:)));
  ID=inpolygon(glon,glat,bound(:,1),bound(:,2));
  if ID==1
    nn=nn+1;
    s.tri(nn,:)=tri(n,:);
  end  
end
figure
plot(bound(:,1),bound(:,2),'r')
hold on
triplot(s.tri,s.lon,s.lat);
s.bound=bound;
s.dep_sub=dep_sub;
end
%====================================================
function [U] = CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds)
% CalcTriDisps.m
%
% Calculates displacements due to slip on a triangular dislocation in an
% elastic half space utilizing the Comninou and Dunders (1975) expressions
% for the displacements due to an angular dislocation in an elastic half
% space.
%
% Arguments
%  sx : x-coordinates of observation points
%  sy : y-coordinates of observation points
%  sz : z-coordinates of observation points
%  x  : x-coordinates of triangle vertices.
%  y  : y-coordinates of triangle vertices.
%  z  : z-coordinates of triangle vertices.
%  pr : Poisson's ratio
%  ss : strike slip displacement
%  ts : tensile slip displacement
%  ds : dip slip displacement
%
% Returns
%  U  : structure containing the displacements (U.x, U.y, U.z)
%
% This paper should and related code should be cited as:
% Brendan J. Meade, Algorithms for the calculation of exact 
% displacements, strains, and stresses for Triangular Dislocation 
% Elements in a uniform elastic half space, Computers & 
% Geosciences (2007), doi:10.1016/j.cageo.2006.12.003.
%
% Use at your own risk and please let me know of any bugs/errors!
%
% Copyright (c) 2006 Brendan Meade
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
%
% Calculate the slip vector in XYZ coordinates
normVec                      = cross([x(2);y(2);z(2)]-[x(1);y(1);z(1)], [x(3);y(3);z(3)]-[x(1);y(1);z(1)]);
normVec                      = normVec./norm(normVec);
if (normVec(3) < 0) % Enforce clockwise circulation
   normVec                   = -normVec;
   [x(2) x(3)]               = swap(x(2), x(3));
   [y(2) y(3)]               = swap(y(2), y(3));
   [z(2) z(3)]               = swap(z(2), z(3));
end
strikeVec                    = [-sin(atan2(normVec(2),normVec(1))) cos(atan2(normVec(2),normVec(1))) 0];
dipVec                       = cross(normVec, strikeVec);
slipComp                     = [ss ds ts];
slipVec                      = [strikeVec(:) dipVec(:) normVec(:)] * slipComp(:);
% Solution vectors
U.x                          = zeros(size(sx));
U.y                          = zeros(size(sx));
U.z                          = zeros(size(sx));
% Add a copy of the first vertex to the vertex list for indexing
x(4)                         = x(1);
y(4)                         = y(1);
z(4)                         = z(1);
%
for iTri = 1:3
   % Calculate strike and dip of current leg
   strike                   = 180/pi*(atan2(y(iTri+1)-y(iTri), x(iTri+1)-x(iTri)));
   segMapLength             = sqrt((x(iTri)-x(iTri+1))^2 + (y(iTri)-y(iTri+1))^2);
   [rx ry]                  = RotateXyVec(x(iTri+1)-x(iTri), y(iTri+1)-y(iTri), -strike);
   dip                      = 180/pi*(atan2(z(iTri+1)-z(iTri), rx));
%   
   if dip >= 0
      beta                  = pi/180*(90-dip);
      if beta > pi/2
         beta               = pi/2-beta;
      end
   else
      beta                  = -pi/180*(90+dip);
      if beta < -pi/2
         beta = pi/2-abs(beta);
      end
   end
%
   ssVec                    = [cos(strike/180*pi) sin(strike/180*pi) 0];
   tsVec                    = [-sin(strike/180*pi) cos(strike/180*pi) 0];
   dsVec                    = cross(ssVec, tsVec);
   lss                      = dot(slipVec, ssVec);
   lts                      = dot(slipVec, tsVec);
   lds                      = dot(slipVec, dsVec);
%
   if (abs(beta) > 0.000001) && (abs(beta-pi) > 0.000001)
      % First angular dislocation
      [sx1 sy1]                 = RotateXyVec(sx-x(iTri), sy-y(iTri), -strike);
      [ux1 uy1 uz1]             = adv(sx1, sy1, sz-z(iTri), z(iTri), beta, pr, lss, lts, lds);
      [ux1_o uy1_o uz1_o]       = adv_opt(sx1, sy1, sz-z(iTri), z(iTri), beta, pr, lss, lts, lds);
       sum(ux1-ux1_o+uy1-uy1_o+uz1-uz1_o)                            
      % Second angular dislocation
      [sx2 sy2]                 = RotateXyVec(sx-x(iTri+1), sy-y(iTri+1), -strike); 
      [ux2 uy2 uz2]             = adv_opt(sx2, sy2, sz-z(iTri+1), z(iTri+1), beta, pr, lss, lts, lds);

      % Rotate vectors to correct for strike
      [uxn uyn]                 = RotateXyVec(ux1-ux2, uy1-uy2, strike);
      uzn                       = uz1-uz2;
 
      % Add the displacements from current leg
      U.x                       = U.x + uxn;
      U.y                       = U.y + uyn;
      U.z                       = U.z + uzn;
   end
end

% Identify indices for stations under current triangle
inPolyIdx                       = find(inpolygon(sx, sy, x, y) == 1);
underIdx = [];
for iIdx = 1 : numel(inPolyIdx)
   d                            = LinePlaneIntersect(x, y, z, sx(inPolyIdx(iIdx)), sy(inPolyIdx(iIdx)), sz(inPolyIdx(iIdx)));
   if d(3)-sz(inPolyIdx(iIdx)) < 0
      underIdx = [underIdx ; inPolyIdx(iIdx)];
   end
end
% Apply static offset to the points that lie underneath the current triangle
U.x(underIdx)                = U.x(underIdx) - slipVec(1);
U.y(underIdx)                = U.y(underIdx) - slipVec(2);
U.z(underIdx)                = U.z(underIdx) - slipVec(3);
%
end
%====================================================
function d = LinePlaneIntersect(x, y, z, sx, sy, sz)
% Calculate the intersection of a line and a plane using a parametric
% representation of the plane.  This is hardcoded for a vertical line.
numerator                       = [1 1 1 1 ; x(1) x(2) x(3) sx ; y(1) y(2) y(3) sy ; z(1) z(2) z(3) sz];
numerator                       = det(numerator);
denominator                     = [1 1 1 0 ; x(1) x(2) x(3) 0 ; y(1) y(2) y(3) 0 ; z(1) z(2) z(3) -sz];
denominator                     = det(denominator);
if denominator == 0;
   denominator                  = eps;
end
t                               = numerator/denominator; % parametric curve parameter
d                               = [sx sy sz]-([sx sy 0]-[sx sy sz])*t;
end
%====================================================
function [a b] = swap(a, b)
% Swap two values
temp                            = a;
a                               = b;
b                               = temp;
end
%====================================================
function [xp yp] = RotateXyVec(x, y, alpha)
% Rotate a vector by an angle alpha
x                             = x(:);
y                             = y(:);
alpha                         = pi/180*alpha;
xp                            = cos(alpha).*x - sin(alpha).*y;
yp                            = sin(alpha).*x + cos(alpha).*y;
end
%====================================================
function [v1 v2 v3] = adv(y1, y2, y3, a, beta, nu, B1, B2, B3)
% These are the displacements in a uniform elastic half space due to slip
% on an angular dislocation (Comninou and Dunders, 1975).  Some of the
% equations for the B2 and B3 cases have been corrected following Thomas
% 1993.  The equations are coded in way such that they roughly correspond
% to each line in original text.  Exceptions have been made where it made 
% more sense because of grouping symbols.

sinbeta           = sin(beta);
cosbeta           = cos(beta);
cotbeta           = cot(beta);
z1                = y1.*cosbeta - y3.*sinbeta;
z3                = y1.*sinbeta + y3.*cosbeta;
R2                = y1.*y1 + y2.*y2 + y3.*y3;
R                 = sqrt(R2);
y3bar             = y3 + 2.*a;
z1bar             = y1.*cosbeta + y3bar.*sinbeta;
z3bar             = -y1.*sinbeta + y3bar.*cosbeta;
R2bar             = y1.*y1 + y2.*y2 + y3bar.*y3bar;
Rbar              = sqrt(R2bar);
F                 = -atan2(y2, y1) + atan2(y2, z1) + atan2(y2.*R.*sinbeta, y1.*z1+(y2.*y2).*cosbeta);
Fbar              = -atan2(y2, y1) + atan2(y2, z1bar) + atan2(y2.*Rbar.*sinbeta, y1.*z1bar+(y2.*y2).*cosbeta);

% Case I: Burgers vector (B1,0,0)
v1InfB1           = 2.*(1-nu).*(F+Fbar) - y1.*y2.*(1./(R.*(R-y3)) + 1./(Rbar.*(Rbar+y3bar))) - ...
                    y2.*cosbeta.*((R.*sinbeta-y1)./(R.*(R-z3)) + (Rbar.*sinbeta-y1)./(Rbar.*(Rbar+z3bar)));
v2InfB1           = (1-2.*nu).*(log(R-y3)+log(Rbar+y3bar) - cosbeta.*(log(R-z3)+log(Rbar+z3bar))) - ...
                    y2.*y2.*(1./(R.*(R-y3))+1./(Rbar.*(Rbar+y3bar)) - cosbeta.*(1./(R.*(R-z3))+1./(Rbar.*(Rbar+z3bar))));
v3InfB1           = y2 .* (1./R - 1./Rbar - cosbeta.*((R.*cosbeta-y3)./(R.*(R-z3)) - (Rbar.*cosbeta+y3bar)./(Rbar.*(Rbar+z3bar))));
v1InfB1           = v1InfB1 ./ (8.*pi.*(1-nu));
v2InfB1           = v2InfB1 ./ (8.*pi.*(1-nu));
v3InfB1           = v3InfB1 ./ (8.*pi.*(1-nu));

v1CB1             = -2.*(1-nu).*(1-2.*nu).*Fbar.*(cotbeta.*cotbeta) + (1-2.*nu).*y2./(Rbar+y3bar) .* ((1-2.*nu-a./Rbar).*cotbeta - y1./(Rbar+y3bar).*(nu+a./Rbar)) + ...
                    (1-2.*nu).*y2.*cosbeta.*cotbeta./(Rbar+z3bar).*(cosbeta+a./Rbar) + a.*y2.*(y3bar-a).*cotbeta./(Rbar.*Rbar.*Rbar) + ...
                    y2.*(y3bar-a)./(Rbar.*(Rbar+y3bar)).*(-(1-2.*nu).*cotbeta + y1./(Rbar+y3bar) .* (2.*nu+a./Rbar) + a.*y1./(Rbar.*Rbar)) + ...
                    y2.*(y3bar-a)./(Rbar.*(Rbar+z3bar)).*(cosbeta./(Rbar+z3bar).*((Rbar.*cosbeta+y3bar) .* ((1-2.*nu).*cosbeta-a./Rbar).*cotbeta + 2.*(1-nu).*(Rbar.*sinbeta-y1).*cosbeta) - a.*y3bar.*cosbeta.*cotbeta./(Rbar.*Rbar));
v2CB1             = (1-2.*nu).*((2.*(1-nu).*(cotbeta.*cotbeta)-nu).*log(Rbar+y3bar) -(2.*(1-nu).*(cotbeta.*cotbeta)+1-2.*nu).*cosbeta.*log(Rbar+z3bar)) - ...
                    (1-2.*nu)./(Rbar+y3bar).*(y1.*cotbeta.*(1-2.*nu-a./Rbar) + nu.*y3bar - a + (y2.*y2)./(Rbar+y3bar).*(nu+a./Rbar)) - ...
                    (1-2.*nu).*z1bar.*cotbeta./(Rbar+z3bar).*(cosbeta+a./Rbar) - a.*y1.*(y3bar-a).*cotbeta./(Rbar.*Rbar.*Rbar) + ...
                    (y3bar-a)./(Rbar+y3bar).*(-2.*nu + 1./Rbar.*((1-2.*nu).*y1.*cotbeta-a) + (y2.*y2)./(Rbar.*(Rbar+y3bar)).*(2.*nu+a./Rbar)+a.*(y2.*y2)./(Rbar.*Rbar.*Rbar)) + ...
                    (y3bar-a)./(Rbar+z3bar).*((cosbeta.*cosbeta) - 1./Rbar.*((1-2.*nu).*z1bar.*cotbeta+a.*cosbeta) + a.*y3bar.*z1bar.*cotbeta./(Rbar.*Rbar.*Rbar) - 1./(Rbar.*(Rbar+z3bar)) .* ((y2.*y2).*(cosbeta.*cosbeta) - a.*z1bar.*cotbeta./Rbar.*(Rbar.*cosbeta+y3bar)));

v3CB1             = 2.*(1-nu).*(((1-2.*nu).*Fbar.*cotbeta) + (y2./(Rbar+y3bar).*(2.*nu+a./Rbar)) - (y2.*cosbeta./(Rbar+z3bar).*(cosbeta+a./Rbar))) + ...
                    y2.*(y3bar-a)./Rbar.*(2.*nu./(Rbar+y3bar)+a./(Rbar.*Rbar)) + ...
                    y2.*(y3bar-a).*cosbeta./(Rbar.*(Rbar+z3bar)).*(1-2.*nu-(Rbar.*cosbeta+y3bar)./(Rbar+z3bar).*(cosbeta + a./Rbar) - a.*y3bar./(Rbar.*Rbar));

v1CB1             = v1CB1 ./ (4.*pi.*(1-nu));
v2CB1             = v2CB1 ./ (4.*pi.*(1-nu));
v3CB1             = v3CB1 ./ (4.*pi.*(1-nu));

v1B1              = v1InfB1 + v1CB1;
v2B1              = v2InfB1 + v2CB1;
v3B1              = v3InfB1 + v3CB1;


% Case II: Burgers vector (0,B2,0)
v1InfB2           = -(1-2.*nu).*(log(R-y3) + log(Rbar+y3bar)-cosbeta.*(log(R-z3)+log(Rbar+z3bar))) + ...
                    y1.*y1.*(1./(R.*(R-y3))+1./(Rbar.*(Rbar+y3bar))) + z1.*(R.*sinbeta-y1)./(R.*(R-z3)) + z1bar.*(Rbar.*sinbeta-y1)./(Rbar.*(Rbar+z3bar));
v2InfB2           = 2.*(1-nu).*(F+Fbar) + y1.*y2.*(1./(R.*(R-y3))+1./(Rbar.*(Rbar+y3bar))) - y2.*(z1./(R.*(R-z3))+z1bar./(Rbar.*(Rbar+z3bar)));
v3InfB2           = -(1-2.*nu).*sinbeta.*(log(R-z3)-log(Rbar+z3bar)) - y1.*(1./R-1./Rbar) + z1.*(R.*cosbeta-y3)./(R.*(R-z3)) - z1bar.*(Rbar.*cosbeta+y3bar)./(Rbar.*(Rbar+z3bar));
v1InfB2           = v1InfB2 ./ (8.*pi.*(1-nu));
v2InfB2           = v2InfB2 ./ (8.*pi.*(1-nu));
v3InfB2           = v3InfB2 ./ (8.*pi.*(1-nu));

v1CB2             = (1-2.*nu).*((2.*(1-nu).*(cotbeta.*cotbeta)+nu).*log(Rbar+y3bar) - (2.*(1-nu).*(cotbeta.*cotbeta)+1).*cosbeta.*log(Rbar+z3bar)) + ...
                    (1-2.*nu)./(Rbar+y3bar).* (-(1-2.*nu).*y1.*cotbeta+nu.*y3bar-a+a.*y1.*cotbeta./Rbar + (y1.*y1)./(Rbar+y3bar).*(nu+a./Rbar)) - ...
                    (1-2.*nu).*cotbeta./(Rbar+z3bar).*(z1bar.*cosbeta - a.*(Rbar.*sinbeta-y1)./(Rbar.*cosbeta)) - a.*y1.*(y3bar-a).*cotbeta./(Rbar.*Rbar.*Rbar) + ...
                    (y3bar-a)./(Rbar+y3bar).*(2.*nu + 1./Rbar.*((1-2.*nu).*y1.*cotbeta+a) - (y1.*y1)./(Rbar.*(Rbar+y3bar)).*(2.*nu+a./Rbar) - a.*(y1.*y1)./(Rbar.*Rbar.*Rbar)) + ...
                    (y3bar-a).*cotbeta./(Rbar+z3bar).*(-cosbeta.*sinbeta+a.*y1.*y3bar./(Rbar.*Rbar.*Rbar.*cosbeta) + (Rbar.*sinbeta-y1)./Rbar.*(2.*(1-nu).*cosbeta - (Rbar.*cosbeta+y3bar)./(Rbar+z3bar).*(1+a./(Rbar.*cosbeta))));
v2CB2             = 2.*(1-nu).*(1-2.*nu).*Fbar.*cotbeta.*cotbeta + (1-2.*nu).*y2./(Rbar+y3bar).*(-(1-2.*nu-a./Rbar).*cotbeta + y1./(Rbar+y3bar).*(nu+a./Rbar)) - ...
                    (1-2.*nu).*y2.*cotbeta./(Rbar+z3bar).*(1+a./(Rbar.*cosbeta)) - a.*y2.*(y3bar-a).*cotbeta./(Rbar.*Rbar.*Rbar) + ...
                    y2.*(y3bar-a)./(Rbar.*(Rbar+y3bar)).*((1-2.*nu).*cotbeta - 2.*nu.*y1./(Rbar+y3bar) - a.*y1./Rbar.*(1./Rbar+1./(Rbar+y3bar))) + ...
                    y2.*(y3bar-a).*cotbeta./(Rbar.*(Rbar+z3bar)).*(-2.*(1-nu).*cosbeta + (Rbar.*cosbeta+y3bar)./(Rbar+z3bar).*(1+a./(Rbar.*cosbeta)) + a.*y3bar./((Rbar.*Rbar).*cosbeta));
v3CB2             = -2.*(1-nu).*(1-2.*nu).*cotbeta .* (log(Rbar+y3bar)-cosbeta.*log(Rbar+z3bar)) - ...
                    2.*(1-nu).*y1./(Rbar+y3bar).*(2.*nu+a./Rbar) + 2.*(1-nu).*z1bar./(Rbar+z3bar).*(cosbeta+a./Rbar) + ...
                   (y3bar-a)./Rbar.*((1-2.*nu).*cotbeta-2.*nu.*y1./(Rbar+y3bar)-a.*y1./(Rbar.*Rbar)) - ...
                   (y3bar-a)./(Rbar+z3bar).*(cosbeta.*sinbeta + (Rbar.*cosbeta+y3bar).*cotbeta./Rbar.*(2.*(1-nu).*cosbeta - (Rbar.*cosbeta+y3bar)./(Rbar+z3bar)) + a./Rbar.*(sinbeta - y3bar.*z1bar./(Rbar.*Rbar) - z1bar.*(Rbar.*cosbeta+y3bar)./(Rbar.*(Rbar+z3bar))));
v1CB2             = v1CB2 ./ (4.*pi.*(1-nu));
v2CB2             = v2CB2 ./ (4.*pi.*(1-nu));
v3CB2             = v3CB2 ./ (4.*pi.*(1-nu));

v1B2              = v1InfB2 + v1CB2;
v2B2              = v2InfB2 + v2CB2;
v3B2              = v3InfB2 + v3CB2;


% Case III: Burgers vector (0,0,B3)
v1InfB3           = y2.*sinbeta.*((R.*sinbeta-y1)./(R.*(R-z3))+(Rbar.*sinbeta-y1)./(Rbar.*(Rbar+z3bar)));
v2InfB3           = (1-2.*nu).*sinbeta.*(log(R-z3)+log(Rbar+z3bar)) - (y2.*y2).*sinbeta.*(1./(R.*(R-z3))+1./(Rbar.*(Rbar+z3bar)));
v3InfB3           = 2.*(1-nu).*(F-Fbar) + y2.*sinbeta.*((R.*cosbeta-y3)./(R.*(R-z3))-(Rbar.*cosbeta+y3bar)./(Rbar.*(Rbar+z3bar)));
v1InfB3           = v1InfB3 ./ (8.*pi.*(1-nu));
v2InfB3           = v2InfB3 ./ (8.*pi.*(1-nu));
v3InfB3           = v3InfB3 ./ (8.*pi.*(1-nu));

v1CB3             = (1-2.*nu).*(y2./(Rbar+y3bar).*(1+a./Rbar) - y2.*cosbeta./(Rbar+z3bar).*(cosbeta+a./Rbar)) - ...
                    y2.*(y3bar-a)./Rbar.*(a./(Rbar.*Rbar) + 1./(Rbar+y3bar)) + ...
                    y2.*(y3bar-a).*cosbeta./(Rbar.*(Rbar+z3bar)).*((Rbar.*cosbeta+y3bar)./(Rbar+z3bar).*(cosbeta+a./Rbar) + a.*y3bar./(Rbar.*Rbar));
v2CB3             = (1-2.*nu).*(-sinbeta.*log(Rbar+z3bar) - y1./(Rbar+y3bar).*(1+a./Rbar) + z1bar./(Rbar+z3bar).*(cosbeta+a./Rbar)) + ...
                    y1.*(y3bar-a)./Rbar.*(a./(Rbar.*Rbar) + 1./(Rbar+y3bar)) - ...
                    (y3bar-a)./(Rbar+z3bar).*(sinbeta.*(cosbeta-a./Rbar) + z1bar./Rbar.*(1+a.*y3bar./(Rbar.*Rbar)) - ...
                    1./(Rbar.*(Rbar+z3bar)).*((y2.*y2).*cosbeta.*sinbeta - a.*z1bar./Rbar.*(Rbar.*cosbeta+y3bar)));
v3CB3             = 2.*(1-nu).*Fbar + 2.*(1-nu).*(y2.*sinbeta./(Rbar+z3bar).*(cosbeta + a./Rbar)) + ...
                    y2.*(y3bar-a).*sinbeta./(Rbar.*(Rbar+z3bar)).*(1 + (Rbar.*cosbeta+y3bar)./(Rbar+z3bar).*(cosbeta+a./Rbar) + a.*y3bar./(Rbar.*Rbar));
v1CB3             = v1CB3 ./ (4.*pi.*(1-nu));
v2CB3             = v2CB3 ./ (4.*pi.*(1-nu));
v3CB3             = v3CB3 ./ (4.*pi.*(1-nu));

v1B3              = v1InfB3 + v1CB3;
v2B3              = v2InfB3 + v2CB3;
v3B3              = v3InfB3 + v3CB3;


% Sum the for each slip component
v1                = B1.*v1B1 + B2.*v1B2 + B3.*v1B3;
v2                = B1.*v2B1 + B2.*v2B2 + B3.*v2B3;
v3                = B1.*v3B1 + B2.*v3B2 + B3.*v3B3;
end
%====================================================
function [v1 v2 v3] = adv_opt(y1, y2, y3, a, beta, nu, B1, B2, B3)
% These are the displacements in a uniform elastic half space due to slip
% on an angular dislocation (Comninou and Dunders, 1975).  Some of the
% equations for the B2 and B3 cases have been corrected following Thomas
% 1993.  The equations are coded in way such that they roughly correspond
% to each line in original text.  Exceptions have been made where it made 
% more sense because of grouping symbols.

sinbeta           = sin(beta);
cosbeta           = cos(beta);
cotbeta           = cot(beta);
z1                = y1.*cosbeta - y3.*sinbeta;
z3                = y1.*sinbeta + y3.*cosbeta;
R2                = y1.*y1 + y2.*y2 + y3.*y3;
R                 = sqrt(R2);
y3bar             = y3 + 2.*a;
z1bar             = y1.*cosbeta + y3bar.*sinbeta;
z3bar             = -y1.*sinbeta + y3bar.*cosbeta;
R2bar             = y1.*y1 + y2.*y2 + y3bar.*y3bar;
Rbar              = sqrt(R2bar);
F                 = -atan2(y2, y1) + atan2(y2, z1)    + atan2(y2.*R.*sinbeta, y1.*z1+(y2.*y2).*cosbeta);
Fbar              = -atan2(y2, y1) + atan2(y2, z1bar) + atan2(y2.*Rbar.*sinbeta, y1.*z1bar+(y2.*y2).*cosbeta);
%
nu2=1-2.*nu;
R3b=1./(R2bar.*Rbar);
acR3=a.*cotbeta.*R3b;
y3bar_a=y3bar-a;
ctb2=cotbeta.*cotbeta;
Ry3b=Rbar+y3bar;
Rz3b=Rbar+z3bar;
lRy3b=log(Ry3b);
lRz3b=log(Rz3b);
lRz3=log(R-z3)+lRz3b;
y23b=y2.*y3bar_a;
aRb=a./Rbar;
RRz=1./(Rbar.*Rz3b);
RRy=1./(Rbar.*Ry3b);
Rcb=Rbar.*cosbeta;
Rcy=Rcb+y3bar;
caRb=cosbeta+aRb;
llcl=nu2.*(log(R-y3)+lRy3b-cosbeta.*lRz3);
rR2b=1./R2bar;
Pn8=1./(8.*pi.*(1-nu));
%
% Case I: Burgers vector (B1,0,0)
v1InfB1           = 2.*(1-nu).*(F+Fbar)-y1.*y2.*(1./(R.*(R-y3))+RRy)-y2.*cosbeta.*((R.*sinbeta-y1)./(R.*(R-z3))+(Rbar.*sinbeta-y1).*RRz);
v2InfB1           = llcl-y2.*y2.*(1./(R.*(R-y3))+RRy-cosbeta.*(1./(R.*(R-z3))+RRz));
v3InfB1           = y2.*(1./R-1./Rbar-cosbeta.*((R.*cosbeta-y3)./(R.*(R-z3))-Rcy.*RRz));

v1CB1             = -2.*(1-nu).*nu2.*Fbar.*ctb2+...
                    nu2.*y2./Ry3b.*((1-2.*nu-aRb).*cotbeta-...
                    y1./Ry3b.*(nu+aRb))+...
                    nu2.*y2.*cosbeta.*cotbeta./Rz3b.*caRb+...
                    y23b.*acR3+...
                    y23b.*RRy.*(-nu2.*cotbeta+y1./Ry3b.*(2.*nu+aRb)+a.*y1./R2bar)+...
                    y23b.*RRz.*(cosbeta./Rz3b.*(Rcy.*(nu2.*cosbeta-aRb).*cotbeta+2.*(1-nu).*(Rbar.*sinbeta-y1).*cosbeta)-a.*y3bar.*cosbeta.*cotbeta./R2bar);
v2CB1             = nu2.*((2.*(1-nu).*ctb2-nu).*lRy3b-(2.*(1-nu).*ctb2+1-2.*nu).*cosbeta.*lRz3b)-...
                    nu2./Ry3b.*(y1.*cotbeta.*(1-2.*nu-aRb)+...
                    nu.*y3bar-a+(y2.*y2)./Ry3b.*(nu+aRb))-...
                    nu2.*z1bar.*cotbeta./Rz3b.*caRb-...
                    y3bar_a.*y1.*acR3+...
                    y3bar_a./Ry3b.*(-2.*nu+1./Rbar.*(nu2.*y1.*cotbeta-a)+(y2.*y2).*RRy.*(2.*nu+aRb)+a.*(y2.*y2).*R3b)+...
                    y3bar_a./Rz3b.*((cosbeta.*cosbeta)-1./Rbar.*(nu2.*z1bar.*cotbeta+a.*cosbeta)+y3bar.*z1bar.*acR3-RRz.*((y2.*y2).*(cosbeta.*cosbeta)-a.*z1bar.*cotbeta./Rbar.*Rcy));
v3CB1             = 2.*(1-nu).*((nu2.*Fbar.*cotbeta)+(y2./Ry3b.*(2.*nu+aRb))-(y2.*cosbeta./Rz3b.*caRb))+...
                    y23b./Rbar.*(2.*nu./Ry3b+a./R2bar)+...
                    y23b.*cosbeta.*RRz.*(1-2.*nu-Rcy./Rz3b.*caRb - a.*y3bar./R2bar);

v1B1              = Pn8.*(v1InfB1 + 2.*v1CB1);
v2B1              = Pn8.*(v2InfB1 + 2.*v2CB1);
v3B1              = Pn8.*(v3InfB1 + 2.*v3CB1);


% Case II: Burgers vector (0,B2,0)
v1InfB2           = -llcl+y1.*y1.*(1./(R.*(R-y3))+RRy)+z1.*(R.*sinbeta-y1)./(R.*(R-z3))+z1bar.*(Rbar.*sinbeta-y1).*RRz;
v2InfB2           = 2.*(1-nu).*(F+Fbar) + y1.*y2.*(1./(R.*(R-y3))+RRy) - y2.*(z1./(R.*(R-z3))+z1bar.*RRz);
v3InfB2           = -nu2.*sinbeta.*lRz3-y1.*(1./R-1./Rbar) + z1.*(R.*cosbeta-y3)./(R.*(R-z3)) - z1bar.*Rcy.*RRz;

v1CB2             = nu2.*((2.*(1-nu).*ctb2+nu).*lRy3b-(2.*(1-nu).*ctb2+1).*cosbeta.*lRz3b) + ...
                    nu2./Ry3b.*(-nu2.*y1.*cotbeta+nu.*y3bar-a+a.*y1.*cotbeta./Rbar+(y1.*y1)./Ry3b.*(nu+aRb))-...
                    nu2.*cotbeta./Rz3b.*(z1bar.*cosbeta-a.*(Rbar.*sinbeta-y1)./Rcb)-y1.*y3bar_a.*acR3 + ...
                    y3bar_a./Ry3b.*(2.*nu + 1./Rbar.*(nu2.*y1.*cotbeta+a) - (y1.*y1).*RRy.*(2.*nu+aRb) - a.*(y1.*y1).*R3b) + ...
                    y3bar_a.*cotbeta./Rz3b.*(-cosbeta.*sinbeta+a.*y1.*y3bar./(R2bar.*Rcb) + (Rbar.*sinbeta-y1)./Rbar.*(2.*(1-nu).*cosbeta - Rcy./Rz3b.*(1+a./Rcb)));
v2CB2             = 2.*(1-nu).*nu2.*Fbar.*ctb2 + nu2.*y2./Ry3b.*(-(1-2.*nu-aRb).*cotbeta + y1./(Rbar+y3bar).*(nu+aRb)) - ...
                    nu2.*y2.*cotbeta./Rz3b.*(1+a./Rcb) - y2.*y3bar_a.*acR3 + ...
                    y23b.*RRy.*(nu2.*cotbeta - 2.*nu.*y1./Ry3b - a.*y1./Rbar.*(1./Rbar+1./Ry3b)) + ...
                    y23b.*cotbeta.*RRz.*(-2.*(1-nu).*cosbeta + Rcy./Rz3b.*(1+a./Rcb) + a.*y3bar./(R2bar.*cosbeta));
v3CB2             = -2.*(1-nu).*nu2.*cotbeta.*(lRy3b-cosbeta.*lRz3b)-2.*(1-nu).*y1./Ry3b.*(2.*nu+aRb)+2.*(1-nu).*z1bar./Rz3b.*caRb + y3bar_a./Rbar.*(nu2.*cotbeta-2.*nu.*y1./Ry3b-a.*y1./R2bar) - ...
                    y3bar_a./Rz3b.*(cosbeta.*sinbeta+Rcy.*cotbeta./Rbar.*(2.*(1-nu).*cosbeta - Rcy./Rz3b)+aRb.*(sinbeta-y3bar.*z1bar./R2bar-z1bar.*Rcy./(RRz)));

v1B2              = Pn8.*(v1InfB2 + 2.*v1CB2);
v2B2              = Pn8.*(v2InfB2 + 2.*v2CB2);
v3B2              = Pn8.*(v3InfB2 + 2.*v3CB2);


% Case III: Burgers vector (0,0,B3)
v1InfB3           = y2.*sinbeta.*((R.*sinbeta-y1)./(R.*(R-z3))+(Rbar.*sinbeta-y1).*RRz);
v2InfB3           = nu2.*sinbeta.*lRz3-(y2.*y2).*sinbeta.*(1./(R.*(R-z3))+RRz);
v3InfB3           = 2.*(1-nu).*(F-Fbar) + y2.*sinbeta.*((R.*cosbeta-y3)./(R.*(R-z3))-Rcy.*RRz);

v1CB3             = nu2.*(y2./Ry3b.*(1+aRb)-y2.*cosbeta./Rz3b.*caRb)-y23b./Rbar.*(a./R2bar + 1./Ry3b) + ...
                    y23b.*cosbeta.*RRz.*(Rcy./Rz3b.*caRb+a.*y3bar./R2bar);
v2CB3             = nu2.*(-sinbeta.*lRz3b-y1./Ry3b.*(1+aRb)+z1bar./Rz3b.*caRb)+...
                    y1.*y3bar_a./Rbar.*(a./R2bar+1./Ry3b)-y3bar_a./Rz3b.*(sinbeta.*(cosbeta-aRb)+z1bar./Rbar.*(1+a.*y3bar./R2bar) - ...
                    RRz.*((y2.*y2).*cosbeta.*sinbeta - a.*z1bar./Rbar.*Rcy));
v3CB3             = 2.*(1-nu).*Fbar+2.*(1-nu).*(y2.*sinbeta./Rz3b.*(cosbeta+aRb))+...
                    y2.*y3bar_a.*sinbeta.*RRz.*(1 + Rcy./Rz3b.*caRb + a.*y3bar./R2bar);

v1B3              = Pn8.*(v1InfB3 + 2.*v1CB3);
v2B3              = Pn8.*(v2InfB3 + 2.*v2CB3);
v3B3              = Pn8.*(v3InfB3 + 2.*v3CB3);


% Sum the for each slip component
v1                = B1.*v1B1 + B2.*v1B2 + B3.*v1B3;
v2                = B1.*v2B1 + B2.*v2B2 + B3.*v2B3;
v3                = B1.*v3B1 + B2.*v3B2 + B3.*v3B3;
end
%====================================================