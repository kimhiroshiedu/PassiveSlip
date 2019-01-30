function [Vneu]=euler2vel(olat,olon,plat,plon,pomega)
%
% function Vneu = euler2vel(olat,olon,plat,plon,pomega)
%
% This function computes the linear velocity Vneu(vx,vy,vz) 
% due to angular rotation (pomega[DEG/MYR]) 
% at a set of input points P(plat,plon).
% See Cox and Hart (1986), p. 155.
%
% INPUT:
%    plat, plon, pomega = euler pole[deg], omega
%    plat, olon         = calc points
% OUTPUT:
%    Vneu = local velocities (mm/yr)
%
% euler pole : convert (lat lon omega) => vector
%               deg/Myr --> rad/yr
%
% JGR 21,20,2191-2194, 1994 NNR-NUVEL-1A 
% Plate  Lat   Lon    Omega(deg/Myr) errors 
% na-pa  48.7   -78.2  0.75   1.3  1.2  -61  0.01
% ri-pa  31.0  -102.4  2.45   3.6  0.6   21  0.57 
% co-pa  36.8  -108.6  2.00   1.0  0.6  -33  0.05 
% ri-na  22.8  -109.4  1.80   1.8  0.6  -57  0.58 
% ri-co   6.8   -83.7  0.54  35.8  1.8  -56  0.52 
% co-na  27.9  -120.7  1.36   1.8  0.7  -67  0.05 
% co-nz   4.8  -124.3  0.91   2.9  1.5  -88  0.05 
% nz-pa  55.6   -90.1  1.36   1.8  0.9   -1  0.02 
% nz-an  40.5   -95.9  0.52   4.5  1.9   -9  0.02 
% nz-sa  56.0   -94.0  0.72   3.6  1.5  -10  0.02 
% an-pa  64.3   -84.0  0.87   1.2  1.0   81  0.01 
% pa-au -60.1  -178.3  1.07   1.0  0.9  -58  0.01 
% eu-pa  61.1   -85.8  0.86   1.3  1.1   90  0.02 
% co-ca  24.1  -119.4  1.31   2.5  1.2  -60  0.05 
% nz-ca  56.2  -104.6  0.55   6.5  2.2  -31  0.03 
% Atlantic  Ocean 
% eu-na  62.4   135.8  0.21   4.1  1.3  -11  0.01 
% af-na  78.8    38.3  0.24   3.7  1.0   77  0.01 
% af-eu  21.0   -20.6  0.12   6.0  0.7   -4  0.02 
% na-sa  16.3   -58.1  0.15   5.9  3.7   -9  0.01 
% af-sa  62.5   -39.4  0.31   2.6  0.8  -11  0.01 
% an-sa  86.4   -40.7  0.26   3.0  1.2  -24  0.01 
% na-ca -74.3   -26.1  0.10  24.7  2.6  -52  0.03 
% ca-sa  50.0   -65.3  0.18  14.9  4.3   -2  0.03 
% Indian  Ocean 
% au-an  13.2    38.2  0.65   1.3  1.0  -63  0.01 
% af-an   5.6   -39.2  0.13   4.4  1.3  -42  0.01 
% au-af  12.4    49.8  0.63   1.2  0.9  -39  0.01 
% au-in  -5.6    77.1  0.30   7.4  3.1  -43  0.07 
% in-af  23.6    28.5  0.41   8.8  1.5  -74  0.06 
% ar-af  24.1    24.0  0.40   4.9  1.3  -65  0.05 
% in-eu  24.4    17.7  0.51   8.8  1.8  -79  0.05 
% ar-eu  24.6    13.7  0.50   5.2  1.7  -72  0.05 
% au-eu  15.1    40.5  0.69   2.1  1.1  -45  0.01 
% in-ar   3.0    91.5  0.03  25.2  2.4  -58  0.04
%
[NN,MM]=size(olat);
olat=olat(:);
olon=olon(:);
%
[px,py,pz]=ell2xyz(plat,plon,0);
pvec=[px, py, pz];
pvec=pomega.*pvec./norm(pvec).*(1e-6.*pi./180);
%
% local observation points : convert (lat lon) = xyz
%                             m --> mm
[ox,oy,oz]=ell2xyz(olat,olon,0);
Oxyz=[ox oy oz];
Oxyz = Oxyz * 1e3;
%
% V = (w E) x (r P), where E is the Euler pole unit vector
Vxyz(:,1) = Oxyz(:,2)*pvec(3) - pvec(2)*Oxyz(:,3);
Vxyz(:,2) = Oxyz(:,3)*pvec(1) - pvec(3)*Oxyz(:,1);
Vxyz(:,3) = Oxyz(:,1)*pvec(2) - pvec(1)*Oxyz(:,2);
%
Vneu=zeros(3,NN*MM);
for N=1:NN*MM
  Vneu(:,N)=xyz2neu([olat(N),olon(N)],Vxyz(N,:));
end
Vneu=reshape(Vneu,3,NN,MM);
%
end
%==============================================================
function [x,y,z]=ell2xyz(lat,lon,h)
% ELL2XYZ  Converts ellipsoidal coordinates to cartesian. Vectorized.
% GRS80
a=6378137.0;
f=1./298.257222101;
e2=1-(1-f)^2;
%
rad=pi/180;
lat=lat.*rad;
lon=lon.*rad;
%
v=a./sqrt(1-e2*sin(lat).*sin(lat));
x=(v+h).*cos(lat).*cos(lon);
y=(v+h).*cos(lat).*sin(lon);
z=(v.*(1-e2)+h).*sin(lat);
end
%==============================================================
function dneu=xyz2neu(bl,dxyz)
%
% TRANSFORMATION FROM (DX,DY,DZ) => (DN,DE,DU)
% CODE BY T.ITO 2006/12/13 ver0.1
% BUG FIX BY T.ITO 2007/01/13 ver0.2
%
deg2rad=pi/180;
bl=bl(1:2).*deg2rad;
dneu=zeros(size(dxyz));
R=[-sin(bl(1)).*cos(bl(2)) -sin(bl(1)).*sin(bl(2)) cos(bl(1)) ; ...
   -sin(bl(2))              cos(bl(2))             0.0        ; ...
    cos(bl(1)).*cos(bl(2))  cos(bl(1)).*sin(bl(2)) sin(bl(1))];
for j=1:length(dxyz(:,1))
    dneu(j,:)=R*dxyz(j,:)';
end
end
