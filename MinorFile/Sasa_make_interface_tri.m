function mk_green_tri(GPS,sub_f,bound_f,n_mesh)
[s]=init_interface_tri(sub_f,bound_f,n_mesh.*10);
[s]=down_tri(s,GPS,n_mesh);
end
%====================================================
function [s]=init_interface_tri(sub_f,bound_f,int_mesh)
%====================================================
Fid=fopen(/home_tmp/sasajima/DATA/PACdw.txt);
dep_sub=textscan(Fid,'%f%f%f');
fclose(Fid);
dep_sub=cell2mat(dep_sub);
%====================================================
Fid=fopen(/home_tmp/sasajima/DATA/PACbound.txt);
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
int_mesh=10000;
while n<int_mesh
  slat=(max_lat-min_lat).*rand(1)+min_lat;
  slon=(max_lon-min_lon).*rand(1)+min_lon;
  ID=inpolygon(slon,slat,bound(:,1),bound(:,2));
  if ID==1
    n=n+1;
    s.lat(n)=slat;
    s.lon(n)=slon;
    s.dep(n)=F(slon,slat);
    %if rem(n,round(int_mesh/10))==1;
      %plot3(s.lon,s.lat,s.dep,'.')
      %pause(.1)
    %end
  end
end
fclose(Fid);
%plot3(s.lon,s.lat,s.dep,'.')
%====================================================
tri = delaunay(s.lon,s.lat);
%====================================================
ntri=length(tri);
nn=0;
Fid=fopen('/home_tmp/sasajima/DATA/PAC_TOHOKU/PACtri.dat','w')
for n=1:ntri
  glon=mean(s.lon(tri(n,:)));  
  glat=mean(s.lat(tri(n,:)));
  gdep=mean(s.dep(tri(n,:)));
  ID=inpolygon(glon,glat,bound(:,1),bound(:,2));

  fprintf(Fid,'%5d %11.5f %11.5f %11.5f\n',n,glon,glat,gdep);

  if ID==1
    nn=nn+1;
    s.tri(nn,:)=tri(n,:);
  end  
end
figure
plot(bound(:,1),bound(:,2),'r')
hold on
triplot(s.tri,s.lon,s.lat);
end
fclose(Fid);
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
%==========================================
%Calculate the tri poin and center point




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
      [ux1 uy1 uz1]             = adv_opt(sx1, sy1, sz-z(iTri), z(iTri), beta, pr, lss, lts, lds);
                                   
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
%Observ point
function [S]=down_tri(S,GPS,n_mesh)
alat0=(mean(GPS.lat)+mean(S.lat))./2;
alon0=(mean(GPS.lon)+mean(S.lon))./2;
[gx,gy]=PLTXY(GPS.lat,GPS.lon,alat0,alon0);
[sx,sy]=PLTXY(S.lat,S.lon,alat0,alon0);
gz=GPS.hight./1000;
sz=S.dep;
Ua=zeros(length(S.tri),1);
for n=1:length(S.tri)
  [U]=CalcTriDisps(gx,gy,gz,sx(S.tri(n)),sy(S.tri(n)),sz(S.tri(n)),0.25,0,1,0);
  Ua(n)=sqrt(U.x.^2+U.y.^2+U.z.^2);
end
Ntri=length(S.tri);
while Ntri > n_mesh
  r_index=ones(length(S.lat,1));
  [~,index]=min(Ua);
  min_tri=S.tri(index);
  f_tri=zeros(length(f_tri),1);
  for n=1:3
    f_tri(n)=sum(find(S.tri,min_tri(n)));
  end
  [~,index]=max(f_tri);
  r_index(min_tri(index))=0;
  S.lat=S.lat(r_index);
  S.lon=S.lon(r_index);
  S.dep=S.dep(r_index);
  S.tri=delaunay(S.lon,S.lat);
  [sx,sy]=PLTXY(S.lat,S.lon,alat0,alon0);
  sz=S.dep;
  Ua=zeros(length(S.tri),1);
  for n=1:length(S.tri)
  tic
    [U]=CalcTriDisps(gx,gy,gz,sx(S.tri(n)),sy(S.tri(n)),sz(S.tri(n)),0.25,0,1,0);
  toc
    Ua(n)=sqrt(U.x.^2+U.y.^2+U.z.^2);
  end
  n_mesh=length(S.tri);
end
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
%====================================================%====================================================
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
