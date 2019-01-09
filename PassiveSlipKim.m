%% PassiveSlip.m
function PassiveSlipKim
% Coded by Ryohei Sasajima final 2013/12/23
% Combined by Hiroshi Kimura 2018/11/12

[triC,tri3,tri,sll]=make_test_trill;
save('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/tri','triC','tri3','tri','sll');
load('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/tri','triC','tri3','tri','sll');

ALON0=142.5;
ALAT0=38.3;
[trixyzC,trixyz3,sxyz]=trill2trixyz(triC,tri3,sll,ALAT0,ALON0);
save('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/trixyz','trixyzC','trixyz3','sxyz');
load('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/trixyz','trixyzC','trixyz3','sxyz');

[sitaS,sitaD,normVec]=strike_dip(trixyzC,trixyz3);
save('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/sita','sitaS','sitaD','normVec');
load('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/sita','sitaS','sitaD','normVec');

% figure(31); quiver(trixyzC(:,1),-trixyzC(:,2),sitaS,sitaD,1,'b')

[xyz]=makexyz;
save('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/xyz','xyz');
load('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/xyz','xyz');

% [si]=makeG_s_Q(trixyzC,trixyz3,sitaS,sitaD,normVec);
% [di]=makeG_d_Q(trixyzC,trixyz3,sitaS,sitaD,normVec);

% [si]=makeG_s_O(xyz,trixyz3);
% [di]=makeG_d_O(xyz,trixyz3);

[gu]=makeGreenDisp(xyz,trixyz3);
[gs]=makeGreenStrain(trixyzC,trixyz3,sitaS,sitaD,normVec);

% [sUxyz,dUxyz]=loadMAT2(xyz);
% [sSsn,sSdn,dSsn,dSdn]=loadMAT(trixyzC);

[Slip]=defineSlipQ(triC,sitaS);
save('/home_tmp/sasajima/DATA/MAT/sdSlip_Q1','Slip');
load('/home_tmp/sasajima/DATA/MAT/sdSlip_Q1','Slip');

[xySlip]=Slip2xyll(Slip,sitaS);
save('/home_tmp/sasajima/DATA/MAT/xySlip_Q1','xySlip');
load('/home_tmp/sasajima/DATA/MAT/xySlip_Q1','xySlip');
figure(103); quiver(triC(:,1),triC(:,2),xySlip(:,1),-xySlip(:,2),'b')
figure(101); quiver(trixyzC(:,1),-trixyzC(:,2),xySlip(:,1),-xySlip(:,2),'g')

ll(:,1:2)=triC(:,1:2);
[vnorm,v,NARFA]=SasaEular2Velo(ll);

%{
[]=relavec2vecVxy(NARFA)
[xyBslip]=BackSlipxy(mcmcC,triC,trixyzC,sitaS,accuV,vecVxy)
%}

[A_Ssn,A_Sdn]=I_CalcStrain(Slip,sSsn,sSdn,dSsn,dSdn);
save('/home_tmp/sasajima/DATA/MAT/A_Ssn_2','A_Ssn','A_Sdn');
load('/home_tmp/sasajima/DATA/MAT/A_Sdn_2','A_Ssn','A_Sdn');
figure(201); quiver(trixyzC(:,1),-trixyzC(:,2),A_Ssn,A_Sdn,5,'b')

[F_Slip]=OutofAsperity01(A_Ssn,A_Sdn,Slip,sSsn,sSdn,dSsn,dSdn);
save('/home_tmp/sasajima/DATA/MAT/FSlip','F_Slip');
load('/home_tmp/sasajima/DATA/MAT/FSlip','F_Slip');

[xyFSlip]=FSlip2xyll(F_Slip,sitaS);
save('/home_tmp/sasajima/DATA/MAT/xyFSlip_Q1','xyFSlip');
load('/home_tmp/sasajima/DATA/MAT/xyFSlip_Q1','xyFSlip');

Fid=fopen('./QxySlip11.dat','w');
n=size(Slip,1);

for i=1:n
    fprintf(Fid, '%11.3f %11.3f %9.5f %9.5f\n',trixyzC(i,1), trixyzC(i,2), xyFSlip(i,1), xyFSlip(i,2));
end
fclose(Fid);

figure(1); quiver(trixyzC(:,1),-trixyzC(:,2),xyFSlip(:,1),-xyFSlip(:,2),2,'g')
hold on
triplot(tri,1000.*sxyz(:,1),1000.*sxyz(:,2));

[F_Ssn,F_Sdn]=I_CalcStrainF(F_Slip,sSsn,sSdn,dSsn,dSdn);
save('/home_tmp/sasajima/DATA/MAT/FStrain','F_Ssn','F_Sdn');
load('/home_tmp/sasajima/DATA/MAT/FStrain','F_Ssn','F_Sdn');
figure(401); quiver(trixyzC(:,1),-trixyzC(:,2),F_Ssn,F_Sdn,5,'b')

absO_Ssn=abs(O_Ssn);
absO_Sdn=abs(O_Sdn);
sumAOSsn=sum(absO_Ssn,1);
sumAOSdn=sum(absO_Sdn,1);

[FO_Ssn,FO_Sdn]=OutofAsperity11(F_Ssn,F_Sdn,Slip);
save('/home_tmp/sasajima/DATA/MAT/FOStrain','FO_Ssn','FO_Sdn');
load('/home_tmp/sasajima/DATA/MAT/FOStrain','FO_Ssn','FO_Sdn');
absFO_Ssn=abs(FO_Ssn);
absFO_Sdn=abs(FO_Sdn);
sumF0Ssn=sum(absFO_Ssn,1);
sumF0Sdn=sum(absFO_Sdn,1);

end

%% CalcTriDisps.m
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
% Implements algorithms described in the journal article:
% Meade, B. J. Algorithms for calculating displacements,
% strains, and stresses for triangular dislocation elements
% in a uniform elastic half space
% Computers and Geosciences, submitted, 2006.
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


% Calculate the slip vector in XYZ coordinates
normVec                      = cross([x(2);y(2);z(2)]-[x(1);y(1);z(1)], [x(3);y(3);z(3)]-[x(1);y(1);z(1)]);
normVec                      = normVec./norm(normVec);

[n1,n2]                      = size(sz);
sz1                          = sz(:);

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

for iTri = 1:3
    % Calculate strike and dip of current leg
    strike                   = 180/pi*(atan2(y(iTri+1)-y(iTri), x(iTri+1)-x(iTri)));
    segMapLength             = sqrt((x(iTri)-x(iTri+1))^2 + (y(iTri)-y(iTri+1))^2);
    [rx ry]                  = RotateXyVec(x(iTri+1)-x(iTri), y(iTri+1)-y(iTri), -strike);
    dip                      = 180/pi*(atan2(z(iTri+1)-z(iTri), rx));
    
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
    
    ssVec                    = [cos(strike/180*pi) sin(strike/180*pi) 0];
    tsVec                    = [-sin(strike/180*pi) cos(strike/180*pi) 0];
    dsVec                    = cross(ssVec, tsVec);
    lss                      = dot(slipVec, ssVec);
    lts                      = dot(slipVec, tsVec);
    lds                      = dot(slipVec, dsVec);
    
    
    if (abs(beta) > 0.000001) && (abs(beta-pi) > 0.000001)
        % First angular dislocation
        [sx1 sy1]                 = RotateXyVec(sx-x(iTri), sy-y(iTri), -strike);
        sx1(abs(sx1)<0.0001)=0.0001;%for eliminate NAN by R. Sasajima and T. Ito ,Nagoya. U. in 2012%
        sy1(abs(sy1)<0.0001)=0.0001;
        [ux1 uy1 uz1]             = adv(sx1, sy1, sz1-z(iTri), z(iTri), beta, pr, lss, lts, lds);
        
        % Second angular dislocation
        [sx2 sy2]                 = RotateXyVec(sx-x(iTri+1), sy-y(iTri+1), -strike);
        sx2(abs(sx2)<0.0001)=0.0001;%for eliminate NAN by R. Sasajima and T. Ito ,Nagoya. U. in 2012%
        sy2(abs(sy2)<0.0001)=0.0001;
        [ux2 uy2 uz2]             = adv(sx2, sy2, sz1-z(iTri+1), z(iTri+1), beta, pr, lss, lts, lds);
        
        % Rotate vectors to correct for strike
        [uxn uyn]                 = RotateXyVec(ux1-ux2, uy1-uy2, strike);
        uzn                       = uz1-uz2;
        
        % Add the displacements from current leg
        U.x                       = U.x + reshape(uxn,n1,n2);
        U.y                       = U.y + reshape(uyn,n1,n2);
        U.z                       = U.z + reshape(uzn,n1,n2);
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
end

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

function [xp yp] = RotateXyVec(x, y, alpha)
% Rotate a vector by an angle alpha
x                             = x(:);
y                             = y(:);
alpha                         = pi/180*alpha;
xp                            = cos(alpha).*x - sin(alpha).*y;
yp                            = sin(alpha).*x + cos(alpha).*y;
end

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
%% CalcTriStrains.m
function [S] = CalcTriStrains(sx, sy, sz, x, y, z, pr, ss, ts, ds)
% CalcTriStrains.m
%
% Calculates strains due to slip on a triangular dislocation in an
% elastic half space utilizing the symbolically differentiated
% displacement gradient tensor derived from the expressions for
% the displacements due to an angular dislocation in an elastic half
% space (Comninou and Dunders, 1975).
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
%  S  : structure containing the strains (S.xx, S.yy, S.zz, S.xy, S.xz, S.yz)
%
% This paper should and related code should be cited as:
% Brendan J. Meade, Algorithms for the calculation of exact
% displacements, strains, and stresses for Triangular Dislocation
% Elements in a uniform elastic half space, Computers &
% Geosciences (2007), doi:10.1016/j.cageo.2006.12.003.
%
% Use at your own risk and please let me know of any bugs/errors.
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
S.xx                         = zeros(size(sx));
S.yy                         = zeros(size(sx));
S.zz                         = zeros(size(sx));
S.xy                         = zeros(size(sx));
S.xz                         = zeros(size(sx));
S.yz                         = zeros(size(sx));

% Add a copy of the first vertex to the vertex list for indexing
x(4)                         = x(1);
y(4)                         = y(1);
z(4)                         = z(1);

for iTri = 1:3
    % Calculate strike and dip of current leg
    strike                   = 180/pi*(atan2(y(iTri+1)-y(iTri), x(iTri+1)-x(iTri)));
    segMapLength             = sqrt((x(iTri)-x(iTri+1))^2 + (y(iTri)-y(iTri+1))^2);
    [rx ry]                  = RotateXyVec(x(iTri+1)-x(iTri), y(iTri+1)-y(iTri), -strike);
    dip                      = 180/pi*(atan2(z(iTri+1)-z(iTri), rx));
    
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
    ssVec                    = [cos(strike/180*pi) sin(strike/180*pi) 0];
    tsVec                    = [-sin(strike/180*pi) cos(strike/180*pi) 0];
    dsVec                    = cross(ssVec, tsVec);
    lss                      = dot(slipVec, ssVec);
    lts                      = dot(slipVec, tsVec);
    lds                      = dot(slipVec, dsVec);
    
    if (abs(beta) > 0.000001) && (abs(beta-pi) > 0.000001)
        % First angular dislocation
        [sx1 sy1]                 = RotateXyVec(sx-x(iTri), sy-y(iTri), -strike);
        sx1(abs(sx1)<0.0001)=0.0001;%for eliminate NAN by R. Sasajima and T. Ito ,Nagoya. U. in 2012%
        sy1(abs(sy1)<0.0001)=0.0001;
        [a11 a22 a33 a12 a13 a23] = advs(sx1, sy1, sz-z(iTri), z(iTri), beta, pr, lss, lts, lds);
        
        % Second angular dislocation
        [sx2 sy2]                 = RotateXyVec(sx-x(iTri+1), sy-y(iTri+1), -strike);
        sx2(abs(sx2)<0.0001)=0.0001;
        sy2(abs(sy2)<0.0001)=0.0001;
        [b11 b22 b33 b12 b13 b23] = advs(sx2, sy2, sz-z(iTri+1), z(iTri+1), beta, pr, lss, lts, lds);
        
        % Rotate tensors to correct for strike
        bxx                       = a11-b11;
        byy                       = a22-b22;
        bzz                       = a33-b33;
        bxy                       = a12-b12;
        bxz                       = a13-b13;
        byz                       = a23-b23;
        
        g                         = pi/180*strike;
        e11n                      = (cos(g)*bxx-sin(g)*bxy)*cos(g)-(cos(g)*bxy-sin(g)*byy)*sin(g);
        e12n                      = (cos(g)*bxx-sin(g)*bxy)*sin(g)+(cos(g)*bxy-sin(g)*byy)*cos(g);
        e13n                      = cos(g)*bxz-sin(g)*byz;
        e22n                      = (sin(g)*bxx+cos(g)*bxy)*sin(g)+(sin(g)*bxy+cos(g)*byy)*cos(g);
        e23n                      = sin(g)*bxz+cos(g)*byz;
        e33n                      = bzz;
        
        % Add the strains from current leg
        S.xx                      = S.xx + e11n;
        S.yy                      = S.yy + e22n;
        S.zz                      = S.zz + e33n;
        S.xy                      = S.xy + e12n;
        S.xz                      = S.xz + e13n;
        S.yz                      = S.yz + e23n;
    end
end
end

function [a b] = swap(a, b)
% Swap two values
temp                            = a;
a                               = b;
b                               = temp;
end

function [e11 e22 e33 e12 e13 e23] = advs(y1, y2, y3, a, b, nu, B1, B2, B3)
% These are the strains in a uniform elastic half space due to slip
% on an angular dislocation.  They were calculated by symbolically
% differentiating the expressions for the displacements (Comninou and
% Dunders, 1975, with typos noted by Thomas 1993) then combining the
% elements of the displacement gradient tensor to form the strain tensor.

e11 = B1.*(1./8.*((2-2.*nu).*(2.*y2./y1.^2./(1+y2.^2./y1.^2)-y2./(y1.*cos(b)-y3.*sin(b)).^2.*cos(b)./(1+y2.^2./(y1.*cos(b)-y3.*sin(b)).^2)+(y2./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).*y1-y2.*(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2.*(2.*y1.*cos(b)-y3.*sin(b)))./(1+y2.^2.*(y1.^2+y2.^2+y3.^2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2)-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*cos(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*y1-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*(2.*y1.*cos(b)+(y3+2.*a).*sin(b)))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))-y2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))-y1.*y2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y1-1./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*y1-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*y1-1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*y1)-y2.*cos(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b).*y1-1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y1-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))+(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).*y1-1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y1-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))))./pi./(1-nu)+1./4.*((-2+2.*nu).*(1-2.*nu).*(y2./y1.^2./(1+y2.^2./y1.^2)-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*cos(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*y1-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*(2.*y1.*cos(b)+(y3+2.*a).*sin(b)))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2)).*cot(b).^2-(1-2.*nu).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((1-2.*nu-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)-y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1+(1-2.*nu).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1.*cot(b)-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))-(1-2.*nu).*y2.*cos(b).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-(1-2.*nu).*y2.*cos(b).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1-3.*a.*y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*y1-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((-1+2.*nu).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2)).*y1-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((-1+2.*nu).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2)).*y1+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)+a./(y1.^2+y2.^2+(y3+2.*a).^2)-2.*a.*y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^2)-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+(2-2.*nu).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1).*cos(b))-a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2)).*y1-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+(2-2.*nu).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1).*cos(b))-a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+(2-2.*nu).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y1.*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1.*cot(b)+(2-2.*nu).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).*y1-1).*cos(b))+2.*a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y1))./pi./(1-nu))+B2.*(1./8.*((-1+2.*nu).*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-cos(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))))+2.*y1.*(1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))+y1.^2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y1-1./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*y1-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*y1-1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*y1)+cos(b).*((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+(y1.*cos(b)-y3.*sin(b)).*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b).*y1-1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y1-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))+cos(b).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+(y1.*cos(b)+(y3+2.*a).*sin(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).*y1-1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y1-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b)))./pi./(1-nu)+1./4.*((1-2.*nu).*(((2-2.*nu).*cot(b).^2+nu)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-((2-2.*nu).*cot(b).^2+1).*cos(b).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((-1+2.*nu).*y1.*cot(b)+nu.*(y3+2.*a)-a+a.*y1.*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1+(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((-1+2.*nu).*cot(b)+a.*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-a.*y1.^2.*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)+2.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-y1.^3./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.^3./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))+(1-2.*nu).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*((y1.*cos(b)+(y3+2.*a).*sin(b)).*cos(b)-a.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-(1-2.*nu).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b).^2-a.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).*y1-1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)+a.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*y1)-a.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)+3.*a.*y1.^2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2)-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*y1.*cot(b)+a)-y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*y1.*cot(b)+a).*y1+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1-2.*nu).*cot(b)-2.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+y1.^3./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+y1.^3./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+y1.^3./(y1.^2+y2.^2+(y3+2.*a).^2).^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a-2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1+3.*a.*y1.^3./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2))-(y3+a).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(-cos(b).*sin(b)+a.*y1.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)))).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+(y3+a).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b)-3.*a.*y1.^2.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2)./cos(b)+(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).*y1-1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))).*y1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*y1)))./pi./(1-nu))+B3.*(1./8.*y2.*sin(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b).*y1-1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y1-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))+(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).*y1-1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y1-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b)))./pi./(1-nu)+1./4.*((1-2.*nu).*(-y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1+y2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+y2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1)+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(a./(y1.^2+y2.^2+(y3+2.*a).^2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)).*y1-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y1-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1)-y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*y1-y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1-2.*a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y1))./pi./(1-nu));
e22 = B1.*(1./8.*((1-2.*nu).*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-cos(b).*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))))-2.*y2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-cos(b).*(1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))))-y2.^2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y2-1./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*y2-cos(b).*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y2-1./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*y2)))./pi./(1-nu)+1./4.*((1-2.*nu).*(((2-2.*nu).*cot(b).^2-nu)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-((2-2.*nu).*cot(b).^2+1-2.*nu).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))+(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(y1.*cot(b).*(1-2.*nu-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+nu.*(y3+2.*a)-a+y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(a.*y1.*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2+2.*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-y2.^3./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y2.^3./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))+(1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+(1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2+3.*a.*y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*y1-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(-2.*nu+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*y1.*cot(b)-a)+y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*y1.*cot(b)-a).*y2+2.*y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-y2.^3./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-y2.^3./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-y2.^3./(y1.^2+y2.^2+(y3+2.*a).^2).^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a+2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2-3.*a.*y2.^3./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2))-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b).^2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)+a.*cos(b))+a.*(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).^2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)+a.*cos(b)).*y2-3.*a.*(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*y2+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).^2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*y2+1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(y2.^2.*cos(b).^2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y2.*cos(b).^2+a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*y2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y2)))./pi./(1-nu))+B2.*(1./8.*((2-2.*nu).*(-2./y1./(1+y2.^2./y1.^2)+1./(y1.*cos(b)-y3.*sin(b))./(1+y2.^2./(y1.*cos(b)-y3.*sin(b)).^2)+((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b))+y2.^2./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b))-2.*y2.^2.*(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2.*cos(b))./(1+y2.^2.*(y1.^2+y2.^2+y3.^2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2)+1./(y1.*cos(b)+(y3+2.*a).*sin(b))./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))+y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))-2.*y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*cos(b))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))+y1.*(1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))+y1.*y2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y2-1./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*y2)-(y1.*cos(b)-y3.*sin(b))./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-y2.*(-(y1.*cos(b)-y3.*sin(b))./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y2-(y1.*cos(b)-y3.*sin(b))./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*y2-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y2-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*y2))./pi./(1-nu)+1./4.*((2-2.*nu).*(1-2.*nu).*(-1./y1./(1+y2.^2./y1.^2)+1./(y1.*cos(b)+(y3+2.*a).*sin(b))./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))+y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))-2.*y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*cos(b))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2)).*cot(b).^2+(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((-1+2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)))-(1-2.*nu).*y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((-1+2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+(1-2.*nu).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2.*cot(b)-y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1)-(1-2.*nu).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+(1-2.*nu).*y2.^2.*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+(1-2.*nu).*y2.^2.*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b)-a.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)+3.*a.*y2.^2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2)+(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)))-y2.^2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)))-y2.^2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)))+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)).*y2-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2))+(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*((-2+2.*nu).*cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./cos(b))-y2.^2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*((-2+2.*nu).*cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./cos(b))-y2.^2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*((-2+2.*nu).*cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./cos(b))+y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*y2-2.*a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2./cos(b).*y2))./pi./(1-nu))+B3.*(1./8.*((1-2.*nu).*sin(b).*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-2.*y2.*sin(b).*(1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-y2.^2.*sin(b).*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y2-1./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*y2))./pi./(1-nu)+1./4.*((1-2.*nu).*(-sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1+y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1-(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2)-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(a./(y1.^2+y2.^2+(y3+2.*a).^2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)).*y1+y1.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y2-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2)+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(sin(b).*(cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).*sin(b)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(sin(b).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(1+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*y2-2.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*a.*(y3+2.*a).*y2+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).*sin(b)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*y2+1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(y2.^2.*cos(b).*sin(b)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y2.*cos(b).*sin(b)+a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*y2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y2)))./pi./(1-nu));
e33 = B1.*(1./8.*y2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2).*y3+1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)-cos(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*cos(b).*y3-1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y3-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))-(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y3+4.*a)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))))./pi./(1-nu)+1./4.*((2-2.*nu).*((1-2.*nu).*(-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*sin(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(1./2.*y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*(2.*y3+4.*a)-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*y1)./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2)).*cot(b)-y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)-1./2.*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)+y2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+1./2.*y2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a))+y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*nu./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)+a./(y1.^2+y2.^2+(y3+2.*a).^2))-1./2.*y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*nu./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)+a./(y1.^2+y2.^2+(y3+2.*a).^2)).*(2.*y3+4.*a)+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-2.*nu./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*(2.*y3+4.*a))+y2.*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1-2.*nu-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-1./2.*y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1-2.*nu-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(2.*y3+4.*a)-y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1-2.*nu-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)-a./(y1.^2+y2.^2+(y3+2.*a).^2)+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*(2.*y3+4.*a)))./pi./(1-nu))+B2.*(1./8.*((-1+2.*nu).*sin(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-y1.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2).*y3+1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a))-sin(b).*((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+(y1.*cos(b)-y3.*sin(b)).*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*cos(b).*y3-1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y3-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))-sin(b).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-(y1.*cos(b)+(y3+2.*a).*sin(b)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+1./2.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y3+4.*a)+(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b)))./pi./(1-nu)+1./4.*((-2+2.*nu).*(1-2.*nu).*cot(b).*((1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-cos(b).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))+(2-2.*nu).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+1./2.*(2-2.*nu).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)+(2-2.*nu).*sin(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-(2-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-1./2.*(2-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2))-1./2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2)).*(2.*y3+4.*a)+(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*(2.*y3+4.*a))-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b).*sin(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(sin(b)-(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2)-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))))+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b).*sin(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(sin(b)-(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2)-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*((1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))).*(2.*y3+4.*a)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b)))-1./2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(sin(b)-(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2)-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))).*(2.*y3+4.*a)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2)-(y3+2.*a).*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2)+(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*(2.*y3+4.*a)-sin(b).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-(y1.*cos(b)+(y3+2.*a).*sin(b)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+1./2.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y3+4.*a)+(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b)))))./pi./(1-nu))+B3.*(1./8.*((2-2.*nu).*(y2./(y1.*cos(b)-y3.*sin(b)).^2.*sin(b)./(1+y2.^2./(y1.*cos(b)-y3.*sin(b)).^2)+(y2./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).*y3+y2.*(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2.*y1)./(1+y2.^2.*(y1.^2+y2.^2+y3.^2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2)+y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*sin(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)-(1./2.*y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*(2.*y3+4.*a)-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*y1)./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))+y2.*sin(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*cos(b).*y3-1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y3-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))-(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y3+4.*a)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))))./pi./(1-nu)+1./4.*((2-2.*nu).*(-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*sin(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(1./2.*y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*(2.*y3+4.*a)-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*y1)./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))-(2-2.*nu).*y2.*sin(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-1./2.*(2-2.*nu).*y2.*sin(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)+y2.*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-1./2.*y2.*(y3+a).*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(2.*y3+4.*a)-y2.*(y3+a).*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+y2.*(y3+a).*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*((1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)+a./(y1.^2+y2.^2+(y3+2.*a).^2)-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*(2.*y3+4.*a)))./pi./(1-nu));
e12 = 1./2.*B1.*(1./8.*((2-2.*nu).*(-2./y1./(1+y2.^2./y1.^2)+1./(y1.*cos(b)-y3.*sin(b))./(1+y2.^2./(y1.*cos(b)-y3.*sin(b)).^2)+((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b))+y2.^2./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b))-2.*y2.^2.*(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2.*cos(b))./(1+y2.^2.*(y1.^2+y2.^2+y3.^2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2)+1./(y1.*cos(b)+(y3+2.*a).*sin(b))./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))+y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))-2.*y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*cos(b))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))-y1.*(1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))-y1.*y2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y2-1./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*y2)-cos(b).*(((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-y2.*cos(b).*(1./(y1.^2+y2.^2+y3.^2).*sin(b).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y2-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*y2+1./(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y2-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*y2))./pi./(1-nu)+1./4.*((-2+2.*nu).*(1-2.*nu).*(-1./y1./(1+y2.^2./y1.^2)+1./(y1.*cos(b)+(y3+2.*a).*sin(b))./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))+y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))-2.*y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*cos(b))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2)).*cot(b).^2+(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((1-2.*nu-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)-y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)))-(1-2.*nu).*y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((1-2.*nu-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)-y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+(1-2.*nu).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2.*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1)+(1-2.*nu).*cos(b).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-(1-2.*nu).*y2.^2.*cos(b).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-(1-2.*nu).*y2.^2.*cos(b).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)+a.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)-3.*a.*y2.^2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2)+(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((-1+2.*nu).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2))-y2.^2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((-1+2.*nu).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2))-y2.^2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((-1+2.*nu).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2))+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1-2.*a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y2)+(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+(2-2.*nu).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1).*cos(b))-a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2))-y2.^2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+(2-2.*nu).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1).*cos(b))-a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2))-y2.^2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+(2-2.*nu).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1).*cos(b))-a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2))+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+(2-2.*nu).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1).*cos(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y2.*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2.*cot(b)+(2-2.*nu)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).*y2.*cos(b))+2.*a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y2))./pi./(1-nu))+1./2.*B2.*(1./8.*((-1+2.*nu).*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-cos(b).*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))))+y1.^2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y2-1./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*y2)+(y1.*cos(b)-y3.*sin(b))./(y1.^2+y2.^2+y3.^2).*sin(b).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y2-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*y2+(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y2-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*y2)./pi./(1-nu)+1./4.*((1-2.*nu).*(((2-2.*nu).*cot(b).^2+nu)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-((2-2.*nu).*cot(b).^2+1).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((-1+2.*nu).*y1.*cot(b)+nu.*(y3+2.*a)-a+a.*y1.*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-a.*y1.*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2-y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2)+(1-2.*nu).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*((y1.*cos(b)+(y3+2.*a).*sin(b)).*cos(b)-a.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-(1-2.*nu).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-a./(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).*y2./cos(b)+a.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*y2)+3.*a.*y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*y1-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*y1.*cot(b)+a)-y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*y1.*cot(b)+a).*y2+y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*y2+y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*y2+y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a.*y2+3.*a.*y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*y2)-(y3+a).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(-cos(b).*sin(b)+a.*y1.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+(y3+a).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-3.*a.*y1.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2)./cos(b).*y2+1./(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).*y2.*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))).*y2+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*y2)))./pi./(1-nu))+1./2.*B3.*(1./8.*sin(b).*(((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))./pi./(1-nu)+1./8.*y2.*sin(b).*(1./(y1.^2+y2.^2+y3.^2).*sin(b).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y2-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*y2+1./(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y2-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*y2)./pi./(1-nu)+1./4.*((1-2.*nu).*(1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)-cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+y2.^2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y2.^2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))-(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(a./(y1.^2+y2.^2+(y3+2.*a).^2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))+y2.^2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(a./(y1.^2+y2.^2+(y3+2.*a).^2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y2-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2)+(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-y2.^2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-y2.^2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))+y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2-2.*a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y2))./pi./(1-nu))+1./2.*B1.*(1./8.*((1-2.*nu).*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-cos(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))))-y2.^2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y1-1./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*y1-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*y1-1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*y1-cos(b).*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y1-1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y1-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b)))))./pi./(1-nu)+1./4.*((1-2.*nu).*(((2-2.*nu).*cot(b).^2-nu)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-((2-2.*nu).*cot(b).^2+1-2.*nu).*cos(b).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))+(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(y1.*cot(b).*(1-2.*nu-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+nu.*(y3+2.*a)-a+y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((1-2.*nu-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+a.*y1.^2.*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)-y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1)-(1-2.*nu).*cos(b).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+(1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+(1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1-a.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)+3.*a.*y1.^2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2)-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(-2.*nu+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*y1.*cot(b)-a)+y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*y1.*cot(b)-a).*y1+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1-2.*nu).*cot(b)-y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*y1-y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*y1-y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a.*y1-3.*a.*y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*y1)-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b).^2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)+a.*cos(b))+a.*(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).^2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a))).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)+a.*cos(b)).*y1-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1-2.*nu).*cos(b).*cot(b)+a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)-3.*a.*(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*y1+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).^2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*y1+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(y2.^2.*cos(b).^2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-a.*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)+a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*y1-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y1)))./pi./(1-nu))+1./2.*B2.*(1./8.*((2-2.*nu).*(2.*y2./y1.^2./(1+y2.^2./y1.^2)-y2./(y1.*cos(b)-y3.*sin(b)).^2.*cos(b)./(1+y2.^2./(y1.*cos(b)-y3.*sin(b)).^2)+(y2./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).*y1-y2.*(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2.*(2.*y1.*cos(b)-y3.*sin(b)))./(1+y2.^2.*(y1.^2+y2.^2+y3.^2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2)-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*cos(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*y1-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*(2.*y1.*cos(b)+(y3+2.*a).*sin(b)))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))+y2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))+y1.*y2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y1-1./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*y1-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*y1-1./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*y1)-y2.*(cos(b)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-(y1.*cos(b)-y3.*sin(b))./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y1-(y1.*cos(b)-y3.*sin(b))./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))+cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y1-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))))./pi./(1-nu)+1./4.*((2-2.*nu).*(1-2.*nu).*(y2./y1.^2./(1+y2.^2./y1.^2)-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*cos(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*y1-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*(2.*y1.*cos(b)+(y3+2.*a).*sin(b)))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2)).*cot(b).^2-(1-2.*nu).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((-1+2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1+(1-2.*nu).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1.*cot(b)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))+(1-2.*nu).*y2.*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+(1-2.*nu).*y2.*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*y1+3.*a.*y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*y1-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))).*y1-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))).*y1+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-2.*nu./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)+2.*nu.*y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))+a.*y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1))-y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*((-2+2.*nu).*cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./cos(b)).*y1-y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*((-2+2.*nu).*cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*y1-2.*a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2./cos(b).*y1))./pi./(1-nu))+1./2.*B3.*(1./8.*((1-2.*nu).*sin(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-y2.^2.*sin(b).*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y1-1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y1-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))))./pi./(1-nu)+1./4.*((1-2.*nu).*(-sin(b).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)+cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1)+(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(a./(y1.^2+y2.^2+(y3+2.*a).^2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))-y1.^2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(a./(y1.^2+y2.^2+(y3+2.*a).^2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))+y1.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y1-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1)+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(sin(b).*(cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).*sin(b)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a))).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(sin(b).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1+cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(1+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*y1-2.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*a.*(y3+2.*a).*y1+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).*sin(b)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*y1+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(y2.^2.*cos(b).*sin(b)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-a.*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)+a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*y1-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y1)))./pi./(1-nu));
e13 = 1./2.*B1.*(1./8.*((2-2.*nu).*(y2./(y1.*cos(b)-y3.*sin(b)).^2.*sin(b)./(1+y2.^2./(y1.*cos(b)-y3.*sin(b)).^2)+(y2./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).*y3+y2.*(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2.*y1)./(1+y2.^2.*(y1.^2+y2.^2+y3.^2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2)-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*sin(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(1./2.*y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*(2.*y3+4.*a)-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*y1)./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))-y1.*y2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y3-1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-1)-1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*y3+4.*a)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1))-y2.*cos(b).*(1./(y1.^2+y2.^2+y3.^2).*sin(b).*y3./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y3-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))+1./2./(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).*(2.*y3+4.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y3+4.*a)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))))./pi./(1-nu)+1./4.*((-2+2.*nu).*(1-2.*nu).*(-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*sin(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(1./2.*y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*(2.*y3+4.*a)-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*y1)./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2)).*cot(b).^2-(1-2.*nu).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((1-2.*nu-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)-y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+(1-2.*nu).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(1./2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+1./2.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a))-(1-2.*nu).*y2.*cos(b).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-1./2.*(1-2.*nu).*y2.*cos(b).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2.*cot(b)-3./2.*a.*y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*(2.*y3+4.*a)+y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((-1+2.*nu).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2))-1./2.*y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((-1+2.*nu).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2)).*(2.*y3+4.*a)-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((-1+2.*nu).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)-1./2.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*(2.*y3+4.*a))+y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+(2-2.*nu).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1).*cos(b))-a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2))-1./2.*y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+(2-2.*nu).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1).*cos(b))-a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(2.*y3+4.*a)-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+(2-2.*nu).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1).*cos(b))-a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+(2-2.*nu).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1).*cos(b)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*((1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1).*((1-2.*nu).*cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a).*cot(b)+1./2.*(2-2.*nu)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).*(2.*y3+4.*a).*cos(b))-a.*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2)+a.*(y3+2.*a).*cos(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*(2.*y3+4.*a)))./pi./(1-nu))+1./2.*B2.*(1./8.*((-1+2.*nu).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-1)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-cos(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))))+y1.^2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y3-1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-1)-1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*y3+4.*a)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1))-sin(b).*((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+(y1.*cos(b)-y3.*sin(b))./(y1.^2+y2.^2+y3.^2).*sin(b).*y3./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y3-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))+sin(b).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+1./2.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).*(2.*y3+4.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-1./2.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y3+4.*a)-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b)))./pi./(1-nu)+1./4.*((1-2.*nu).*(((2-2.*nu).*cot(b).^2+nu).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-((2-2.*nu).*cot(b).^2+1).*cos(b).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((-1+2.*nu).*y1.*cot(b)+nu.*(y3+2.*a)-a+a.*y1.*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu-1./2.*a.*y1.*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)-y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)-1./2.*y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a))+(1-2.*nu).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*((y1.*cos(b)+(y3+2.*a).*sin(b)).*cos(b)-a.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-(1-2.*nu).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b).*sin(b)-1./2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).*(2.*y3+4.*a)./cos(b)+1./2.*a.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*(2.*y3+4.*a))-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1.*cot(b)+3./2.*a.*y1.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*(2.*y3+4.*a)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*y1.*cot(b)+a)-y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*y1.*cot(b)+a)-y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*y1.*cot(b)+a).*(2.*y3+4.*a)+1./2.*y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(2.*y3+4.*a)+y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+1./2.*y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a.*(2.*y3+4.*a)+3./2.*a.*y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*(2.*y3+4.*a))+cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-cos(b).*sin(b)+a.*y1.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))))-(y3+a).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(-cos(b).*sin(b)+a.*y1.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)))).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+(y3+a).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*y1-3./2.*a.*y1.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2)./cos(b).*(2.*y3+4.*a)+1./2./(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).*(2.*y3+4.*a).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)))-1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))).*(2.*y3+4.*a)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*(2.*y3+4.*a))))./pi./(1-nu))+1./2.*B3.*(1./8.*y2.*sin(b).*(1./(y1.^2+y2.^2+y3.^2).*sin(b).*y3./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y3-((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))+1./2./(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).*(2.*y3+4.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y3+4.*a)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)-y1)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b)))./pi./(1-nu)+1./4.*((1-2.*nu).*(-y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)-1./2.*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)+y2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+1./2.*y2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a))-y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(a./(y1.^2+y2.^2+(y3+2.*a).^2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))+1./2.*y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(a./(y1.^2+y2.^2+(y3+2.*a).^2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)).*(2.*y3+4.*a)-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-a./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*(2.*y3+4.*a)-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1))+y2.*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-1./2.*y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(2.*y3+4.*a)-y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*((1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)+a./(y1.^2+y2.^2+(y3+2.*a).^2)-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*(2.*y3+4.*a)))./pi./(1-nu))+1./2.*B1.*(1./8.*y2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2).*y1+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1-cos(b).*(1./(y1.^2+y2.^2+y3.^2).*cos(b).*y1./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y1-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))-1./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))))./pi./(1-nu)+1./4.*((2-2.*nu).*((1-2.*nu).*(y2./y1.^2./(1+y2.^2./y1.^2)-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*cos(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*y1-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*(2.*y1.*cos(b)+(y3+2.*a).*sin(b)))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2)).*cot(b)-y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1+y2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+y2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1)-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*nu./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)+a./(y1.^2+y2.^2+(y3+2.*a).^2)).*y1+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-2.*nu./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y1)-y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1-2.*nu-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*y1-y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1-2.*nu-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1+2.*a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y1))./pi./(1-nu))+1./2.*B2.*(1./8.*((-1+2.*nu).*sin(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-1./(y1.^2+y2.^2+y3.^2).^(1./2)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2).*y1+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1)+cos(b).*((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+(y1.*cos(b)-y3.*sin(b))./(y1.^2+y2.^2+y3.^2).*cos(b).*y1./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y1-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))-cos(b).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y1+(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b)))./pi./(1-nu)+1./4.*((-2+2.*nu).*(1-2.*nu).*cot(b).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-cos(b).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-(2-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+(2-2.*nu).*y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+(2-2.*nu).*y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)+(2-2.*nu).*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-(2-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-(2-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1-(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2)).*y1+(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-2.*nu./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)+2.*nu.*y1.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-a./(y1.^2+y2.^2+(y3+2.*a).^2)+2.*a.*y1.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^2)+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b).*sin(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(sin(b)-(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2)-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y1.*cot(b).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))).*y1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b)))-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(sin(b)-(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2)-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))).*y1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-(y3+2.*a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2)+2.*(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y1-cos(b).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y1+(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b)))))./pi./(1-nu))+1./2.*B3.*(1./8.*((2-2.*nu).*(-y2./(y1.*cos(b)-y3.*sin(b)).^2.*cos(b)./(1+y2.^2./(y1.*cos(b)-y3.*sin(b)).^2)+(y2./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).*y1-y2.*(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2.*(2.*y1.*cos(b)-y3.*sin(b)))./(1+y2.^2.*(y1.^2+y2.^2+y3.^2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2)+y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*cos(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)-(y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*y1-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*(2.*y1.*cos(b)+(y3+2.*a).*sin(b)))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))+y2.*sin(b).*(1./(y1.^2+y2.^2+y3.^2).*cos(b).*y1./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y1-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y1-sin(b))-1./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))))./pi./(1-nu)+1./4.*((2-2.*nu).*(y2./y1.^2./(1+y2.^2./y1.^2)-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*cos(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*y1-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*(2.*y1.*cos(b)+(y3+2.*a).*sin(b)))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))-(2-2.*nu).*y2.*sin(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-(2-2.*nu).*y2.*sin(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1-y2.*(y3+a).*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*y1-y2.*(y3+a).*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))+y2.*(y3+a).*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y1-sin(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1-2.*a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y1))./pi./(1-nu));
e23 = 1./2.*B1.*(1./8.*((1-2.*nu).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-1)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3)+(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-cos(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))))-y2.^2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y3-1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-1)-1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*y3+4.*a)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)-cos(b).*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y3-1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))-1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y3+4.*a)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b)))))./pi./(1-nu)+1./4.*((1-2.*nu).*(((2-2.*nu).*cot(b).^2-nu).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-((2-2.*nu).*cot(b).^2+1-2.*nu).*cos(b).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))+(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(y1.*cot(b).*(1-2.*nu-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+nu.*(y3+2.*a)-a+y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)-(1-2.*nu)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(1./2.*a.*y1.*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)+nu-y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)-1./2.*y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a))-(1-2.*nu).*sin(b).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+(1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+1./2.*(1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y1.*cot(b)+3./2.*a.*y1.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*(2.*y3+4.*a)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-2.*nu+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*y1.*cot(b)-a)+y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(-2.*nu+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*y1.*cot(b)-a)+y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*y1.*cot(b)-a).*(2.*y3+4.*a)-1./2.*y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(2.*y3+4.*a)-y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)-1./2.*y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a.*(2.*y3+4.*a)-3./2.*a.*y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*(2.*y3+4.*a))+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b).^2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)+a.*cos(b))+a.*(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).^2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)))-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b).^2-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)+a.*cos(b))+a.*(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).^2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a))).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)+a.*cos(b)).*(2.*y3+4.*a)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1-2.*nu).*sin(b).*cot(b)+a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)+a.*(y3+2.*a).*sin(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)-3./2.*a.*(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*(2.*y3+4.*a)+1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).^2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*(2.*y3+4.*a)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(y2.^2.*cos(b).^2-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-a.*sin(b).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)+1./2.*a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*(2.*y3+4.*a)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b)).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1))))./pi./(1-nu))+1./2.*B2.*(1./8.*((2-2.*nu).*(y2./(y1.*cos(b)-y3.*sin(b)).^2.*sin(b)./(1+y2.^2./(y1.*cos(b)-y3.*sin(b)).^2)+(y2./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).*y3+y2.*(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2.*y1)./(1+y2.^2.*(y1.^2+y2.^2+y3.^2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2)-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*sin(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(1./2.*y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*(2.*y3+4.*a)-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*y1)./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))+y1.*y2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).*y3-1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y3).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-1)-1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*y3+4.*a)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1))-y2.*(-sin(b)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-(y1.*cos(b)-y3.*sin(b))./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y3-(y1.*cos(b)-y3.*sin(b))./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))+sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))-1./2.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y3+4.*a)-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))))./pi./(1-nu)+1./4.*((2-2.*nu).*(1-2.*nu).*(-y2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2.*sin(b)./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+(1./2.*y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).*(2.*y3+4.*a)-y2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*y1)./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2)).*cot(b).^2-(1-2.*nu).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((-1+2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*cot(b)+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+(1-2.*nu).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(-1./2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a).*cot(b)-y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)-1./2.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a))+(1-2.*nu).*y2.*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+1./2.*(1-2.*nu).*y2.*cot(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*(2.*y3+4.*a)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2.*cot(b)+3./2.*a.*y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(5./2).*(2.*y3+4.*a)+y2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)))-1./2.*y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))).*(2.*y3+4.*a)-y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+1./2.*a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)).*(2.*y3+4.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)))+y2.*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*((-2+2.*nu).*cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./cos(b))-1./2.*y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*((-2+2.*nu).*cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./cos(b)).*(2.*y3+4.*a)-y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*((-2+2.*nu).*cos(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./cos(b)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))+y2.*(y3+a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*((1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./cos(b)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-1./2.*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./cos(b).*(2.*y3+4.*a)+a./(y1.^2+y2.^2+(y3+2.*a).^2)./cos(b)-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2./cos(b).*(2.*y3+4.*a)))./pi./(1-nu))+1./2.*B3.*(1./8.*((1-2.*nu).*sin(b).*((1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))+(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-y2.^2.*sin(b).*(-1./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y3-1./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y3-cos(b))-1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(2.*y3+4.*a)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))))./pi./(1-nu)+1./4.*((1-2.*nu).*(-sin(b).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1)+1./2.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)+sin(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-1./2.*(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a))+y1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(a./(y1.^2+y2.^2+(y3+2.*a).^2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a))-1./2.*y1.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(a./(y1.^2+y2.^2+(y3+2.*a).^2)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)).*(2.*y3+4.*a)+y1.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-a./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*(2.*y3+4.*a)-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+1))-1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(sin(b).*(cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).*sin(b)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)))+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(sin(b).*(cos(b)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).*sin(b)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a))).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./2.*sin(b).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*y3+4.*a)+sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-1./2.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(1+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)).*(2.*y3+4.*a)+(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(a./(y1.^2+y2.^2+(y3+2.*a).^2)-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*(2.*y3+4.*a))+1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(y2.^2.*cos(b).*sin(b)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*(2.*y3+4.*a)+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(y2.^2.*cos(b).*sin(b)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*y3+4.*a)+cos(b))-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-a.*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)+1./2.*a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*(2.*y3+4.*a)-a.*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(1./2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*(2.*y3+4.*a)+1))))./pi./(1-nu))+1./2.*B1.*(1./8.*(1./(y1.^2+y2.^2+y3.^2).^(1./2)-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-cos(b).*(((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))))./pi./(1-nu)+1./8.*y2.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2).*y2+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2-cos(b).*(1./(y1.^2+y2.^2+y3.^2).*cos(b).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y2-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y2+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*y2))./pi./(1-nu)+1./4.*((2-2.*nu).*((1-2.*nu).*(-1./y1./(1+y2.^2./y1.^2)+1./(y1.*cos(b)+(y3+2.*a).*sin(b))./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))+y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))-2.*y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*cos(b))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2)).*cot(b)+1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y2.^2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)-cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+y2.^2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y2.^2.*cos(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2))+(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*nu./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)+a./(y1.^2+y2.^2+(y3+2.*a).^2))-y2.^2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(2.*nu./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)+a./(y1.^2+y2.^2+(y3+2.*a).^2))+y2.*(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-2.*nu./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-2.*a./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y2)+(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1-2.*nu-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-y2.^2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1-2.*nu-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-y2.^2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1-2.*nu-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))+y2.*(y3+a).*cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2+2.*a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y2))./pi./(1-nu))+1./2.*B2.*(1./8.*((-1+2.*nu).*sin(b).*(1./(y1.^2+y2.^2+y3.^2).^(1./2).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-y1.*(-1./(y1.^2+y2.^2+y3.^2).^(3./2).*y2+1./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2)+(y1.*cos(b)-y3.*sin(b))./(y1.^2+y2.^2+y3.^2).*cos(b).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y2-(y1.*cos(b)-y3.*sin(b)).*((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*y2-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y2+(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*y2)./pi./(1-nu)+1./4.*((-2+2.*nu).*(1-2.*nu).*cot(b).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))+(2-2.*nu).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2.*(2.*nu+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+(2-2.*nu).*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2-(2-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-(2-2.*nu).*(y1.*cos(b)+(y3+2.*a).*sin(b))./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2-(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((1-2.*nu).*cot(b)-2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a)-a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2)).*y2+(y3+a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*nu.*y1./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)+y3+2.*a).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2+2.*a.*y1./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y2)+(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b).*sin(b)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(sin(b)-(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2)-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-(y3+a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y2.*cot(b).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*((2-2.*nu).*cos(b)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))).*y2+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a).*cot(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(-cos(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2)-a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*(sin(b)-(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2)-(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))).*y2+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*(2.*(y3+2.*a).*(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y2-(y1.*cos(b)+(y3+2.*a).*sin(b))./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y2+(y1.*cos(b)+(y3+2.*a).*sin(b)).*((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*y2)))./pi./(1-nu))+1./2.*B3.*(1./8.*((2-2.*nu).*(1./(y1.*cos(b)-y3.*sin(b))./(1+y2.^2./(y1.*cos(b)-y3.*sin(b)).^2)+((y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b))+y2.^2./(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b))-2.*y2.^2.*(y1.^2+y2.^2+y3.^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2.*cos(b))./(1+y2.^2.*(y1.^2+y2.^2+y3.^2).*sin(b).^2./(y1.*(y1.*cos(b)-y3.*sin(b))+y2.^2.*cos(b)).^2)-1./(y1.*cos(b)+(y3+2.*a).*sin(b))./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))+y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))-2.*y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*cos(b))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))+sin(b).*(((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(1./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)))+y2.*sin(b).*(1./(y1.^2+y2.^2+y3.^2).*cos(b).*y2./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b))-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2).^(3./2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).*y2-((y1.^2+y2.^2+y3.^2).^(1./2).*cos(b)-y3)./(y1.^2+y2.^2+y3.^2)./((y1.^2+y2.^2+y3.^2).^(1./2)-y1.*sin(b)-y3.*cos(b)).^2.*y2-1./(y1.^2+y2.^2+(y3+2.*a).^2).*cos(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b))+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*y2+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*y2))./pi./(1-nu)+1./4.*((2-2.*nu).*(-1./y1./(1+y2.^2./y1.^2)+1./(y1.*cos(b)+(y3+2.*a).*sin(b))./(1+y2.^2./(y1.*cos(b)+(y3+2.*a).*sin(b)).^2)+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))+y2.^2./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b))-2.*y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*sin(b)./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2.*cos(b))./(1+y2.^2.*(y1.^2+y2.^2+(y3+2.*a).^2).*sin(b).^2./(y1.*(y1.*cos(b)+(y3+2.*a).*sin(b))+y2.^2.*cos(b)).^2))+(2-2.*nu).*sin(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-(2-2.*nu).*y2.^2.*sin(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-(2-2.*nu).*y2.^2.*sin(b)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)+(y3+a).*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-y2.^2.*(y3+a).*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))-y2.^2.*(y3+a).*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(1+((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))+a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2))+y2.*(y3+a).*sin(b)./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(1./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b).*y2./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).^2.*(cos(b)+a./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2))./(y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*y2-((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2).*cos(b)+y3+2.*a)./((y1.^2+y2.^2+(y3+2.*a).^2).^(1./2)-y1.*sin(b)+(y3+2.*a).*cos(b)).*a./(y1.^2+y2.^2+(y3+2.*a).^2).^(3./2).*y2-2.*a.*(y3+2.*a)./(y1.^2+y2.^2+(y3+2.*a).^2).^2.*y2))./pi./(1-nu));
end
%====================================================
%% I_CalcStrain.m
function [A1_Ssn,A1_Sdn]=I_CalcStrain(Slip,sSsn,sSdn,dSsn,dSdn)
n=size(Slip,1);
SSlip=Slip(1:n,1:2);
slip=reshape(SSlip,2*n,1);
sdn=[sSsn dSsn; sSdn dSdn];

whos

A1_S=sdn*slip;
A1_Ssn=A1_S(1:n);
A1_Sdn=A1_S(n+1:2*n);
end

%% I_CalcStrainF.m
function [F_Ssn,F_Sdn]=I_CalcStrainF(F_Slip,sSsn,sSdn,dSsn,dSdn)

n=size(F_Slip,1);
FF_Slip=F_Slip(1:n,1:2);
slip=reshape(FF_Slip,2*n,1);
sdn=[sSsn dSsn; sSdn dSdn];
F_S=sdn*slip;
F_Ssn=F_S(1:n);
F_Sdn=F_S(n+1:2*n);
end

%====================================================
%% OutofAsperity01.m
function [A1_Slip]=OutofAsperity01(A1_Ssn,A1_Sdn,Slip,sSsn,sSdn,dSsn,dSdn)
index_out_asp=(Slip(:,1).^2+Slip(:,2).^2)==0;
%O_sSsn=sSsn(~index_asp,~index_asp);
%O_sSdn=sSdn(~index_asp,~index_asp);
%O_dSsn=dSsn(~index_asp,~index_asp);
%O_dSdn=dSdn(~index_asp,~index_asp);
%O_Ssn=A_Ssn(~index_asp);
%O_Sdn=A_Sdn(~index_asp);
%
k=sum(index_out_asp);
d=[A1_Ssn(index_out_asp);A1_Sdn(index_out_asp)];
G=[sSsn(index_out_asp,index_out_asp),dSsn(index_out_asp,index_out_asp);...
   sSdn(index_out_asp,index_out_asp),dSdn(index_out_asp,index_out_asp)];
m=G\d;
A1_Slip=Slip;
A1_Slip(index_out_asp,1:2)=-[m(1:k) m(k+1:end)];

%{
k=0;
for i=1:n
normSlip=norm(Slip(i,:));
 if normSlip==0
  k=k+1;
  k_sSsn(k,:)=sSsn(i,:);
  k_sSdn(k,:)=sSdn(i,:);
  k_dSsn(k,:)=dSsn(i,:);
  k_dSdn(k,:)=dSdn(i,:);
  O_Ssn(k,1)=A_Ssn(i);
  O_Sdn(k,1)=A_Sdn(i);
 else
 end
end

m=0;

for j=1:n
normSlip=norm(Slip(j,:));
 if normSlip==0
  m=m+1;
  O_sSsn(:,m)=k_sSsn(:,j);
  O_sSdn(:,m)=k_sSdn(:,j);
  O_dSsn(:,m)=k_dSsn(:,j);
  O_dSdn(:,m)=k_dSdn(:,j);
 else
 end
end
%}
end
%% OutofAsperity01.m
function [FO_Ssn,FO_Sdn]=OutofAsperity11(F_Ssn,F_Sdn,Slip)
n=length(Slip);
k=0;
for i=1:n
normSlip=norm(Slip(i,:));
 if normSlip==0
  k=k+1;
 
  FO_Ssn(k,1)=F_Ssn(i);
  FO_Sdn(k,1)=F_Sdn(i);
 else
 end
end
end

%====================================================
%% SasaEular2Velo.m
function [vnorm,v,NARFA]=SasaEular2Velo(ll)

Fid=fopen('/home_tmp/sasajima/DATA/Pole_data/PAC-OKH_MORVEL56.dat','r');
Epole=textscan(Fid,'%f %f %f');
fclose(Fid);
Epole=cell2mat(Epole);

%Fid2=fopen('/home_tmp/sasajima/DATA/boso/IA-KA_vll.txt','w');

R(1,1)=6371000;%[m]
%a=6378137.0;
%b=6356752.3;
%aabb=(a.^2)./(b.^2);

diglonP(1,1)=Epole(1,1);%longutitude of Eular pole [digree]%
diglatP(1,1)=Epole(1,2);%latitude of Eular pole [digree]%
digwMa(1,1)=Epole(1,3);%Rotation velocity [digree/Ma]%
radlonP=diglonP./180.*pi;
radlatP(1,1)=diglatP(1,1)./180.*pi;
radwMa=digwMa./180.*pi;

radlonX(:,1)=ll(:,1)./180.*pi;%lon of observed point[radian]
radlatX(:,1)=ll(:,2)./180.*pi;%lat of observed point[radian]
diglonX(:,1)=ll(:,2);%lon of observed poiint [digree]

NARFA(:,1)=radlonX(:,1);
NARFA(:,2)=radlatX(:,1);

CE(1,1)=0.5.*pi-radlatP(1,1);%[radian],scolar
CB(:,1)=0.5.*pi-radlatX(:,1);%[radian],vector
RMD(:,1)=radlonP-radlonX(:,1);%[radian],vector

cCE(1,1)=cos(CE);
sCE(1,1)=sin(CE);

n=length(ll);
v=zeros(n,2);
vnorm=zeros(n,1);

for i=1:n
    cCB(1,1)=cos(CB(i,1));
    sCB(1,1)=sin(CB(i,1));
    cRMD(1,1)=cos(RMD(i,1));
    
    cDLT(1,1)=cCB(1,1).*cCE(1,1)+sCB(1,1).*sCE(1,1).*cRMD(1,1);
    DLT(1,1)=acos(cDLT(1,1));
    
    sDLT(1,1)=sin(DLT(1,1));
    L(1,1)=sDLT(1,1).*R(1,1);%distance between obs point and Eular Vector [m]
    
    vmMa(1,1)=L(1,1).*radwMa(1,1);%velocity[m/Ma]
    vnorm(i,1)=vmMa(1,1)./10000;%velocity[cm/year]
    
    cARFA(1,1)=(cCE(1,1)-(cCB(1,1).*cDLT(1,1)))./(sCB(1,1).*sDLT(1,1));
    ARFA(1,1)=acos(cARFA(1,1));%[radian]
    
    dltdiglon(1,1)=diglonP(1,1)-diglonX(i,1);%[digree]
    
    if dltdiglon>0
        NARFA(i,3)=0.5.*pi-ARFA(1,1);%countour? clock wise from North [radian]
    else
        NARFA(i,3)=0.5.*pi+ARFA(1,1);%countour? clock wise from North[radian]
    end
    
    v(i,1)=vnorm(i,1).*cos(NARFA(i,3));%vlat;N=plus[cm/year]
    v(i,2)=vnorm(i,1).*sin(NARFA(i,3));%vlon;E=plus[cm/year]
    
    %digNARFA(i,1)=NARFA(i,3).*180./pi;%clock wise from North[digree]
    
    %fprintf(Fid2,'%15.8f %15.8f %15.8f %15.8f\n',v(i,1),v(i,2),vnorm(i,1),digNARFA(i,1));
    
end

%fclose(Fid2);
end

%% Sasa_ll2xyz.m
function [trixyzC,trixyz3,sxyz]=trill2trixyz(triC,tri3,sll,ALAT0,ALON0)
% Convert triC, tri3, and sll from LonLat to local XY cordinate.
% Original coded by Sasajima.
% Modified by H. Kimura in 2018/11/13.

nn=size(triC,1);
% mm=size( sll,1);
% trixyz3=zeros(nn,3,3);
% trixyzC=zeros(nn,3);

%*** Convert of sll
[sxyz(:,2),sxyz(:,1)]=PLTXY(sll(:,1),sll(:,2),ALAT0,ALON0);

%*** Convert of tri3
tri3lon=[tri3(:,1);tri3(:,4);tri3(:,7)];
tri3lat=[tri3(:,2);tri3(:,5);tri3(:,8)];
tri3dep=[tri3(:,3);tri3(:,6);tri3(:,9)];
[tmpxyz3(:,2),tmpxyz3(:,1)]=PLTXY(tri3lat,tri3lon,ALAT0,ALON0);
tmpxyz3(:,3)=1e3.*tri3dep;
trixyz3=[tmpxyz3(1:nn,:) tmpxyz3(1+nn:2*nn,:) tmpxyz3(1+2*nn:3*nn,:)];

%*** Convert of triC
[trixyzC(:,2),trixyzC(:,1)]=PLTXY(triC(:,2),triC(:,1),ALAT0,ALON0);
trixyzC(:,3)=1e3.*triC(:,3);

end
%====================================================

%% FSlip2xyll.m
function [xyFSlip]=FSlip2xyll(F_Slip,sitaS)
xyFSlip(:,1)=F_Slip(:,1).*cos(sitaS(:))+F_Slip(:,2).*sin(sitaS(:));
xyFSlip(:,2)=-F_Slip(:,1).*sin(sitaS(:))+F_Slip(:,2).*cos(sitaS(:));
end
%====================================================
%% defineslipQ.m
function [Slip]=defineSlipQ(triC,sitaS)
%A1=Rectangle asperity_1%
%p=minimum X and max Y, q=max X and minimum Y%
%y=+ wa south%

%2011 Tohoku-oki%
A1SV=1.8810;
A1plon=143;
A1plat=41;
A1qlon=145;
A1qlat=38;
A1Slip=0.090;
A1xSlip=-A1Slip.*sin(A1SV);
A1ySlip=A1Slip.*cos(A1SV);

n=size(triC,1);
Slip=zeros(n,2);
define=zeros(n,1);

for i=1:n
  if triC(i,1)<A1plon,Slip(i,:)=[0,0];
  elseif triC(i,1)>A1qlon,Slip(i,:)=[0,0];
  elseif triC(i,2)<A1qlat,Slip(i,:)=[0,0];
  elseif triC(i,2)>A1plat,Slip(i,:)=[0,0];
  else
    A1sSlip=cos(sitaS(i)).*A1xSlip-sin(sitaS(i)).*A1ySlip;
    A1dSlip=sin(sitaS(i)).*A1xSlip+cos(sitaS(i)).*A1ySlip;
    Slip(i,:)=[A1sSlip,A1dSlip];
    define(i,1)=1;
  end
end
%{
%Hokkaido 500year%
A2SV=;
A2plon=145.0;
A2plat=42.2;
A2qlon=146.5;
A2qlat=41.25;
A2Slip=9.2;
A2xSlip=-A2Slip.*sin(A2SV)
A2ySlip=A2Slip.*cos(A2SV)

n=size(triC,1);
Slip=zeros(n,2);
define=zeros(n,1);

 for i=1:n
  if     define(i,1)==1;
  elseif triC(i,1)<A2plon,Slip(i,:)=[0,0];
  elseif triC(i,1)>A2qlon,Slip(i,:)=[0,0];
  elseif triC(i,2)<A2qlat,Slip(i,:)=[0,0];
  elseif triC(i,2)>A2plat,Slip(i,:)=[0,0];
  else
   A2sSlip=cos(sitaS(i)).*A2xSlip-sin(sitaS(i)).*A2ySlip;
   A2dSlip=sin(sitaS(i)).*A2xSlip+cos(sitaS(i)).*A2ySlip;
   Slip(i,:)=[A2sSlip,A2dSlip];
   define(i,1)=[1];
  end
 end

%1896 Meiji Sanriku-oki%
A3SV=;
A3plon=143.7;
A3plat=39.6;
A3qlon=144.15;
A3qlat=38.8;
A3Slip=8.85;
A3xSlip=-A3Slip.*sin(A3SV)
A3ySlip=A3Slip.*cos(A3SV)

n=size(triC,1);
Slip=zeros(n,2);
define=zeros(n,1);

 for i=1:n
  if     define(i,1)==1;
  elseif triC(i,1)<A3plon,Slip(i,:)=[0,0];
  elseif triC(i,1)>A3qlon,Slip(i,:)=[0,0];
  elseif triC(i,2)<A3qlat,Slip(i,:)=[0,0];
  elseif triC(i,2)>A3plat,Slip(i,:)=[0,0];
  else
   A3sSlip=cos(sitaS(i)).*A3xSlip-sin(sitaS(i)).*A3ySlip;
   A3dSlip=sin(sitaS(i)).*A3xSlip+cos(sitaS(i)).*A3ySlip;
   Slip(i,:)=[A3sSlip,A3dSlip];
   define(i,1)=[1];
  end
 end

%2003 Tokachi-oki%
A4SV=;
A4plon=143.65;
A4plat=42.2;
A4qlon=144.2;
A4qlat=41.7;
A4Slip=8.5;
A4xSlip=-A4Slip.*sin(A3SV)
A4ySlip=A4Slip.*cos(A3SV)

n=size(triC,1);
Slip=zeros(n,2);
define=zeros(n,1);

 for i=1:n
  if     define(i,1)==1;
  elseif triC(i,1)<A4plon,Slip(i,:)=[0,0];
  elseif triC(i,1)>A4qlon,Slip(i,:)=[0,0];
  elseif triC(i,2)<A4qlat,Slip(i,:)=[0,0];
  elseif triC(i,2)>A4plat,Slip(i,:)=[0,0];
  else
   A4sSlip=cos(sitaS(i)).*A4xSlip-sin(sitaS(i)).*A4ySlip;
   A4dSlip=sin(sitaS(i)).*A4xSlip+cos(sitaS(i)).*A4ySlip;
   Slip(i,:)=[A4sSlip,A4dSlip];
   define(i,1)=[1];
  end
 end 

%}
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

%% ll2xy.m
function [xy]=ll2xy(ll)
%-------------------
%  PLTXY TRANSFORMS (ALON,ALAT) TO (X,Y)
%  WHEN ICORD.NE.0  PLTXY MAKES NO CHANGE IN
%  TRANSFORMATION BETWEEN (X,Y) AND (ALAT,ALONG).
%-------------------
A=6.378160e3;
E2=6.6944541e-3;
E12=6.7395719e-3;
D=5.72958e1;
RD=1.0/D;
ALON0=144;
ALAT0=40;

mm=length(ll);
xy=zeros(mm,2);

for i=1:mm
    ALON=ll(i,1);
    ALAT=ll(i,2);
    
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
    xy(i,1:2)=[X,Y];
end

end
%====================================================

%% loadMAT.m
function [sSsn,sSdn,dSsn,dSdn]=loadMAT(trixyzC)

rootname='/home_tmp/sasajima/DATA/GreenF/PACTohoku2PACTohoku/';
ss='sSsn/sSsn';
sd='sSdn/sSdn';
ds='dSsn/dSsn';
dd='dSdn/dSdn';
extension='.dat';

n=length(trixyzC);

for i=1:n
    w=num2str(i);
    
    filename1= [rootname,ss,w,extension];
    filename2= [rootname,sd,w,extension];
    filename3= [rootname,ds,w,extension];
    filename4= [rootname,dd,w,extension];
    
    load(filename1,'sSsni','-mat');
    load(filename2,'sSdni','-mat');
    load(filename3,'dSsni','-mat');
    load(filename4,'dSdni','-mat');
    
    sSsn(:,i)=sSsni(:,1);
    sSdn(:,i)=sSdni(:,1);
    dSsn(:,i)=dSsni(:,1);
    dSdn(:,i)=dSdni(:,1);
end
end
%% loadMAT2.m
function [sUxyz,dUxyz]=loadMAT2(xyz,trixyz3)

rootname='/home_tmp/sasajima/DATA/GreenF/PAC2test/';
sU='sU';
dU='dU';
extension='.dat';

m=size(xyz,1);
n=size(trixyz3,1);

for i=1:n
    w=num2str(i);
    
    filename1= [rootname,sU,w,extension];
    filename2= [rootname,dU,w,extension];
    
    load(filename1,'sUxyzi','-mat');
    load(filename2,'dUxyzi','-mat');
    
    sUxyz(:,i,1:3)=sSsni(:,1:3);
    dUxyz(:,i,1:3)=dSdni(:,1);
end
end
%% make Green's function of displacement
function [gu]=makeGreenDisp(xyz,trixyz3)
% This function make Green's function of displacement (gu).
% Input
%  xyz     : site location.
%  trixyz3 : trimesh coordianate.
% Output
%  gu.st   : surface displacement due to strike slip on a fault.
%  gu.dp   : surface displacement due to dip slip on a fault.
%  gu.ts   : surface displacement due to tensile slip on a fault.
%  *Each matrix (gu.*) contain 3 x NOBS by NFLT elements.

pr  =0.25; % Poisson's ratio
nobs=size(xyz,1);
nflt=size(trixyz3,1);
sx  =xyz(:,1);
sy  =xyz(:,2);
sz  =zeros(nobs,1);

for n=1:nflt
    trix=[trixyz3(n,1),trixyz3(n,4),trixyz3(n,7)];
    triy=[trixyz3(n,2),trixyz3(n,5),trixyz3(n,8)];
    triz=[trixyz3(n,3),trixyz3(n,6),trixyz3(n,9)];
%     [U] = CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds);
    [U] = CalcTriDisps(sx, sy, sz, trix, triy, triz, pr, 1, 0, 0); % Strike
    gu.st(1:3:3*nobs,n)= U.x;
    gu.st(2:3:3*nobs,n)= U.y;
    gu.st(3:3:3*nobs,n)=-U.z;
    [U] = CalcTriDisps(sx, sy, sz, trix, triy, triz, pr, 0, 1, 0); % Tensile
    gu.ts(1:3:3*nobs,n)= U.x;
    gu.ts(2:3:3*nobs,n)= U.y;
    gu.ts(3:3:3*nobs,n)=-U.z;
    [U] = CalcTriDisps(sx, sy, sz, trix, triy, triz, pr, 0, 0, 1); % Dip
    gu.dp(1:3:3*nobs,n)= U.x;
    gu.dp(2:3:3*nobs,n)= U.y;
    gu.dp(3:3:3*nobs,n)=-U.z;
end

end
%% make Green's function of strain
function [gs]=makeGreenStrain(trixyzC,trixyz3,sitaS,sitaD,normVec)
% This function make Green's function of strain (gs).
% Input
%  trixyzC : center of trimesh.
%  trixyz3 : trimesh coordianate.
%  sitaS   : angle of fault strike from X-axis.
%  sitaD   : angle of fault dip from XY-plane.
%  normVec : normalized normal vectors of faults.
% Output
%  gs.st   : strain due to strike slip on a fault.
%  gs.dp   : strain due to dip slip on a fault.
%  gs.ts   : strain due to tensile slip on a fault.
%  *Each matrix (gs.*) contain 6 x NFLT by NFLT elements.
pr  =0.25; % Poisson's ratio
nflt=size(trixyz3,1);
sx=trixyzC(:,1)+normVec(:,1);
sy=trixyzC(:,2)+normVec(:,2);
sz=trixyzC(:,3)+normVec(:,3);

for n=1:nflt
    trix=[trixyz3(n,1),trixyz3(n,4),trixyz3(n,7)];
    triy=[trixyz3(n,2),trixyz3(n,5),trixyz3(n,8)];
    triz=[trixyz3(n,3),trixyz3(n,6),trixyz3(n,9)];
%     [S] = CalcTriStrains(sx, sy, sz, x, y, z, pr, ss, ts, ds);
    [S] = CalcTriStrains(sx, sy, sz, trix, triy, triz, pr, 1, 0, 0); % Strike
    gs.st(1:6:6*nflt,n)=S.xx;
    gs.st(2:6:6*nflt,n)=S.xy;
    gs.st(3:6:6*nflt,n)=S.xz;
    gs.st(4:6:6*nflt,n)=S.yy;
    gs.st(5:6:6*nflt,n)=S.yz;
    gs.st(6:6:6*nflt,n)=S.zz;
    [S] = CalcTriStrains(sx, sy, sz, trix, triy, triz, pr, 0, 1, 0); % Tensile
    gs.ts(1:6:6*nflt,n)=S.xx;
    gs.ts(2:6:6*nflt,n)=S.xy;
    gs.ts(3:6:6*nflt,n)=S.xz;
    gs.ts(4:6:6*nflt,n)=S.yy;
    gs.ts(5:6:6*nflt,n)=S.yz;
    gs.ts(6:6:6*nflt,n)=S.zz;
    [S] = CalcTriStrains(sx, sy, sz, trix, triy, triz, pr, 0, 0, 1); % Dip
    gs.dp(1:6:6*nflt,n)=S.xx;
    gs.dp(2:6:6*nflt,n)=S.xy;
    gs.dp(3:6:6*nflt,n)=S.xz;
    gs.dp(4:6:6*nflt,n)=S.yy;
    gs.dp(5:6:6*nflt,n)=S.yz;
    gs.dp(6:6:6*nflt,n)=S.zz;
end
[gs]=trans_xyz2strdip(gs,sitaS,sitaD);
end
%% Transform tensor from xyz to strike-dip
function [gs]=trans_xyz2strdip(gs,sitaS,sitaD)
% This function transforms a strain tensor from xyz to fault strike-dip.
% 
% Output
% gs.stst : strain of strike direction on the fault due to strike slip.
% gs.stdp : strain of dip direction on the fault due to strike slip.
% gs.stts : strain of tensile direction on the fault due to strike slip.
% gs.dpst : strain of strike direction on the fault due to dip slip.
% gs.dpdp : strain of dip direction on the fault due to dip slip.
% gs.dpts : strain of tensile direction on the fault due to dip slip.
% gs.tsst : strain of strike direction on the fault due to tensile slip.
% gs.tsdp : strain of dip direction on the fault due to tensile slip.
% gs.tsts : strain of tensile direction on the fault due to tensile slip.
% 
% Note that strain of strike direction on the fault corresponds to ezx',
% strain of dip direction on the fault corresponds to ezy', and
% strain of tensile direction on the fault corresponds to ezz'.

% c1=cos(sitaS);
% s1=sin(sitaS);
% c2=cos(sitaD);
% s2=sin(sitaD);
c1=repmat(cos(sitaS)',1,size(sitaS,2));
s1=repmat(sin(sitaS)',1,size(sitaS,2));
c2=repmat(cos(sitaD)',1,size(sitaD,2));
s2=repmat(sin(sitaD)',1,size(sitaD,2));

[sst,sdp,sts]=calctrans(gs.st,c1,s1,c2,s2); % response to strike slip
gs.stst=sst;
gs.stdp=sdp;
gs.stts=sts;
[sst,sdp,sts]=calctrans(gs.ts,c1,s1,c2,s2); % response to tensile slip
gs.tsst=sst;
gs.tsdp=sdp;
gs.tsts=sts;
[sst,sdp,sts]=calctrans(gs.dp,c1,s1,c2,s2); % response to dip slip
gs.dpst=sst;
gs.dpdp=sdp;
gs.dpts=sts;

end
%%
function [sst,sdp,sts]=calctrans(sxyz,c1,s1,c2,s2)
% E_stdp = Tx * Tz * E_xyz * Tz' * Tx'
% <->
% |exx' exy' exz'| |1   0  0| | c1 s1 0| |exx exy exz| | c1 s1 0|T |1   0  0|T
% |eyx' eyy' eyz'|=|0  c2 s2|*|-s1 c1 0|*|eyx eyy eyz|*|-s1 c1 0| *|0  c2 s2|
% |ezx' ezy' ezz'| |0 -s2 c2| |  0  0 1| |ezx ezy ezz| |  0  0 1|  |0 -s2 c2|
% where 
% c1 : cos(strike)
% s1 : sin(strike)
% c2 : cos(dip)
% s2 : sin(dip)
% Note strike is the counter clock wise angle from x-axis, dip is angle
% from xy-plane. "T" indicates the transpose of matrix.

c1_2 = c1.^2;
c2_2 = c2.^2;
s1_2 = s1.^2;
s2_2 = s2.^2;
c1s1 = c1.*s1;
c2s2 = c2.*s2;
% ezx' = -s2*( c1*s1*( eyy - exx ) + ( c1^2 - c2^2 )*exy ) + c2( c1*exz + s1*eyz );
% ezy' = c2*s2*( ezz - s1^2*exx + 2*c1*s1*exy - c1^2*eyy ) + ( c2^2 - s2^2 )*( c1*eyz - s1*exz );
% ezz' = s2^2*( s1^2*exx -2*c1*s1*exy + c1^2*eyy ) + 2*c2*s2*( s1*exz - c1*eyz ) + c2^2*ezz;
exx = sxyz(1:6:end,:);
exy = sxyz(2:6:end,:);
exz = sxyz(3:6:end,:);
eyy = sxyz(4:6:end,:);
eyz = sxyz(5:6:end,:);
ezz = sxyz(6:6:end,:);
sst = -s2.*( c1s1.*( eyy - exx ) + ( c1_2 - c2_2 ).*exy )...
         + c2.*( c1.*exz + s1.*eyz );
sdp =  c2.*s2.*( ezz - s1_2.*exx + 2.*c1s1.*exy - c1_2.*eyy )...
         + ( c2_2 - s2_2 ).*( c1.*eyz - s1.*exz );
sts = s2_2.*( s1_2.*exx -2*c1s1.*exy + c1_2.*eyy )...
         + 2.*c2s2.*( s1.*exz - c1.*eyz )...
         + c2_2.*ezz;
end
%% makeG_d_O.m
function[di]=makeG_d_O(xyz,trixyz3)

%change the value%
pr=0.25;
ss= 0;
ds= 1;
ts= 0;
%==================

m=size(xyz,1);
n=size(trixyz3,1);

sx=xyz(:,1);
sy=xyz(:,2);
sz=zeros(m,1);

rootname='/home_tmp/sasajima/DATA/GreenF/PAC2test/';
dU='dU';
extension='.dat';

for i=1:n%numbers of fault%
    
    di=i;
    w=num2str(i);
    
    x=[trixyz3(i,1),trixyz3(i,4),trixyz3(i,7)];
    y=[trixyz3(i,2),trixyz3(i,5),trixyz3(i,8)];
    z=[trixyz3(i,3),trixyz3(i,6),trixyz3(i,9)];
    
    [U] = CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds);
    dUxyz(:,1)=U.x;
    dUxyz(:,2)=U.y;
    dUxyz(:,3)=-U.z;
    
    %[S] = CalcTriStrains(sx, sy, sz, x, y, z, pr, ss, ts, ds);
    %Sxx(:,i)=S.xx;
    %Sxy(:,i)=S.xy;
    %Sxz(:,i)=S.xz;
    %Syy(:,i)=S.yy;
    %Syz(:,i)=S.yz;
    %Szz(:,i)=S.zz;
    
    filename1= [rootname,dU,w,extension];
    
    dUxyzi(:,1:3)=dUxyz(:,1:3);
    
    save(filename1,'dUxyzi','-v7.3');
end
end
%% makeG_d_Q.m
function[di]=makeG_d_Q(trixyzC,trixyz3,sitaS,sitaD,normVec)

%change the value%
pr=0.25;
ss= 0;
ds= 1;
ts= 0;
%==================

n=length(trixyzC);

sx=trixyzC(:,1)+normVec(:,1);
sy=trixyzC(:,2)+normVec(:,2);
sz=trixyzC(:,3)+normVec(:,3);
c1=cos(sitaS);
s1=sin(sitaS);
c2=cos(sitaD);
s2=sin(sitaD);

rootname='/home_tmp/sasajima/DATA/GreenF/PACtest2PACtest/';
dds='dSsn/dSsn';
ddn='dSdn/dSdn';
extension='.dat';

for i=1:n%numbers of fault%
    
    di=i;
    w=num2str(i);
    
    x=[trixyz3(i,1),trixyz3(i,4),trixyz3(i,7)];
    y=[trixyz3(i,2),trixyz3(i,5),trixyz3(i,8)];
    z=[trixyz3(i,3),trixyz3(i,6),trixyz3(i,9)];
    
    %[U] = CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds);
    %dUxyz(:,i,1)=U.x;
    %dUxyz(:,i,2)=U.y;
    %dUxyz(:,i,3)=-U.z;
    
    [S] = CalcTriStrains(sx, sy, sz, x, y, z, pr, ss, ts, ds);
    Sxx(:,1)=S.xx;
    Sxy(:,1)=S.xy;
    Sxz(:,1)=S.xz;
    Syy(:,1)=S.yy;
    Syz(:,1)=S.yz;
    Szz(:,1)=S.zz;
    
    %dSss(:,i)=c1^2*Sxx+2*c1*s1*Sxy+s1^2*Syy;
    %dSsd(:,i)=c2*(c1*s1*(Syy-Sxx)+(c1^2-s1^2)*Sxy)+s2*(c1*Sxz+s1*Syz);
    dSsn(:,i)=-s2(i).*(c1(i).*s1(i).*(Syy(:,1)-Sxx(:,1))+(c1(i).^2-s1(i).^2).*Sxy(:,1))...
             + c2(i).*(c1(i).*Sxz(:,1)+s1(i).*Syz(:,1));
    %dSdd(:,i)=c2^2*(s1^2*Sxx-2*s1*c1*Sxy+c1^2*Syy)+2*c2*s2*(c1*Syz-s1*Sxz)+s2^2*Szz;
    dSdn(:,i)=c2(i).*s2(i).*(Szz(:,1)-(s1(i)).^2.*Sxx(:,1)+2.*s1(i).*c1(i).*Sxy(:,1)-(c1(i)).^2.*Syy(:,1))+((c2(i)).^2-(s2(i)).^2).*(c1(i).*Syz(:,1)-s1(i).*Sxz(:,1));
    %dSnn(:,i)=s2^2*(s1^2*Sxx-2*s1*c1*Sxy+c1^2*Syy)-2*c2*s2*(c1*Syz-s1*Sxz)+c2^2*Szz;
    
    filename1= [rootname,dds,w,extension];
    filename2= [rootname,ddn,w,extension];
    
    dSsni(:,1)=dSsn(:,1);
    dSdni(:,1)=dSdn(:,1);
    
    save(filename1,'dSsni','-v7.3');
    save(filename2,'dSdni','-v7.3');
end
end
%% makeG_s_O.m
function[si]=makeG_s_O(xyz,trixyz3)

%change the value%
pr=0.25;
ss= 1;
ds= 0;
ts= 0;
%==================

m=size(xyz,1);
n=size(trixyz3,1);

sx=xyz(:,1);
sy=xyz(:,2);
sz=zeros(m,1);

rootname='/home_tmp/sasajima/DATA/GreenF/PAC2test/';
sU='sU';
extension='.dat';

for i=1:n%numbers of fault%
    
    si=i;
    w=num2str(i);
    
    x=[trixyz3(i,1),trixyz3(i,4),trixyz3(i,7)];
    y=[trixyz3(i,2),trixyz3(i,5),trixyz3(i,8)];
    z=[trixyz3(i,3),trixyz3(i,6),trixyz3(i,9)];
    
    [U] = CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds);
    sUxyz(:,1)=U.x;
    sUxyz(:,2)=U.y;
    sUxyz(:,3)=-U.z;
    
    %[S] = CalcTriStrains(sx, sy, sz, x, y, z, pr, ss, ts, ds);
    %Sxx(:,i)=S.xx;
    %Sxy(:,i)=S.xy;
    %Sxz(:,i)=S.xz;
    %Syy(:,i)=S.yy;
    %Syz(:,i)=S.yz;
    %Szz(:,i)=S.zz;
    
    filename1= [rootname,sU,w,extension];
    
    sUxyzi(:,1:3)=sUxyz(:,1:3);
    
    save(filename1,'sUxyzi','-v7.3');
end
end
%% makeG_s_Q.m
function[si]=makeG_s_Q(trixyzC,trixyz3,sitaS,sitaD,normVec)

%change the value%
pr=0.25;
ss= 1;
ds= 0;
ts= 0;
%==================

n=length(trixyzC);

sx=trixyzC(:,1)+normVec(:,1);
sy=trixyzC(:,2)+normVec(:,2);
sz=trixyzC(:,3)+normVec(:,3);
c1=cos(sitaS);
s1=sin(sitaS);
c2=cos(sitaD);
s2=sin(sitaD);


rootname='/home_tmp/sasajima/DATA/GreenF/PACtest2PACtest/';
ssn='sSsn/sSsn';
sdn='sSdn/sSdn';
extension='.dat';

for i=1:n%numbers of fault%
    
    si=i;
    w=num2str(i);
    
    x=[trixyz3(i,1),trixyz3(i,4),trixyz3(i,7)];
    y=[trixyz3(i,2),trixyz3(i,5),trixyz3(i,8)];
    z=[trixyz3(i,3),trixyz3(i,6),trixyz3(i,9)];
    
    %[U] = CalcTriDisps(sx, sy, sz, x, y, z, pr, ss, ts, ds);
    %sUxyz(:,i,1)=U.x;
    %sUxyz(:,i,2)=U.y;
    %sUxyz(:,i,3)=-U.z;
    [S] = CalcTriStrains(sx, sy, sz, x, y, z, pr, ss, ts, ds);
    Sxx(:,1)=S.xx;
    Sxy(:,1)=S.xy;
    Sxz(:,1)=S.xz;
    Syy(:,1)=S.yy;
    Syz(:,1)=S.yz;
    Szz(:,1)=S.zz;
    
    %sSss(:,i)=c1^2*Sxx+2*c1*s1*Sxy+s1^2*Syy;
    %sSsd(:,i)=c2*(c1*s1*(Syy-Sxx)+(c1^2-s1^2)*Sxy)+s2*(c1*Sxz+s1*Syz);
    sSsn(:,1)=-s2(i).*(c1(i).*s1(i).*(Syy(:,1)-Sxx(:,1))+((c1(i)).^2-(s1(i)).^2).*Sxy(:,1))+c2(i).*(c1(i).*Sxz(:,1)+s1(i).*Syz(:,1));
    %sSdd(:,i)=c2^2*(s1^2*Sxx-2*s1*c1*Sxy+c1^2*Syy)+2*c2*s2*(c1*Syz-s1*Sxz)+s2^2*Szz;
    sSdn(:,1)=c2(i).*s2(i).*(Szz(:,1)-(s1(i)).^2.*Sxx(:,1)+2.*s1(i).*c1(i).*Sxy(:,1)-(c1(i)).^2.*Syy(:,1))+((c2(i)).^2-(s2(i)).^2).*(c1(i).*Syz(:,1)-s1(i).*Sxz(:,1));
    %sSnn(:,i)=s2^2*(s1^2*Sxx-2*s1*c1*Sxy+c1^2*Syy)-2*c2*s2*(c1*Syz-s1*Sxz)+c2^2*Szz;
    
    filename1= [rootname,ssn,w,extension];
    filename2= [rootname,sdn,w,extension];
    
    sSsni(:,1)=sSsn(:,1);
    sSdni(:,1)=sSdn(:,1);
    
    save(filename1,'sSsni','-v7.3');
    save(filename2,'sSdni','-v7.3');
end
end
%% make_test_trill.m
%====================================================
function [triC,tri3,tri,sll]=make_test_trill
%by Ryohei Sasajima 2013/12/23
%====================================================
%====================================================
Fid0=fopen('/home/sasajima/Dropbox/yellow/PACdepth201312.txt','r');
dep_main=textscan(Fid0,'%f %f %f');
fclose(Fid0);
dep_main=cell2mat(dep_main);
%====================================================
Fid00=fopen('/home/sasajima/Dropbox/yellow/depth0_1402.txt','r');
dep_sub=textscan(Fid00,'%f %f %f');
fclose(Fid00);
dep_sub=cell2mat(dep_sub);
%===================================================
%====================================================
Fid1=fopen('/home/sasajima/Dropbox/yellow/PAC_blue.txt','r');
bound=textscan(Fid1,'%f %f');
fclose(Fid1);
bound=cell2mat(bound);

Fid2=fopen('/home/sasajima/Dropbox/yellow/PAC_blue.txt','r');
PACblue=textscan(Fid2,'%f %f');
fclose(Fid2);
PACblue=cell2mat(PACblue);

Fid3=fopen('/home/sasajima/Dropbox/yellow/PAC_green.txt','r');
PACgreen=textscan(Fid3,'%f %f');
fclose(Fid3);
PACgreen=cell2mat(PACgreen);

Fid4=fopen('/home_tmp/sasajima/DATA/blue.dat','r');
mblue=textscan(Fid4,'%f %f');
fclose(Fid4);
mblue=cell2mat(mblue);
bb=size(mblue,1);

Fid5=fopen('/home_tmp/sasajima/DATA/green.dat','r');
mgreen=textscan(Fid5,'%f %f');
fclose(Fid5);
mgreen=cell2mat(mgreen);
gg=size(mgreen,1);

Fid7=fopen('/home/sasajima/Dropbox/yellow/PAC_red.txt','r');
PACred=textscan(Fid7,'%f %f');
fclose(Fid7);
PACred=cell2mat(PACred);

Fid8=fopen('/home/sasajima/Dropbox/yellow/PAC_black.txt','r');
PACblack=textscan(Fid8,'%f %f');
fclose(Fid8);
PACblack=cell2mat(PACblack);

Fid9=fopen('/home_tmp/sasajima/DATA/red.dat','r');
mred=textscan(Fid9,'%f %f');
fclose(Fid9);
mred=cell2mat(mred);
rr=size(mred,1);

Fid10=fopen('/home_tmp/sasajima/DATA/black.dat','r');
mblack=textscan(Fid10,'%f %f');
fclose(Fid10);
mblack=cell2mat(mblack);
blbl=size(mblack,1);

Fid11=fopen('/home_tmp/sasajima/DATA/m_all.dat','r');
mall=textscan(Fid11,'%f %f');
fclose(Fid11);
mall=cell2mat(mall);
alal=size(mall,1);


s=0;
for bl=1:blbl
    Blue=inpolygon(mblack(bl,1),mblack(bl,2),PACblue(:,1),PACblue(:,2));
    if Blue==0
    else
        s=s+1;
        slon(s,1)=mblack(bl,1);
        slat(s,1)=mblack(bl,2);
    end
end

%{

 for b=1:bb;

 Blue=inpolygon(mblue(b,1),mblue(b,2),PACblue(:,1),PACblue(:,2));
 BRed=inpolygon(mblue(b,1),mblue(b,2),PACred(:,1),PACred(:,2));

  if (Blue==1)&&(BRed==0)
   s=s+1;
   slon(s,1)=mblue(b,1);
   slat(s,1)=mblue(b,2);
  else
  end
 end

 sblue=s-sred

 for g=1:gg;

  Green=inpolygon(mgreen(g,1),mgreen(g,2),PACgreen(:,1),PACgreen(:,2));
 GBlue=inpolygon(mgreen(g,1),mgreen(g,2),PACblue(:,1),PACblue(:,2));

  if (Green==1)&&(GBlue==0)
   s=s+1;
   slon(s,1)=mgreen(g,1);
   slat(s,1)=mgreen(g,2);
  else
  end
 end

 sgreen=s-sred-sblue

 for bl=1:blbl;

  Black=inpolygon(mblack(bl,1),mblack(bl,2),PACblack(:,1),PACblack(:,2));
 BlGreen=inpolygon(mblack(bl,1),mnlack(bl,2),PACgreen(:,1),PACgreen(:,2));

  if (Black==1)&&(BlGreen==0)
   s=s+1;
   slon(s,1)=mblack(bl,1);
   slat(s,1)=mblack(bl,2);
  else
  end
 end

 sblack=s-sred-sblue-sgreen


 for al=1:alal;

  All=inpolygon(mall(al,1),mall(al,2),bound(:,1),bound(:,2));
 AlBlack=inpolygon(mall(al,1),mall(al,2),PACblack(:,1),PACblack(:,2));

  if (All==1)&&(AlBlack==0)
   s=s+1;
   slon(s,1)=mall(al,1);
   slat(s,1)=mall(al,2);
  else
  end
 end

 soutblack=s-sred-sblue-sgreen-sblack
 sall=s
 pause

%}

sll(:,1:2)=[slon(:,1),slat(:,1)];

%ll=[s.lon,s.lat];
%plot(s.lon,s.lat,'.')
%====================================================
tri = delaunay(slon,slat);
%====================================================
ntri=length(tri);
E=scatteredInterpolant(dep_main(:,1),dep_main(:,2),dep_main(:,3),'natural');
F=scatteredInterpolant(dep_sub(:,1),dep_sub(:,2),dep_sub(:,3),'natural');

tri3=zeros(ntri,3,3);
triC=zeros(ntri,3);
nn=0;
for n=1:ntri
  lon1=slon(tri(n,1));
  lat1=slat(tri(n,1));
  dep1=-F(lon1,lat1)-E(lon1,lat1);
  lon2=slon(tri(n,2));
  lat2=slat(tri(n,2));
  dep2=-F(lon2,lat2)-E(lon2,lat2);
  lon3=slon(tri(n,3));
  lat3=slat(tri(n,3));
  dep3=-F(lon3,lat3)-E(lon3,lat3);
  
  glon=(lon1+lon2+lon3)/3;
  glat=(lat1+lat2+lat3)/3;
  gdep=(dep1+dep2+dep3)/3;
  
  IN=inpolygon(glon,glat,bound(:,1),bound(:,2));
  
  if gdep==0
  elseif IN==0
  else
    nn=nn+1;
    triC(nn,1:3)=[glon,glat,gdep];
    tri3(nn,1)=lon1;
    tri3(nn,2)=lat1;
    tri3(nn,3)=dep1;
    tri3(nn,4)=lon2;
    tri3(nn,5)=lat2;
    tri3(nn,6)=dep2;
    tri3(nn,7)=lon3;
    tri3(nn,8)=lat3;
    tri3(nn,9)=dep3;
  end
end
tri3(nn+1:ntri,:)=[];
triC(nn+1:ntri,:)=[];

Fid6=fopen('/home_tmp/sasajima/DATA/PAC_tri3.txt','w');
for w=1:nn
  for u=1:3
    fprintf(Fid6,'%10.4f %9.4f %10.4f\n',tri3(w,1,u),tri3(w,2,u),tri3(w,3,u));
  end
end
fclose(Fid6);

end
%% make_trill.m
%====================================================
function [triC,tri3,tri,sll]=make_trill
%by Ryohei Sasajima 2013/12/23
%====================================================
%====================================================
Fid0=fopen('/home/sasajima/Dropbox/yellow/PACdepth201312.txt','r');
dep_main=textscan(Fid0,'%f %f %f');
fclose(Fid0);
dep_main=cell2mat(dep_main);
%====================================================
Fid00=fopen('/home/sasajima/Dropbox/yellow/depth0_1402.txt','r');
dep_sub=textscan(Fid00,'%f %f %f');
fclose(Fid00);
dep_sub=cell2mat(dep_sub);
%===================================================
%====================================================
Fid1=fopen('/home/sasajima/Dropbox/yellow/PAC_all.txt','r');
bound=textscan(Fid1,'%f %f');
fclose(Fid1);
bound=cell2mat(bound);

Fid2=fopen('/home/sasajima/Dropbox/yellow/PAC_blue.txt','r');
PACblue=textscan(Fid2,'%f %f');
fclose(Fid2);
PACblue=cell2mat(PACblue);

Fid3=fopen('/home/sasajima/Dropbox/yellow/PAC_green.txt','r');
PACgreen=textscan(Fid3,'%f %f');
fclose(Fid3);
PACgreen=cell2mat(PACgreen);

Fid4=fopen('/home_tmp/sasajima/DATA/blue.dat','r');
mblue=textscan(Fid4,'%f %f');
fclose(Fid4);
mblue=cell2mat(mblue);
bb=size(mblue,1);

Fid5=fopen('/home_tmp/sasajima/DATA/green.dat','r');
mgreen=textscan(Fid5,'%f %f');
fclose(Fid5);
mgreen=cell2mat(mgreen);
gg=size(mgreen,1);

Fid7=fopen('/home/sasajima/Dropbox/yellow/PAC_red.txt','r');
PACred=textscan(Fid7,'%f %f');
fclose(Fid7);
PACred=cell2mat(PACred);

Fid8=fopen('/home/sasajima/Dropbox/yellow/PAC_black.txt','r');
PACblack=textscan(Fid8,'%f %f');
fclose(Fid8);
PACblack=cell2mat(PACblack);

Fid9=fopen('/home_tmp/sasajima/DATA/red.dat','r');
mred=textscan(Fid9,'%f %f');
fclose(Fid9);
mred=cell2mat(mred);
rr=size(mred,1);

Fid10=fopen('/home_tmp/sasajima/DATA/black.dat','r');
mblack=textscan(Fid10,'%f %f');
fclose(Fid10);
mblack=cell2mat(mblack);
blbl=size(mblack,1);

Fid11=fopen('/home_tmp/sasajima/DATA/m_all.dat','r');
mall=textscan(Fid11,'%f %f');
fclose(Fid11);
mall=cell2mat(mall);
alal=size(mall,1);

s=0;
for r=1:rr
    Red=inpolygon(mred(r,1),mred(r,2),PACred(:,1),PACred(:,2));
    if Red==0
    else
        s=s+1;
        slon(s,1)=mred(r,1);
        slat(s,1)=mred(r,2);
    end
end

sred=s;

for b=1:bb
    Blue=inpolygon(mblue(b,1),mblue(b,2),PACblue(:,1),PACblue(:,2));
    BRed=inpolygon(mblue(b,1),mblue(b,2),PACred(:,1),PACred(:,2));
    if (Blue==1)&&(BRed==0)
        s=s+1;
        slon(s,1)=mblue(b,1);
        slat(s,1)=mblue(b,2);
    else
    end
end

sblue=s-sred;

for g=1:gg
    Green=inpolygon(mgreen(g,1),mgreen(g,2),PACgreen(:,1),PACgreen(:,2));
    GBlue=inpolygon(mgreen(g,1),mgreen(g,2),PACblue(:,1),PACblue(:,2));
    if (Green==1)&&(GBlue==0)
        s=s+1;
        slon(s,1)=mgreen(g,1);
        slat(s,1)=mgreen(g,2);
    else
    end
end

sgreen=s-sred-sblue;

for bl=1:blbl
    Black=inpolygon(mblack(bl,1),mblack(bl,2),PACblack(:,1),PACblack(:,2));
    BlGreen=inpolygon(mblack(bl,1),mnlack(bl,2),PACgreen(:,1),PACgreen(:,2));
    if (Black==1)&&(BlGreen==0)
        s=s+1;
        slon(s,1)=mblack(bl,1);
        slat(s,1)=mblack(bl,2);
    else
    end
end

sblack=s-sred-sblue-sgreen;


for al=1:alal
    All=inpolygon(mall(al,1),mall(al,2),bound(:,1),bound(:,2));
    AlBlack=inpolygon(mall(al,1),mall(al,2),PACblack(:,1),PACblack(:,2));
    if (All==1)&&(AlBlack==0)
        s=s+1;
        slon(s,1)=mall(al,1);
        slat(s,1)=mall(al,2);
    else
    end
end

soutblack=s-sred-sblue-sgreen-sblack;
sall=s;
pause

sll(:,1:2)=[slon(:,1),slat(:,1)];

%====================================================
tri = delaunay(slon,slat);
%====================================================
ntri=length(tri);
E=scatteredInterpolant(dep_main(:,1),dep_main(:,2),dep_main(:,3),'natural');%depth of interplate -km
F=scatteredInterpolant(dep_sub(:,1),dep_sub(:,2),dep_sub(:,3),'natural');%depth of seafloor -km

tri3=zeros(ntri,3,3);
triC=zeros(ntri,3);
nn=0;
for n=1:ntri
  lon1=slon(tri(n,1));
  lat1=slat(tri(n,1));
  dep1=-F(lon1,lat1)-E(lon1,lat1);%+km
  lon2=slon(tri(n,2));
  lat2=slat(tri(n,2));
  dep2=-F(lon2,lat2)-E(lon2,lat2);
  lon3=slon(tri(n,3));
  lat3=slat(tri(n,3));
  dep3=-F(lon3,lat3)-E(lon3,lat3);
  
  glon=(lon1+lon2+lon3)/3;
  glat=(lat1+lat2+lat3)/3;
  gdep=(dep1+dep2+dep3)/3;
  
  IN=inpolygon(glon,glat,bound(:,1),bound(:,2));
  
  if gdep==0
  elseif IN==0
  else
    nn=nn+1;
    triC(nn,1:3)=[glon,glat,gdep];
    tri3(nn,1,1)=lon1;
    tri3(nn,1,2)=lon2;
    tri3(nn,1,3)=lon3;
    tri3(nn,2,1)=lat1;
    tri3(nn,2,2)=lat2;
    tri3(nn,2,3)=lat3;
    tri3(nn,3,1)=dep1;
    tri3(nn,3,2)=dep2;
    tri3(nn,3,3)=dep3;
  end
end
tri3(nn+1:ntri,:)=[];
triC(nn+1:ntri,:)=[];

Fid6=fopen('/home_tmp/sasajima/DATA/PAC_tri3.txt','w');
for w=1:nn
    for u=1:3
        fprintf(Fid6,'%10.4f %9.4f %10.4f\n',tri3(w,1,u),tri3(w,2,u),tri3(w,3,u));
    end
end
fclose(Fid6);

end
%% makexyz.m
function  [xyz]=makexyz
%by Ryohei Sasajima
%last 2014/02/02
%====================================================
%{
Fid1=fopen('/home_tmp/sasajima/DATA/GEODETIC_DATA/coordinates_F3/sitelocate.txt','r');
geonet=textscan(Fid1,'%f %f %f');
fclose(Fid1);
geonet=cell2mat(geonet);
%}

Fid2=fopen('/home_tmp/sasajima/DATA/xyz01.dat','r');
xyz01=textscan(Fid2,'%f %f');
fclose(Fid2);
xyz=cell2mat(xyz01);

end
%====================================================
%% Slip2xyll.m
function [xySlip]=Slip2xyll(Slip,sitaS)
xySlip(:,1)=Slip(:,1).*cos(sitaS(:))+Slip(:,2).*sin(sitaS(:));
xySlip(:,2)=-Slip(:,1).*sin(sitaS(:))+Slip(:,2).*cos(sitaS(:));
end
%====================================================

%% strike_dip.m
function[sitaS,sitaD,normVec]=strike_dip(trixyzC,trixyz3)
%Define strike and dip of fault%
n=size(trixyzC,1);
ca=zeros(n,3);
cb=zeros(n,3);
InormVec=zeros(n,3);
normVec=zeros(n,3);
strikeVec=zeros(n,3);
dipVec=zeros(n,3);
sitaS=zeros(n,1);
sitaD=zeros(n,1);

for i=1:n
    ca(i,:)=[trixyz3(i,1)-trixyz3(i,7),trixyz3(i,2)-trixyz3(i,8),trixyz3(i,3)-trixyz3(i,9)];
    cb(i,:)=[trixyz3(i,1)-trixyz3(i,4),trixyz3(i,2)-trixyz3(i,5),trixyz3(i,3)-trixyz3(i,6)];
    InormVec(i,:)              = cross(ca(i,:),cb(i,:),2);% direct to deeper %
    normVec(i,:)              = InormVec(i,:)./(sqrt((InormVec(i,1)).^2+(InormVec(i,2)).^2+(InormVec(i,3).^2)));
    if (normVec(i,3) < 0) % Enforce clockwise circulation
        normVec(i,:)               = -normVec(i,:);
    end
    strikeVec(i,:)              = [-sin(atan2(normVec(i,2),normVec(i,1))),cos(atan2(normVec(i,2),normVec(i,1))),0];%direct to left hand of who toward dip direction
    dipVec(i,:)                 = cross(normVec(i,:), strikeVec(i,:),2);
    %if normVec(1)==0 and normVec(2)==0
    % strikeVec=[1 0 0];
    % dipVec   =[0 1 0];
    % normVec  =[0 0 1];
    %end
    sitaS(i)=atan2(strikeVec(i,2),strikeVec(i,1)); %from x-axis to y-axis rotation [rad]%
    sitaD(i)=asin(dipVec(i,3));
end
end

%% PLTXY
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
