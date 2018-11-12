%% PassiveSlip.m
function PassiveSlipKim
%by Ryohei Sasajima
%final 2013/12/23


%[triC,tri3,tri,sll]=make_test_trill;
%save('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/tri','triC','tri3','tri','sll');
load('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/tri','triC','tri3','tri','sll');


%[trixyzC,trixyz3,sxyz]=trill2trixyz(triC,tri3,sll);
%save('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/trixyz','trixyzC','trixyz3','sxyz');
load('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/trixyz','trixyzC','trixyz3','sxyz');

%Gtrixyz3=trixyz3;


%[sitaS,sitaD,normVec]=strike_dip(trixyzC,trixyz3);
%save('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/sita','sitaS','sitaD','normVec');
load('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/sita','sitaS','sitaD','normVec');

%trixyz3=Gtrixyz3;

n=size(sitaS,1)

%figure(31); quiver(trixyzC(:,1),-trixyzC(:,2),sitaS,sitaD,1,'b')

[si]=makeG_s_Q(trixyzC,trixyz3,sitaS,sitaD,normVec);

[di]=makeG_d_Q(trixyzC,trixyz3,sitaS,sitaD,normVec);

%{

 [xyz]=makexyz;
 save('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/xyz','xyz');
 %load('/home_tmp/sasajima/DATA/PassiveSlip/PAC_test/xyz','xyz');

 [si]=makeG_s_O(xyz,trixyz3);

 [di]=makeG_d_O(xyz,trixyz3);

 [sUxyz,dUxyz]=loadMAT2(xyz);

 [sSsn,sSdn,dSsn,dSdn]=loadMAT(trixyzC);

 finish_loadMAT=0
 
 zeroV=zeros(n,1);
 %tri3(1243,3,:)
 whos
 %figure(22); quiver(trixyzC(:,1),-trixyzC(:,2),dSdn(:,1155),zeroV(:,1),5,'b')
 
 %
 %triC(1155,:)
 %
 [Slip]=defineSlipQ(triC,sitaS);
 save('/home_tmp/sasajima/DATA/MAT/sdSlip_Q1','Slip');
 %load('/home_tmp/sasajima/DATA/MAT/sdSlip_Q1','Slip');
 
 [xySlip]=Slip2xyll(Slip,sitaS);
 save('/home_tmp/sasajima/DATA/MAT/xySlip_Q1','xySlip');
 %load('/home_tmp/sasajima/DATA/MAT/xySlip_Q1','xySlip');
 %figure(103); quiver(triC(:,1),triC(:,2),xySlip(:,1),-xySlip(:,2),'b')
 figure(101); quiver(trixyzC(:,1),-trixyzC(:,2),xySlip(:,1),-xySlip(:,2),'g')
 
 ll(:,1:2)=triC(:,1:2);
 %[vnorm,v,NARFA]=SasaEular2Velo(ll);

 []=relavec2vecVxy(NARFA)

 [xyBslip]=BackSlipxy(mcmcC,triC,trixyzC,sitaS,accuV,vecVxy)

 [A_Ssn,A_Sdn]=I_CalcStrain(Slip,sSsn,sSdn,dSsn,dSdn);
 save('/home_tmp/sasajima/DATA/MAT/A_Ssn_2','A_Ssn','A_Sdn');
 %load('/home_tmp/sasajima/DATA/MAT/A_Sdn_2','A_Ssn','A_Sdn');
 figure(201); quiver(trixyzC(:,1),-trixyzC(:,2),A_Ssn,A_Sdn,5,'b')
 
 %whos
 
 [F_Slip]=OutofAsperity01(A_Ssn,A_Sdn,Slip,sSsn,sSdn,dSsn,dSdn);
 save('/home_tmp/sasajima/DATA/MAT/FSlip','F_Slip');
 %load('/home_tmp/sasajima/DATA/MAT/FSlip','F_Slip');

 [xyFSlip]=FSlip2xyll(F_Slip,sitaS);
 save('/home_tmp/sasajima/DATA/MAT/xyFSlip_Q1','xyFSlip');
 %load('/home_tmp/sasajima/DATA/MAT/xyFSlip_Q1','xyFSlip');
 
 
 %Fid=fopen('./QxySlip11.dat','w');

 
 %n=size(Slip,1);
 
 %for i=1:n
 %fprintf(Fid, '%11.3f %11.3f %9.5f %9.5f\n',trixyzC(i,1), trixyzC(i,2), xyFSlip(i,1), xyFSlip(i,2));
 %end
 %fclose(Fid);
 
 %
 figure(1); quiver(trixyzC(:,1),-trixyzC(:,2),xyFSlip(:,1),-xyFSlip(:,2),2,'g')
 %hold on
 %triplot(tri,1000.*sxyz(:,1),1000.*sxyz(:,2));
 
 [F_Ssn,F_Sdn]=I_CalcStrainF(F_Slip,sSsn,sSdn,dSsn,dSdn);
 save('/home_tmp/sasajima/DATA/MAT/FStrain','F_Ssn','F_Sdn');
 %load('/home_tmp/sasajima/DATA/MAT/FStrain','F_Ssn','F_Sdn');
 figure(401); quiver(trixyzC(:,1),-trixyzC(:,2),F_Ssn,F_Sdn,5,'b')
 
 %whos
 
 %absO_Ssn=abs(O_Ssn);
 %absO_Sdn=abs(O_Sdn);
 %sumAOSsn=sum(absO_Ssn,1)
 %sumAOSdn=sum(absO_Sdn,1)
 
 [FO_Ssn,FO_Sdn]=OutofAsperity11(F_Ssn,F_Sdn,Slip);
 save('/home_tmp/sasajima/DATA/MAT/FOStrain','FO_Ssn','FO_Sdn');
 %load('/home_tmp/sasajima/DATA/MAT/FOStrain','FO_Ssn','FO_Sdn');
 absFO_Ssn=abs(FO_Ssn);
 absFO_Sdn=abs(FO_Sdn);
 sumF0Ssn=sum(absFO_Ssn,1)
 sumF0Sdn=sum(absFO_Sdn,1)
%}
end
%% AnormaryVecCalc.m
function [AnormaryVecxy]=AnormaryVecCalc(vnorm,SVxy,xyFSlip)

n=length(vnorm);

xyplateV(:,1)=vnorm(:,1).*sin(SVxy(:,1));%x-component of Velocity [m/y]
xyplateV(:,2)=vnorm(:,1).*cos(SVxy(:,1));%y-component of Velocity [m/y]
xyCreep(:,1:2)=xyFSlip(:,1:2)+xyplateV(:,1:2);%[m/y]

for i=1:n
    XinstSVxy(i,1)=atan2(xyCreep(i,2),xyCreep(i,1));%countour clock wise from X [radian]
    instSVxy(i,1)=2.*pi-XinstSVxy(i,1)+pi.*0.5;%clock wise from Y [radian]
    
    if instSVxy(i,1)>=(2.*pi)
        instSVxy(i,1)=instSVxy(i,1)-2.*pi;
    else
        continue
    end
    
end

AnormaryVecxy(:,1)=(SVxy(:,1)-instSVxy(:,1)).*180./pi;%[digree]

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

%% DropValue.m
function DropValue
load('/home_tmp/sasajima/DATA/MAT/tri_23','triC');
load('/home_tmp/sasajima/DATA/MAT/sita_23','sitaD');
%load('/home_tmp/sasajima/DATA/MAT/sRSxz_1','sR2Sxz');
%{
 load('/home/sasajima/DATA/MAT/sRSyz_1','sR2Syz');
 load('/home/sasajima/DATA/MAT/dRSxz_1','dR2Sxz');
 load('/home/sasajima/DATA/MAT/dRSyz_1','dR2Syz');
 load('/home/sasajima/DATA/MAT/Slip_1','Slip');
 load('/home/sasajima/DATA/MAT/define_1','define');
 load('/home/sasajima/DATA/MAT/I1_R2Sxz_1','I1_R2Sxz');
 load('/home/sasajima/DATA/MAT/I1_R2Syz_1','I1_R2Syz');
 load('/home/sasajima/DATA/MAT/I1_Slip_1','I1_Slip');
 load('/home/sasajima/DATA/MAT/I2_R2Sxz_1','I2_R2Sxz');
 load('/home/sasajima/DATA/MAT/I2_R2Syz_1','I2_R2Syz');
 load('/home/sasajima/DATA/MAT/I2_Slip_1','I2_Slip');
%}
sx=triC(:,1);
sy=triC(:,2);
%Slip(:,1)
%Slip(:,2)
%I1_R2Syz
%I1_R2Syz
%dR2Syz
%sR2Syz
%dR2Sxz
%sR2Sxz
%Slip(:,2)
%I1_Slip(:,2)
%I2_Slip(:,2)
%sumSlip=Slip+I1_Slip+I2_Slip;
whos
GdR2Syz=sitaD(:,1)

min_s=0;
max_s=pi./4;
index_max=GdR2Syz>max_s;
index_min=GdR2Syz<min_s;
index_nan=isnan(GdR2Syz);
color_index=jet(128);
RRR=fix((GdR2Syz-min_s)./((max_s-min_s)./128)+1);
RRR(index_max)=128;
RRR(index_min)=1;
RRR(index_nan)=1;
RR(:,1)=color_index(RRR,1);
RR(:,2)=color_index(RRR,2);
RR(:,3)=color_index(RRR,3);
scatter(sx,sy,20,RR,'filled')
%scatter(sx,sy,RSyz,'filled')
%fclose(Fid2);
%figure(15); quiver3(sx,sy,sz,Usumx,Usumy,Usumz,'r')
%figure(16); quiver(sx,sy,Slip(:,1),Slip(:,2),'r'); grid on;
end
%% DropValueQ.m
function DropValueQ
load('/home_tmp/sasajima/DATA/PassiveSlip/Tohoku/tri_Tohoku','triC','tri3','tri','sll');
load('/home_tmp/sasajima/DATA/PassiveSlip/Tohoku/trixyz_Tohoku','trixyzC','trixyz3','sxyz');
load('/home_tmp/sasajima/DATA/PassiveSlip/Tohoku/sita_Tohoku','sitaS','sitaD','normVec');

%[sSsn,sSdn,dSsn,dSdn]=loatMAT(trixyzC);

%load('/home_tmp/sasajima/DATA/MAT/sSsn_Boso','sSsn');
%load('/home_tmp/sasajima/DATA/MAT/sSdn_Q07','sSdn');
%load('/home_tmp/sasajima/DATA/MAT/dSsn_Q07','dSsn');
%load('/home_tmp/sasajima/DATA/MAT/dSdn_Boso','dSdn');

%load('/home_tmp/sasajima/DATA/MAT/sxyz_Q07','sUxyz');
%load('/home_tmp/sasajima/DATA/MAT/dUxyz_Q07','dUxyz');

sx=triC(:,1);
sy=triC(:,2);
sz=triC(:,3);
%Slip(:,1)
%Slip(:,2)
%I1_R2Syz
%I1_R2Syz
%dR2Syz
%sR2Syz
%dR2Sxz
%sR2Sxz
%Slip(:,2)
%I1_Slip(:,2)
%I2_Slip(:,2)
%sumSlip=Slip+I1_Slip+I2_Slip;
whos
%A_Ssn(:);
%dSdn(101,:)
%Value=triC(:,3);
Value=sitaS(:,1);

min_s=pi./2
max_s=pi
index_max=Value>max_s;
index_min=Value<min_s;
index_nan=isnan(Value);
color_index=jet(128);
RRR=fix((Value-min_s)./((max_s-min_s)./128)+1);
RRR(index_max)=128;
RRR(index_min)=1;
RRR(index_nan)=1;
RR(:,1)=color_index(RRR,1);
RR(:,2)=color_index(RRR,2);
RR(:,3)=color_index(RRR,3);
figure(1);
scatter(sx,sy,20,RR,'filled')
%hold on
%figure(2)
%triplot(tri,1000.*sll(:,1),1000.*sll(:,2));
%fclose(Fid2);
%figure(15); quiver3(sx,sy,sz,dUxyz(:,1,1),dUxyz(:,1,2),dUxyz(:,1,3),'r')
%figure(16); quiver(sx,sy,normVec(:,1),-normVec(:,2),'r');
%grid on;
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
k=sum(index_out_asp)
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
%% SV_ll2xy_2boso.m
function [SVxy]=SV_ll2xy_2boso(Vec,NARFA)
%-------------------
%  PLTXY TRANSFORMS (ALAT,ALONG) TO (X,Y)
%  WHEN ICORD.NE.0  PLTXY MAKES NO CHANGE IN
%  TRANSFORMATION BETWEEN (X,Y) AND (ALAT,ALONG).
%-------------------
%constant%
A=6.378160e3;
E2=6.6944541e-3;
E12=6.7395719e-3;
D=5.72958e1;
RD=1.0/D;
R=6371;
%Change Value%
ALON0=142.5;
ALAT0=38.3;

whos

PlateV(:,1:2)=Vec(:,1:2);%velocity of oceanic plate clock wise from North (degree)%

n=length(Vec);

dLAT=0.1./180.*pi;

for j=1:n
    
    %calculation%
    RPV(j,1)=PlateV(j,3)/180*pi;
    
    if cos(RPV(j,1))>=0 %NS component of plateV is N;
        
        AftLAT=Vec(j,2)+dLAT.*180./pi;%digree
        
        if sin(RPV(j,1))>=0 %EW component of plateV is E;
            
            radVec(j,1)=Vec(j,3)./180.*pi;%[radian]
            dPSY=sin(dLAT).*cos(radVec(j,1));%[radian]
            dDISTANCE=2.*R.*sin(dPSY./2);%linear distance[km]
            appaR=R.*sin(pi./2-AftLAT);%[km]
            dLON=asin(dDISTANCE./appaR);%[radian]
            
        else %EW component of plateV is W;
            radVec(j,1)=(360-Vec(j,3))./180.*pi;%radian
            dPSY=sin(dLAT).*cos(radVec(j,1));%radian
            dDISTANCE=2.*R.*sin(dPSY./2);%linear distance[km]
            appaR=R.*sin(pi./2-AftLAT);%[km]
            dLON=asin(dDISTANCE./appaR).*(-1);%[radian]
            
        end
        
        AftLON=Vec(j,1)+dLON.*180./pi;%[digree]
        
    else %NS component of plateV is S;
        
        AftLAT=Vec(j,2)-dLAT.*180./pi;%digree
        
        if sin(PRV(j,1))>=0 %EW component of plateV is E;
            
            radVec(j,1)=(180-Vec(j,3))./180.*pi;%[radian]
            dPSY=sin(dLAT).*cos(radVec(j,1));%[radian]
            dDISTANCE=2.*R.*sin(dPSY./2);%linear distance[km]
            appaR=R.*sin(pi./2-AftLAT);%[km]
            dLON=asin(dDISTANCE./appaR);%[radian]
            
        else %EW component of plateV is W;
            radVec(j,1)=(Vec(j,3)-180)./180.*pi;%radian
            dPSY=sin(dLAT).*cos(radVec(j,1));%radian
            dDISTANCE=2.*R.*sin(dPSY./2);%linear distance[km]
            appaR=R.*sin(pi./2-AftLAT);%[km]
            dLON=asin(dDISTANCE./appaR).*(-1);%[radian]
            
        end
        
        AftLON=Vec(j,1)+dLON.*180./pi;%digree
        
    end
    
    %lon-lat to xyz-sphere-coordinate to xy-coordinate
    Vlon=[PlateV(j,1);AftLON];
    Vlat=[PlateV(j,2);AftLAT];
    
    for i=1:2;
        ALON=Vlon(i);
        ALAT=Vlat(i);
        
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
        
        VX(i)=X;
        VY(i)=Y;
    end
    SVxy(j,1)=atan2(VY(2)-VY(1),VX(2)-VX(1));%counter clock wise from Y-axis (radian)%
    SVxyzdig=SVxy.*180./pi;
    
end
end


%% SV_ll2xy_boso.m
function [SVxy]=SV_ll2xy_boso(Vec)
%-------------------
%  PLTXY TRANSFORMS (ALAT,ALONG) TO (X,Y)
%  WHEN ICORD.NE.0  PLTXY MAKES NO CHANGE IN
%  TRANSFORMATION BETWEEN (X,Y) AND (ALAT,ALONG).
%-------------------
%constant%
A=6.378160e3;
E2=6.6944541e-3;
E12=6.7395719e-3;
D=5.72958e1;
RD=1.0/D;
R=6371;
%Change Value%
ALON0=142.5;
ALAT0=38.3;

whos

PlateV(:,1:3)=Vec(:,1:3);%velocity of oceanic plate clock wise from North (lon,lat,degree)%

n=length(Vec);

dLAT=0.1./180.*pi;

for j=1:n
    
    %calculation%
    RPV(j,1)=PlateV(j,3)/180*pi;
    
    if cos(RPV(j,1))>=0 %NS component of plateV is N;
        
        AftLAT=Vec(j,2)+dLAT.*180./pi;%digree
        
        if sin(RPV(j,1))>=0 %EW component of plateV is E;
            
            radVec(j,1)=Vec(j,3)./180.*pi;%[radian]
            dPSY=sin(dLAT).*cos(radVec(j,1));%[radian]
            dDISTANCE=2.*R.*sin(dPSY./2);%linear distance[km]
            appaR=R.*sin(pi./2-AftLAT);%[km]
            dLON=asin(dDISTANCE./appaR);%[radian]
            
        else %EW component of plateV is W;
            radVec(j,1)=(360-Vec(j,3))./180.*pi;%radian
            dPSY=sin(dLAT).*cos(radVec(j,1));%radian
            dDISTANCE=2.*R.*sin(dPSY./2);%linear distance[km]
            appaR=R.*sin(pi./2-AftLAT);%[km]
            dLON=asin(dDISTANCE./appaR).*(-1);%[radian]
            
        end
        
        AftLON=Vec(j,1)+dLON.*180./pi;%[digree]
        
    else %NS component of plateV is S;
        
        AftLAT=Vec(j,2)-dLAT.*180./pi;%digree
        
        if sin(PRV(j,1))>=0 %EW component of plateV is E;
            
            radVec(j,1)=(180-Vec(j,3))./180.*pi;%[radian]
            dPSY=sin(dLAT).*cos(radVec(j,1));%[radian]
            dDISTANCE=2.*R.*sin(dPSY./2);%linear distance[km]
            appaR=R.*sin(pi./2-AftLAT);%[km]
            dLON=asin(dDISTANCE./appaR);%[radian]
            
        else %EW component of plateV is W;
            radVec(j,1)=(Vec(j,3)-180)./180.*pi;%radian
            dPSY=sin(dLAT).*cos(radVec(j,1));%radian
            dDISTANCE=2.*R.*sin(dPSY./2);%linear distance[km]
            appaR=R.*sin(pi./2-AftLAT);%[km]
            dLON=asin(dDISTANCE./appaR).*(-1);%[radian]
            
        end
        
        AftLON=Vec(j,1)+dLON.*180./pi;%digree
        
    end
    
    %lon-lat to xyz-sphere-coordinate to xy-coordinate
    Vlon=[PlateV(j,1);AftLON];
    Vlat=[PlateV(j,2);AftLAT];
    
    for i=1:2;
        ALON=Vlon(i);
        ALAT=Vlat(i);
        
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
        
        VX(i)=X;
        VY(i)=Y;
    end
    SVxy(j,1)=atan2(VY(2)-VY(1),VX(2)-VX(1));%counter clock wise from Y-axis (radian)%
    SVxyzdig=SVxy.*180./pi;
    
end
end

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

for i=1:n;
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
    
    if dltdiglon>0;
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

%% SasaTriDislocaBoso.m
function SasaTriDislocaBoso

%[triC,tri3,tri,sll]=Sasa_make_trill;
%save('/home_tmp/sasajima/DATA/MAT/tri_Boso','triC','tri3','tri','sll');
load('/home_tmp/sasajima/DATA/MAT/tri_Boso','triC','tri3','tri','sll');


%[trixyzC,trixyz3,sxyz]=Sasa_ll2xyz(triC,tri3,sll);
%save('/home_tmp/sasajima/DATA/MAT/trixyz_Boso','trixyzC','trixyz3','sxyz');
load('/home_tmp/sasajima/DATA/MAT/trixyz_Boso','trixyzC','trixyz3','sxyz');


%[sitaS,sitaD,normVec]=strike_dip(trixyzC,trixyz3);
%save('/home_tmp/sasajima/DATA/MAT/sitaBoso','sitaS','sitaD','normVec');
load('/home_tmp/sasajima/DATA/MAT/sitaBoso','sitaS','sitaD','normVec');

n=size(sitaS,1)

%figure(31); quiver(trixyzC(:,1),-trixyzC(:,2),sitaS,sitaD,1,'b')
%
%[sSsn,sSdn]=makeG_s_Q(trixyzC,trixyz3,sitaS,sitaD,normVec);
%save('/home_tmp/sasajima/DATA/MAT/sSsn_Boso','sSsn','-v7.3');
%save('/home_tmp/sasajima/DATA/MAT/sSdn_Boso','sSdn','-v7.3');
%save('/home_tmp/sasajima/DATA/MAT/sxyz_Q11','sUxyz');
%[dSsn,dSdn]=makeG_d_Q(trixyzC,trixyz3,sitaS,sitaD,normVec);
%save('/home_tmp/sasajima/DATA/MAT/dSsn_Boso','dSsn','-v7,3');
%save('/home_tmp/sasajima/DATA/MAT/dSdn_Boso','dSdn','-v7.3');
%save('/home_tmp/sasajima/DATA/MAT/dUxyz_Q11','dUxyz');

%{
 load('/home_tmp/sasajima/DATA/MAT/sSsn_Boso','sSsn');
 load('/home_tmp/sasajima/DATA/MAT/sSdn_Boso','sSdn');
 load('/home_tmp/sasajima/DATA/MAT/dSsn_Boso','dSsn');
 load('/home_tmp/sasajima/DATA/MAT/dSdn_Boso','dSdn');
%}

%
[sSsn,sSdn,dSsn,dSdn]=loatMAT(trixyzC);

finish_loatMAT=0
zeroV=zeros(n,1);
%tri3(1243,3,:)
whos
%figure(22); quiver(trixyzC(:,1),-trixyzC(:,2),dSdn(:,1155),zeroV(:,1),5,'b')

%
%triC(1155,:)
%
[llF]=defineEllipsll;

Velo(1,1)=0.0290;%(m/y)
Velo(2,1)=0.0267;
Velo(3,1)=0.0279;
Velo(4,1)=0.0263;
Velo(5,1)=0.0265;
Velo(6,1)=0.0255;
Velo(7,1)=0;

%{
 Vec(1,1:3)=[139.350,35.300,325.0];%epilon,epilat,PlateVec crockwise from North [digree]
 Vec(2,1:3)=[139.680,35.140,336.5];
 Vec(3,1:3)=[139.950,35.660,329.4];
 Vec(4,1:3)=[140.160,34.680,333.6];
 Vec(5,1:3)=[139.970,34.920,335.0];
 Vec(6,1:3)=[140.595,34.700,332.5];
 Vec(7,1:3)=[140.500,35.075,346.0];
%}

%triC(1:10,:)
%llF(1,:,:)
%
%[SVxy]=SV_ll2xy_boso(Vec)
SVxy(1,1)=215.0./180.*pi;
SVxy(2,1)=215.0./180.*pi;
SVxy(3,1)=215.0./180.*pi;
SVxy(4,1)=215.0./180.*pi;
SVxy(5,1)=215.0./180.*pi;
SVxy(6,1)=215.0./180.*pi;
SVxy(7,1)=215.0./180.*pi;

%
[Slip]=defineSlipBoso(triC,sitaS,llF,Velo,SVxy);
%finish_defineSlip=0
%
save('/home_tmp/sasajima/DATA/MAT/bososlip/sdSlip_boso','Slip');
%load('/home_tmp/sasajima/DATA/MAT/bososlip/sdSlip_boso','Slip');
%Slip
figure(1); quiver(trixyzC(:,1),-trixyzC(:,2),Slip(:,1),Slip(:,2),5,'b')
%
whos

%
[A_Ssn,A_Sdn]=I_CalcStrain(Slip,sSsn,sSdn,dSsn,dSdn);
finish_1stStrain=0
%save('/home_tmp/sasajima/DATA/MAT/bososlip/A_Ssn_boso','A_Ssn','A_Sdn');
%load('/home_tmp/sasajima/DATA/MAT/bososlip/A_Sdn_boso','A_Ssn','A_Sdn');
figure(201); quiver(trixyzC(:,1),-trixyzC(:,2),A_Ssn,A_Sdn,5,'b')

%whos

[F_Slip]=OutofAsperity01(A_Ssn,A_Sdn,Slip,sSsn,sSdn,dSsn,dSdn);
finish_F_Slip=0
%save('/home_tmp/sasajima/DATA/MAT/bososlip/FSlip_boso','F_Slip');
%load('/home_tmp/sasajima/DATA/MAT/bososlip/FSlip_boso','F_Slip');

[xyFSlip]=FSlip2xyll(F_Slip,sitaS);
finish_xyFSlip=0
save('/home_tmp/sasajima/DATA/MAT/bososlip/xyFSlip_boso','xyFSlip');
%
load('/home_tmp/sasajima/DATA/MAT/bososlip/xyFSlip_boso','xyFSlip');
figure(4); quiver(trixyzC(:,1),-trixyzC(:,2),xyFSlip(:,1),xyFSlip(:,2),5,'b')

%Fid=fopen('./QxySlip11.dat','w');

%n=size(Slip,1);
%

%{
[FSlipll]=vecxy2ll(trixyzC,xyFSlip)
 finish_FSlipll=0
 save('/home_tmp/sasajima/DATA/MAT/bososlip/FSlipll_boso','FSlipll');
 %load('/home_tmp/sasajima/DATA/MAT/bososlip/FSlipll_boso','FSlipll');
 figure(2); quiver(triC(:,1),triC(:,2),FSlipll(:,1),FSlipll(:,2),5,'b')
%}

%

%
ll=triC;

[vnorm,v,NARFA]=SasaEular2Velo(ll)

Vec=NARFA;
[SVxy]=SV_ll2xy_boso(Vec)

%{
Value=SVxy(:,1).*180./pi;
sx=trixyzC(:,1);
sy=trixyzC(:,2);

min_s=35
max_s=36
index_max=Value>max_s;
index_min=Value<min_s;
index_nan=isnan(Value);
color_index=jet(128);
RRR=fix((Value-min_s)./((max_s-min_s)./128)+1);
RRR(index_max)=128;
RRR(index_min)=1;
RRR(index_nan)=1;
RR(:,1)=color_index(RRR,1);
RR(:,2)=color_index(RRR,2);
RR(:,3)=color_index(RRR,3);
figure(19);
scatter(sx,-sy,20,RR,'filled')
%}

%
[AnormaryVec]=AnormaryVecCalc(vnorm,SVxy,xyFSlip)

Fid2=fopen('/home_tmp/sasajima/DATA/AnormaryVec2.txt','w');
m=length(ll);
for i=1:m
    fprintf(Fid2, '%12.5f %12.5f %12.5f\n',ll(i,1),ll(i,2),AnormaryVec(i,1));
end
fclose(Fid2);
%

%for i=1:n
%fprintf(Fid, '%11.3f %11.3f %9.5f %9.5f\n',trixyzC(i,1), trixyzC(i,2), xyFSlip(i,1), xyFSlip(i,2));
%end
%fclose(Fid);

%{
 figure(1); quiver(trixyzC(:,1),-trixyzC(:,2),xyFSlip(:,1),-xyFSlip(:,2),2,'g')
 %hold on
 %triplot(tri,1000.*sxyz(:,1),1000.*sxyz(:,2));
 
 [F_Ssn,F_Sdn]=I_CalcStrainF(F_Slip,sSsn,sSdn,dSsn,dSdn);
 %save('/home_tmp/sasajima/DATA/MAT/FStrain','F_Ssn','F_Sdn');
 %load('/home_tmp/sasajima/DATA/MAT/FStrain','F_Ssn','F_Sdn');
 figure(401); quiver(trixyzC(:,1),-trixyzC(:,2),F_Ssn,F_Sdn,5,'b')
 
 %whos
 
 %absO_Ssn=abs(O_Ssn);
 %absO_Sdn=abs(O_Sdn);
 %sumAOSsn=sum(absO_Ssn,1)
 %sumAOSdn=sum(absO_Sdn,1)
 
 [FO_Ssn,FO_Sdn]=OutofAsperity11(F_Ssn,F_Sdn,Slip);
 %save('/home_tmp/sasajima/DATA/MAT/FOStrain','FO_Ssn','FO_Sdn');
 %load('/home_tmp/sasajima/DATA/MAT/FOStrain','FO_Ssn','FO_Sdn');
 absFO_Ssn=abs(FO_Ssn);
 absFO_Sdn=abs(FO_Sdn);
 sumF0Ssn=sum(absFO_Ssn,1)
 sumF0Sdn=sum(absFO_Sdn,1)
%}
end
%% SasaTriDislocaQ.m
function SasaTriDislocaQ

%[triC,tri3,tri,sll]=Sasa_make_trill;
%save('/home_tmp/sasajima/DATA/MAT/tri_Qalld1','triC','tri3','tri','sll');
load('/home_tmp/sasajima/DATA/MAT/tri_Qalld1','triC','tri3','tri','sll');


%[trixyzC,trixyz3,sxyz]=Sasa_ll2xyz(triC,tri3,sll);
%save('/home_tmp/sasajima/DATA/MAT/trixyz_Qalld1','trixyzC','trixyz3','sxyz');
load('/home_tmp/sasajima/DATA/MAT/trixyz_Qalld1','trixyzC','trixyz3','sxyz');


%[sitaS,sitaD,normVec]=strike_dip(trixyzC,trixyz3);
%save('/home_tmp/sasajima/DATA/MAT/sitaQalld1','sitaS','sitaD','normVec');
load('/home_tmp/sasajima/DATA/MAT/sitaQalld1','sitaS','sitaD','normVec');

n=size(sitaS,1)

%figure(31); quiver(trixyzC(:,1),-trixyzC(:,2),sitaS,sitaD,1,'b')
%
[sSsn,sSdn]=makeG_s_Q(trixyzC,trixyz3,sitaS,sitaD,normVec);
save('/home_tmp/sasajima/DATA/MAT/sSsn_Qalld1','sSsn');
save('/home_tmp/sasajima/DATA/MAT/sSdn_Qalld1','sSdn');
%save('/home_tmp/sasajima/DATA/MAT/sxyz_Q11','sUxyz');
[dSsn,dSdn]=makeG_d_Q(trixyzC,trixyz3,sitaS,sitaD,normVec);
save('/home_tmp/sasajima/DATA/MAT/dSsn_Qalld1','dSsn');
save('/home_tmp/sasajima/DATA/MAT/dSdn_Qalld1','dSdn');
%save('/home_tmp/sasajima/DATA/MAT/dUxyz_Q11','dUxyz');

%{
 load('/home_tmp/sasajima/DATA/MAT/sSsn_Qalld1','sSsn');
 load('/home_tmp/sasajima/DATA/MAT/sSdn_Qalld1','sSdn');
 load('/home_tmp/sasajima/DATA/MAT/dSsn_Qalld1','dSsn');
 load('/home_tmp/sasajima/DATA/MAT/dSdn_Qalld1','dSdn');
 
 
 zeroV=zeros(n,1);
 %tri3(1243,3,:)
 whos
 %figure(22); quiver(trixyzC(:,1),-trixyzC(:,2),dSdn(:,1155),zeroV(:,1),5,'b')
 
 %
 %triC(1155,:)
 %
 [Slip]=defineSlipQ(triC,sitaS);
 save('/home_tmp/sasajima/DATA/MAT/sdSlip_Q1','Slip');
 %load('/home_tmp/sasajima/DATA/MAT/sdSlip_Q1','Slip');
 
 [xySlip]=Slip2xyll(Slip,sitaS);
 save('/home_tmp/sasajima/DATA/MAT/xySlip_Q1','xySlip');
 %load('/home_tmp/sasajima/DATA/MAT/xySlip_Q1','xySlip');
 %figure(103); quiver(triC(:,1),triC(:,2),xySlip(:,1),-xySlip(:,2),'b')
 figure(101); quiver(trixyzC(:,1),-trixyzC(:,2),xySlip(:,1),-xySlip(:,2),'g')
 
 [A_Ssn,A_Sdn]=I_CalcStrain(Slip,sSsn,sSdn,dSsn,dSdn);
 save('/home_tmp/sasajima/DATA/MAT/A_Ssn_2','A_Ssn','A_Sdn');
 %load('/home_tmp/sasajima/DATA/MAT/A_Sdn_2','A_Ssn','A_Sdn');
 figure(201); quiver(trixyzC(:,1),-trixyzC(:,2),A_Ssn,A_Sdn,5,'b')
 
 %whos
 
 [F_Slip]=OutofAsperity01(A_Ssn,A_Sdn,Slip,sSsn,sSdn,dSsn,dSdn);
 save('/home_tmp/sasajima/DATA/MAT/FSlip','F_Slip');
 %load('/home_tmp/sasajima/DATA/MAT/FSlip','F_Slip');

 [xyFSlip]=FSlip2xyll(F_Slip,sitaS);
 save('/home_tmp/sasajima/DATA/MAT/xyFSlip_Q1','xyFSlip');
 %load('/home_tmp/sasajima/DATA/MAT/xyFSlip_Q1','xyFSlip');
 
 
 %Fid=fopen('./QxySlip11.dat','w');

 
 %n=size(Slip,1);
 
 %for i=1:n
 %fprintf(Fid, '%11.3f %11.3f %9.5f %9.5f\n',trixyzC(i,1), trixyzC(i,2), xyFSlip(i,1), xyFSlip(i,2));
 %end
 %fclose(Fid);
 
 %
 figure(1); quiver(trixyzC(:,1),-trixyzC(:,2),xyFSlip(:,1),-xyFSlip(:,2),2,'g')
 %hold on
 %triplot(tri,1000.*sxyz(:,1),1000.*sxyz(:,2));
 
 [F_Ssn,F_Sdn]=I_CalcStrainF(F_Slip,sSsn,sSdn,dSsn,dSdn);
 save('/home_tmp/sasajima/DATA/MAT/FStrain','F_Ssn','F_Sdn');
 %load('/home_tmp/sasajima/DATA/MAT/FStrain','F_Ssn','F_Sdn');
 figure(401); quiver(trixyzC(:,1),-trixyzC(:,2),F_Ssn,F_Sdn,5,'b')
 
 %whos
 
 %absO_Ssn=abs(O_Ssn);
 %absO_Sdn=abs(O_Sdn);
 %sumAOSsn=sum(absO_Ssn,1)
 %sumAOSdn=sum(absO_Sdn,1)
 
 [FO_Ssn,FO_Sdn]=OutofAsperity11(F_Ssn,F_Sdn,Slip);
 save('/home_tmp/sasajima/DATA/MAT/FOStrain','FO_Ssn','FO_Sdn');
 %load('/home_tmp/sasajima/DATA/MAT/FOStrain','FO_Ssn','FO_Sdn');
 absFO_Ssn=abs(FO_Ssn);
 absFO_Sdn=abs(FO_Sdn);
 sumF0Ssn=sum(absFO_Ssn,1)
 sumF0Sdn=sum(absFO_Sdn,1)
%}
end
%% Sasa_ll2xyz.m
function [trixyzC,trixyz3,sxyz]=trill2trixyz(triC,tri3,sll)
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
ALON0=142.5;
ALAT0=38.3;
nn=length(triC);

trixyz3=zeros(nn,3,3);

mm=length(sll);

for i=1:mm
    ALON=sll(i,1);
    ALAT=sll(i,2);
    
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
    sxyz(i,1:2)=[X,Y];
end

for i=1:nn
    ALON=tri3(i,1,1);
    ALAT=tri3(i,2,1);
    
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
    trixyz3(i,1,1)=[X*1000];
    trixyz3(i,2,1)=[-Y*1000];
    trixyz3(i,3,1)=[tri3(i,3,1)*1000];
end

for i=1:nn
    ALON=tri3(i,1,2);
    ALAT=tri3(i,2,2);
    
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
    trixyz3(i,1,2)=[X*1000];
    trixyz3(i,2,2)=[-Y*1000];
    trixyz3(i,3,2)=[tri3(i,3,2)*1000];
end

for i=1:nn
    ALON=tri3(i,1,3);
    ALAT=tri3(i,2,3);
    
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
    trixyz3(i,1,3)=[X*1000];
    trixyz3(i,2,3)=[-Y*1000];
    trixyz3(i,3,3)=[tri3(i,3,3)*1000];
end

trixyzC=zeros(nn,3);

for i=1:nn
    ALON=triC(i,1);
    ALAT=triC(i,2);
    
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
    trixyzC(i,1)=[X*1000];
    trixyzC(i,2)=[-Y*1000];
    trixyzC(i,3)=[triC(i,3)*1000];
end

end
%====================================================

%% Sasa_ll2xyz02.m
function [xy]=Sasa_ll2xyz02(sll)
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
ALON0=140.0;
ALAT0=35.0;
nn=length(triC);

trixyz3=zeros(nn,3,3);

mm=length(sll);

for i=1:mm
    ALON=sll(i,1);
    ALAT=sll(i,2);
    
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
    sxyz(i,1:2)=[X,Y];
end

for i=1:nn
    ALON=tri3(i,1,1);
    ALAT=tri3(i,2,1);
    
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
    trixyz3(i,1,1)=[X*1000];
    trixyz3(i,2,1)=[-Y*1000];
    trixyz3(i,3,1)=[tri3(i,3,1)*1000];
end

for i=1:nn
    ALON=tri3(i,1,2);
    ALAT=tri3(i,2,2);
    
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
    trixyz3(i,1,2)=[X*1000];
    trixyz3(i,2,2)=[-Y*1000];
    trixyz3(i,3,2)=[tri3(i,3,2)*1000];
end

for i=1:nn
    ALON=tri3(i,1,3);
    ALAT=tri3(i,2,3);
    
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
    trixyz3(i,1,3)=[X*1000];
    trixyz3(i,2,3)=[-Y*1000];
    trixyz3(i,3,3)=[tri3(i,3,3)*1000];
end

trixyzC=zeros(nn,3);

for i=1:nn
    ALON=triC(i,1);
    ALAT=triC(i,2);
    
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
    trixyzC(i,1)=[X*1000];
    trixyzC(i,2)=[-Y*1000];
    trixyzC(i,3)=[triC(i,3)*1000];
end

end
%====================================================

%% Sasa_make_trill.m
%====================================================
function [triC,tri3,tri,sll]=Sasa_make_trill
%====================================================
Fid=fopen('/home/sasajima/Dropbox/yellow/Kanto_bound2.txt','r');
bound=textscan(Fid,'%f %f');
fclose(Fid);
bound=cell2mat(bound);
%====================================================
Fid=fopen('/home/sasajima/Dropbox/yellow/PHdepth20130906.txt','r');
dep_main=textscan(Fid,'%f %f %f');
fclose(Fid);
dep_main=cell2mat(dep_main);
%====================================================
Fid=fopen('/home/sasajima/Dropbox/yellow/sg_collect.txt','r');
dep_sub=textscan(Fid,'%f %f %f');
fclose(Fid);
dep_sub=cell2mat(dep_sub);
%===================================================
%E=TriScatteredInterp(dep_main(:,1),dep_main(:,2),dep_main(:,3),'natural');
%F=TriScatteredInterp(dep_sub(:,1),dep_sub(:,2),dep_sub(:,3),'natural');
%min_lon=min(bound(:,1)); max_lon=max(bound(:,1));
%min_lat=min(bound(:,2)); max_lat=max(bound(:,2));
%figure
%plot(bound(:,1),bound(:,2),'r')
%hold on
%n=0;
%int_mesh=3000;
%while n<int_mesh
%  slat=(max_lat-min_lat).*rand(1)+min_lat;
%  slon=(max_lon-min_lon).*rand(1)+min_lon;
%  ID=inpolygon(slon,slat,bound(:,1),bound(:,2));
%  if ID==1
%    n=n+1;
%    s.lat(n)=slat;
%    s.lon(n)=slon;
%s.dep(n)=F(slon,slat)-E(slon,slat);
%if rem(n,round(int_mesh/10))==1;
%  plot3(s.lon,s.lat,s.dep,'.')
%  pause(.1)
%end
%  end
%end
%====================================================
Fid=fopen('/home/sasajima/Dropbox/yellow/Kanto_mesh2.dat','r');
Tll=textscan(Fid,'%f %f %f');
fclose(Fid);
Tll=cell2mat(Tll);
s.lon(:,1)=Tll(:,1);
s.lat(:,1)=Tll(:,2);
sll=Tll;
%ll=[s.lon,s.lat];
%plot(s.lon,s.lat,'.')
%====================================================
tri = delaunay(s.lon,s.lat);
%====================================================
ntri=length(tri);
%nn=0;
%Fid=fopen('/home/sasajima/DATA/XYtriC.dat','w');
%Fid2=fopen('/home/sasajima/DATA/XYtri3.dat','w');
E=TriScatteredInterp(dep_main(:,1),dep_main(:,2),dep_main(:,3),'natural');
F=TriScatteredInterp(dep_sub(:,1),dep_sub(:,2),dep_sub(:,3),'natural');



%tri3=zeros(ntri,3,3);
%triC=zeros(ntri,3);
nn=0;
for n=1:ntri
    
    slon1=s.lon(tri(n,1));
    slat1=s.lat(tri(n,1));
    sdep1=F(slon1,slat1)+E(slon1,slat1);
    
    slon2=s.lon(tri(n,2));
    slat2=s.lat(tri(n,2));
    sdep2=F(slon2,slat2)+E(slon2,slat2);
    
    slon3=s.lon(tri(n,3));
    slat3=s.lat(tri(n,3));
    sdep3=F(slon3,slat3)+E(slon3,slat3);
    
    if sdep1<0.1
        sdep1=0;
    else
    end
    
    if sdep2<0.1
        sdep2=0;
    else
    end
    
    if sdep3<0.1
        sdep3=0;
    else
    end
    
    glon=(slon1+slon2+slon3)/3;
    glat=(slat1+slat2+slat3)/3;
    gdep=(sdep1+sdep2+sdep3)/3;
    
    IN=inpolygon(glon,glat,bound(:,1),bound(:,2));
    
    if gdep==0
    elseif IN==0
    else
        nn=nn+1;
        triC(nn,:)=[glon,glat,gdep];
        tri3(nn,1,1)=slon1;
        tri3(nn,1,2)=slon2;
        tri3(nn,1,3)=slon3;
        tri3(nn,2,1)=slat1;
        tri3(nn,2,2)=slat2;
        tri3(nn,2,3)=slat3;
        tri3(nn,3,1)=sdep1;
        tri3(nn,3,2)=sdep2;
        tri3(nn,3,3)=sdep3;
    end
end

%fprintf(Fid,'%5d %11.5f %11.5f %11.5f\n',n,gx,gy,gdep);
%fprintf(Fid2,'%5d %11.5f %11.5f %11.5f\n',n,sx1,sy1,sz1,sx2,sy2,sz2,sx3,sy3,sz3);
%if ID==1
%  nn=nn+1;
%  s.tri(nn,:)=tri(n,:);
%end
%whos

%whos
%figure(47);
%plot(bound(:,1),bound(:,2),'r')
%hold on
%fclose(Fid);
%fclose(Fid2);
end
%% calcUxyz.m
function calcUxyz(xyz,sUxyz,dUxyz,Slip)
%by Ryohei Sasajima
%last 2014/02/02
m=size(xyz,1);
n=size(Slip,1);

G(1:m,1:n)=sUxyz(1:m,1:n,1);
G(1:m,n+1:2*n)=dUxyz(1:m,1:n,1);

G(m+1:2*m,1:n)=sUxyz(1:m,1:n,2);
G(m+1:2*m,n+1:2*n)=sUxyz(1:m,1:n,2);

G(2*m+1:3*m,1:n)=sUxyz(1:m,1:n,3);
G(2*m+1:3*m,n+1:2*n)=sUxyz(1:m,1:n,3);

SSlip(1:n,1)=Slip(1:n,1);
SSlip(n+1:2*n,1)=Slip(1:n,2);

Uxyz(1:3*m,1)=G.*SSlip;

Fid1=fopen
Fid2=fopen

for i=1:m
    fprintf(Fid1,'%12.6f %12.6f %12.6f %12.6f',xyz(i,1),xyz(i,2),Uxyz(i,1),Uxyz(m+i,1));
    fprintf(Fid2,'%12.6f %12.6f',xyz(i,1),xyz(i,2),Uxyz(2*m+i,1));
end

end
%% circleSV.m
function circleSV

Fid3=fopen('/home/sasajima/Dropbox/yellow/addPAC_aft311.txt','r');
tmeca3=textscan(Fid3,'%f %f %f %f %f %f %f %f %f %f %f %f %f');
fclose(Fid3);
meca3=cell2mat(tmeca3);

Fid4=fopen('/home/sasajima/Dropbox/yellow/addPAC_before311.txt','r');
tmeca4=textscan(Fid4,'%f %f %f %f %f %f %f %f %f %f %f %f %f');
fclose(Fid4);
meca4=cell2mat(tmeca4);

Fid1=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_th_befoMw9.txt','r');
tmeca1=textscan(Fid1,'%f %f %f %f %f %f %f %f %f %f %f %f %f');
fclose(Fid1);
meca1=cell2mat(tmeca1);

Fid2=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_th_aftMw9.txt','r');
tmeca2=textscan(Fid2,'%f %f %f %f %f %f %f %f %f %f %f %f %f');
fclose(Fid2);
meca2=cell2mat(tmeca2);

Fid5=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEwSVall.txt','w');

Fid6=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEnwSVall.txt','w');

Fid10=fopen('/home_tmp/sasajima/DATA/SVmesh.dat','r');
tmblack=textscan(Fid10,'%f %f');
fclose(Fid10);
mblack=cell2mat(tmblack);

Fid11=fopen('/home/sasajima/Dropbox/yellow/PAC_OKH_thrust.txt','r');
tbound=textscan(Fid11,'%f %f');
fclose(Fid11);
bound=cell2mat(tbound);

rr=size(mblack,1);
s=0;

for r=1:rr;
    
    Red=inpolygon(mblack(r,1),mblack(r,2),bound(:,1),bound(:,2));
    
    if Red==0
    else
        s=s+1
        meshll(s,1)=mblack(r,1);
        meshll(s,2)=mblack(r,2);
    end
    
end%for r=1:rr

si1=size(meca1,1);
si2=size(meca2,1);
si3=size(meca3,1);
si4=size(meca4,1);

meca(1:si1,1:13)=meca1;
meca((si1+1):(si1+si2),1:13)=meca2;
meca((si1+si2+1):(si1+si2+si3),1:13)=meca3;
meca((si1+si2+si3+1):(si1+si2+si3+si4),1:13)=meca4;

n=size(meca,1)

nj=size(meshll,1)

ll(1:nj,1:2)=meshll(1:nj,1:2);
[xy]=ll2xy(ll);
meshXY=xy;

clearvars ll xy;

ll(1:n,1:2)=meca(1:n,1:2);
[xy]=ll2xy(ll);
mecaXY=xy;

thsh=50;%[km]

for jj=1:nj
    jj
    cumNM0=0;
    cumSV=0;
    NWcumSV=0;
    NWcount=0;
    
    for i=1:n
        
        dist=sqrt((meshXY(jj,1)-mecaXY(i,1))^2+(meshXY(jj,2)-mecaXY(i,2))^2);
        
        if dist<thsh
            
            M0=meca(i,10)^(meca(i,11)-7);
            cumNM0(1,1)=cumNM0(1,1)+M0;
            cumSV(1,1)=cumSV(1,1)+(meca(i,4)+90)*M0;
            NWcumSV(1,1)=NWcumSV(1,1)+meca(i,4)+90;
            NWcount(1,1)=NWcount(1,1)+1;
            
        else
        end
        
    end
    
    if cumNM0==0;
        cumNM0=100;
    else
    end
    
    if NWcount==0;
        NWcount=1;
    else
    end
    
    avewSV(1,1)=cumSV(1,1)/cumNM0(1,1);%clockwise from North, (toward trench-->);
    avenwSV(1,1)=NWcumSV(1,1)/NWcount(1,1);%non weighted
    
    fprintf(Fid5,'%11.3f %11.3f %11.3f\n',meshll(jj,1),meshll(jj,2),avewSV(1,1));
    fprintf(Fid6,'%11.3f %11.3f %11.3f\n',meshll(jj,1),meshll(jj,2),avenwSV(1,1));
    
    
end

end

%% cmtmake_trill.m
%====================================================
function [triC,tri3,tri,sll]=cmtmake_trill
%by Ryohei Sasajima 2013/12/23
%====================================================
%====================================================
Fid0=fopen('/home/sasajima/Dropbox/yellow/PACdepth201312.txt','r');
dep_main=textscan(Fid0,'%f %f %f');
fclose(Fid0);
dep_main=cell2mat(dep_main);
%====================================================
%====================================================
Fid1=fopen('/home/sasajima/Dropbox/yellow/PAC_all.txt','r');
bound=textscan(Fid1,'%f %f');
fclose(Fid1);
bound=cell2mat(bound);


Fid10=fopen('/home_tmp/sasajima/DATA/green.dat','r');
mgreen=textscan(Fid10,'%f %f');
fclose(Fid10);
mgreen=cell2mat(mgreen);
gg=size(mgreen,1);

s=0;

for g=1:gg;
    
    Green=inpolygon(mgreen(g,1),mgreen(g,2),bound(:,1),bound(:,2));
    
    if Green==0
    else
        s=s+1;
        slon(s,1)=mgreen(g,1);
        slat(s,1)=mgreen(g,2);
    end
end

sall=s

sll(:,1:2)=[slon(:,1),slat(:,1)];

%====================================================
tri = delaunay(slon,slat);
%====================================================
ntri=length(tri);

E=TriScatteredInterp(dep_main(:,1),dep_main(:,2),dep_main(:,3),'natural');%depth of interplate -km

nn=0;
for n=1:ntri
    
    lon1=slon(tri(n,1));
    lat1=slat(tri(n,1));
    dep1=-E(lon1,lat1);%+km
    
    lon2=slon(tri(n,2));
    lat2=slat(tri(n,2));
    dep2=-E(lon2,lat2);
    
    lon3=slon(tri(n,3));
    lat3=slat(tri(n,3));
    dep3=-E(lon3,lat3);
    
    glon=(lon1+lon2+lon3)/3;
    glat=(lat1+lat2+lat3)/3;
    gdep=(dep1+dep2+dep3)/3;
    
    if gdep==0
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

end
%% cmtstrike_dip.m
function[sitaS,sitaD,normVec]=cmtstrike_dip(trixyzC,trixyz3)
%whos
%Define strike and dip of fault%
n=size(trixyzC,1);
ca=zeros(n,3);
cb=zeros(n,3);
InormVec=zeros(n,3);
normVec=zeros(n,3);
strikeVec=zeros(n,3);
dipVec=zeros(n,3);
RsitaS=zeros(n,1);
Rdip=zeros(n,1);
sitaS=zeros(n,1);
sitaD=zeros(n,1);

for i=1:n
    ca(i,:)=[trixyz3(i,1,1)-trixyz3(i,1,3),trixyz3(i,2,1)-trixyz3(i,2,3),trixyz3(i,3,1)-trixyz3(i,3,3)];
    cb(i,:)=[trixyz3(i,1,2)-trixyz3(i,1,3),trixyz3(i,2,2)-trixyz3(i,2,3),trixyz3(i,3,2)-trixyz3(i,3,3)];
    InormVec(i,:)              = cross(ca(i,:),cb(i,:),2);% direct to deeper %
    normVec(i,:)              = InormVec(i,:)./(sqrt((InormVec(i,1)).^2+(InormVec(i,2)).^2+(InormVec(i,3).^2)));
    if (normVec(i,3) < 0) % Enforce clockwise circulation
        normVec(i,:)               = -normVec(i,:);
        % [x(2) x(3)]               = swap(x(2), x(3));%???
        % [y(2) y(3)]               = swap(y(2), y(3));%???
        % [z(2) z(3)]               = swap(z(2), z(3));%???
    end
    strikeVec(i,:)              = [-sin(atan2(normVec(i,2),normVec(i,1))),cos(atan2(normVec(i,2),normVec(i,1))),0];%direct to left hand of who toward dip direction
    dipVec(i,:)                 = cross(normVec(i,:), strikeVec(i,:),2);
    
    %if normVec(1)==0 and normVec(2)==0
    % strikeVec=[1 0 0];
    % dipVec   =[0 1 0];
    % normVec  =[0 0 1];
    %end
    
    RsitaS(i)=atan2(strikeVec(i,2),strikeVec(i,1)); %from x-axis to y-axis rotation [rad]%
    %RsitaSY=2.5*pi-RsitaS; %from y-axis clockwise[rad]%
    %sitaSY=RsitaSY*180/pi;
    sitaS(i)=RsitaS(i)-0.5*pi; %from Y-axis(=North) to clockwise [rad]
    
    Rdip(i)=asin(dipVec(i,3));
    %dip=Rdip*180/pi;
    sitaD(i)=Rdip(i);%dip angle[rad]
end
end
%% cmttrill2trixyz.m
function [trixyzC,trixyz3]=cmttrill2trixyz(triC,tri3,sll)
%-------------------
%  PLTXY TRANSFORMS (ALAT,ALONG) TO (X,Y)
%  WHEN ICORD.NE.0  PLTXY MAKES NO CHANGE IN
%  TRANSFORMATION BETWEEN (X,Y) AND (ALAT,ALONG).
%-------------------%lat==>Y, lon==>X,
A=6.378160e3;
E2=6.6944541e-3;
E12=6.7395719e-3;
D=5.72958e1;
RD=1.0/D;
nn=length(triC);

ALON0=144.0
ALAT0=40.0

trixyz3=zeros(nn,3,3);

mm=length(sll);

%{
for i=1:mm
 ALON=sll(i,1);
 ALAT=sll(i,2);

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
 sxyz(i,1:2)=[X,Y];
end
%}

for i=1:nn
    ALON=tri3(i,1,1);
    ALAT=tri3(i,2,1);
    
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
    trixyz3(i,1,1)=[X*1000];
    trixyz3(i,2,1)=[-Y*1000];
    trixyz3(i,3,1)=[tri3(i,3,1)*1000];
end

for i=1:nn
    ALON=tri3(i,1,2);
    ALAT=tri3(i,2,2);
    
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
    trixyz3(i,1,2)=[X*1000];
    trixyz3(i,2,2)=[-Y*1000];
    trixyz3(i,3,2)=[tri3(i,3,2)*1000];
end

for i=1:nn
    ALON=tri3(i,1,3);
    ALAT=tri3(i,2,3);
    
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
    trixyz3(i,1,3)=[X*1000];
    trixyz3(i,2,3)=[-Y*1000];
    trixyz3(i,3,3)=[tri3(i,3,3)*1000];
end

trixyzC=zeros(nn,3);

for i=1:nn
    ALON=triC(i,1);
    ALAT=triC(i,2);
    
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
    trixyzC(i,1)=[X*1000];
    trixyzC(i,2)=[-Y*1000];
    trixyzC(i,3)=[triC(i,3)*1000];
end

end
%====================================================
%% FSlip2xyll.m
function [xyFSlip]=FSlip2xyll(F_Slip,sitaS)
xyFSlip(:,1)=F_Slip(:,1).*cos(sitaS(:))+F_Slip(:,2).*sin(sitaS(:));
xyFSlip(:,2)=-F_Slip(:,1).*sin(sitaS(:))+F_Slip(:,2).*cos(sitaS(:));
end
%====================================================
%% copygroupCMT.m
function [n]=groupCMT(NARFA,sitaD,sitaS,triC)


Fid0=fopen('/home/sasajima/Dropbox/yellow/PACdepth201312.txt','r');
dep_main=textscan(Fid0,'%f %f %f');
fclose(Fid0);
dep_main=cell2mat(dep_main);

E=TriScatteredInterp(dep_main(:,1),dep_main(:,2),dep_main(:,3),'linear');%depth of interplate -km
F=TriScatteredInterp(triC(:,1),triC(:,2),sitaS(:,1),'linear');
G=TriScatteredInterp(triC(:,1),triC(:,2),sitaD(:,1),'linear');
H=TriScatteredInterp(triC(:,1),triC(:,2),NARFA(:,3),'linear');

Fid1=fopen('/home/sasajima/Dropbox/yellow/NEall.txt','r');
tmeca=textscan(Fid1,'%f %f %f %f %f %f %f %f %f %f %f %f %f');
fclose(Fid1);
meca=cell2mat(tmeca);

Fid00=fopen('/home/sasajima/Dropbox/yellow/NEallDate.txt','r');

Fid2=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_thrust97_11Mw9.txt','w');
Fid3=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEhang97_11Mw9.txt','w');
Fid4=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_intra_upp97_11Mw9.txt','w');
Fid5=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_intra_bot97_11Mw9.txt','w');

Fid7=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_thrust11Mw9_140215.txt','w');
Fid8=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEhang11Mw9_140215.txt','w');
Fid9=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_intra_upp11Mw9_140215.txt','w');
Fid10=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_intra_bot11Mw9_140215.txt','w');

n=size(meca,1);
tt=0;
RD=0;

for i=1:n
    i
    whos
    tline=fgetl(Fid00);
    tt=tt+1;
    
    RD=RD+1;
    serialTime(1,1)=datenum(tline(1:16),'yyyy/mm/dd,HH:MM');
    
    d311='2011 03 11 14:46:00';
    t311=datenum(d311,'yyyy mm dd HH:MM:SS');
    
    lon=meca(i,1);
    lat=meca(i,2);
    focaldep=meca(i,6);
    focalstrike=meca(i,7);
    focaldip=meca(i,8);
    focalSV=meca(i,4)+90;
    
    thrustdep=-E(lon,lat);
    thruststrike=F(lon,lat);
    thrustdip=G(lon,lat);
    thrustSV=H(lon,lat)*180/pi;
    
    if serialTime(1,1)<t311;
        
        if (focaldep>(thrustdep-10))&&(focaldep<(thrustdep+10))&&(focalstrike>(thruststrike-45))&&(focalstrike<(thruststrike+45))&&(focaldip>(thrustdip-12.5))&&(focaldip<(thrustdip+12.5))&&(focalSV>(thrustSV-22.5))&&(focalSV<(thrustSV+22.5));
            
            fprintf(Fid2,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));
            
        elseif (focaldep<thrustdep);
            
            fprintf(Fid3,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));
            
        elseif (focaldep>thrustdep)&&(focaldep<thrustdep+25);
            
            fprintf(Fid4,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));
            
        else
            
            fprintf(Fid5,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));
            
        end
        
        
        
    else
        
        if (focaldep>(thrustdep-10))&&(focaldep<(thrustdep+10))&&(focalstrike>(thruststrike-45))&&(focalstrike<(thruststrike+45))&&(focaldip>(thrustdip-12.5))&&(focaldip<(thrustdip+12.5))&&(focalSV>(thrustSV-22.5))&&(focalSV<(thrustSV+22.5));
            
            fprintf(Fid7,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));
            
        elseif (focaldep<thrustdep);
            
            fprintf(Fid8,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));
            
        elseif (focaldep>thrustdep)&&(focaldep<thrustdep+25);
            
            fprintf(Fid9,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));
            
        else
            
            fprintf(Fid10,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));
            
        end
        
    end
    
    fclose(Fid2);
    fclose(Fid3);
    fclose(Fid4);
    fclose(Fid5);
    fclose(Fid7);
    fclose(Fid8);
    fclose(Fid9);
    fclose(Fid10);
    fclose(Fid00);
    
end
end
%% defineEllipsXY.m
function [xyF]=defineEllipsXY

%change value
Pcenter=[,];%[x,y] [m]
a=0;%first radius (half radius)[m]
b=0;%second radius (half radius)[m]
digArfa=0;%amplitude, clockwise [digree]
n=90;
%===================================
radArfa=digArfa./180.*pi;%clockwise [radian]

theta=0-2.*pi./n;
%make ellips line
for i=1:n
    theta=theta+2.*pi./n;
    xy1(i,1)=a.*cos(theta);
    xy1(i,2)=b.*cos(theta);
    
    xy2(i,1)=xy1(i,1).*cos(radArfa)+xy1(i,2).*sin(radArfa);
    xy2(i,2)=-xy1(i,1).*sin(radArfa)+xy1(i,2).*cos(radArfa);
    
    xyF(i,1)=xy2(i,1)+Pcenter(1,1);
    xyF(i,2)=xy2(i,2)+Pcenter(1,2);
end

end
%% defineEllipsll.m
function [llF]=defineEllipsll

n=30;%30ko/360digree

theta=0-2.*pi./n;
llF=zeros(7,n,2);

%Odawara%
Pcenter=[139.350,35.300];%[lon,lat] [digree]
a=0.110;%first radius (half radius)[digree(lon)]
b=0.080;%second radius (half radius)[digree(lat)]
digArfa=50.0;%amplitude, clockwise from Notrh [digree]

radArfa=digArfa./180.*pi;%clockwise [radian]

theta=0-2.*pi./n;

%make ellips line
for i=1:n
    theta=theta+2.*pi./n;
    ll1(i,1)=a.*cos(theta);
    ll1(i,2)=b.*sin(theta);
    
    ll2(i,1)=ll1(i,1).*cos(radArfa)+ll1(i,2).*sin(radArfa);
    ll2(i,2)=-ll1(i,1).*sin(radArfa)+ll1(i,2).*cos(radArfa);
    
    llF(1,i,1)=ll2(i,1)+Pcenter(1,1);
    llF(1,i,2)=ll2(i,2)+Pcenter(1,2);
end

%Uraga-Straight%
Pcenter=[139.680,35.140];%[lon,lat] [digree]
a=0.110;%first radius (half radius)[m]
b=0.090;%second radius (half radius)[m]
digArfa=0.0;%amplitude, clockwise [digree]

radArfa=digArfa./180.*pi;%clockwise [radian]

theta=0-2.*pi./n;

%make ellips line
for i=1:n
    theta=theta+2.*pi./n;
    ll1(i,1)=a.*cos(theta);
    ll1(i,2)=b.*sin(theta);
    
    ll2(i,1)=ll1(i,1).*cos(radArfa)+ll1(i,2).*sin(radArfa);
    ll2(i,2)=-ll1(i,1).*sin(radArfa)+ll1(i,2).*cos(radArfa);
    
    llF(2,i,1)=ll2(i,1)+Pcenter(1,1);
    llF(2,i,2)=ll2(i,2)+Pcenter(1,2);
end


%N-Tokyo Bay%
Pcenter=[139.950,35.660];%[lon,lat] [digree]
a=0.050;%first radius (half radius)[m]
b=0.045;%second radius (half radius)[m]
digArfa=0.0;%amplitude, clockwise [digree]

radArfa=digArfa./180.*pi;%clockwise [radian]

theta=0-2.*pi./n;

%make ellips line
for i=1:n
    theta=theta+2.*pi./n;
    ll1(i,1)=a.*cos(theta);
    ll1(i,2)=b.*sin(theta);
    
    ll2(i,1)=ll1(i,1).*cos(radArfa)+ll1(i,2).*sin(radArfa);
    ll2(i,2)=-ll1(i,1).*sin(radArfa)+ll1(i,2).*cos(radArfa);
    
    llF(3,i,1)=ll2(i,1)+Pcenter(1,1);
    llF(3,i,2)=ll2(i,2)+Pcenter(1,2);
end


%off-Tateyama%
Pcenter=[140.160,34.680];%[lon,lat] [digree]
a=0.090;%first radius (half radius)[m]
b=0.055;%second radius (half radius)[m]
digArfa=0.000;%amplitude, clockwise [digree]

radArfa=digArfa./180.*pi;%clockwise [radian]

theta=0-2.*pi./n;

%make ellips line
for i=1:n
    theta=theta+2.*pi./n;
    ll1(i,1)=a.*cos(theta);
    ll1(i,2)=b.*sin(theta);
    
    ll2(i,1)=ll1(i,1).*cos(radArfa)+ll1(i,2).*sin(radArfa);
    ll2(i,2)=-ll1(i,1).*sin(radArfa)+ll1(i,2).*cos(radArfa);
    
    llF(4,i,1)=ll2(i,1)+Pcenter(1,1);
    llF(4,i,2)=ll2(i,2)+Pcenter(1,2);
end


%S-Boso%
Pcenter=[139.970,34.920];%[lon,lat] [digree]
a=0.200;%first radius (half radius)[m]
b=0.130;%second radius (half radius)[m]
digArfa=15.0;%amplitude, clockwise [digree]

radArfa=digArfa./180.*pi;%clockwise [radian]

theta=0-2.*pi./n;

%make ellips line
for i=1:n
    theta=theta+2.*pi./n;
    ll1(i,1)=a.*cos(theta);
    ll1(i,2)=b.*sin(theta);
    
    ll2(i,1)=ll1(i,1).*cos(radArfa)+ll1(i,2).*sin(radArfa);
    ll2(i,2)=-ll1(i,1).*sin(radArfa)+ll1(i,2).*cos(radArfa);
    
    llF(5,i,1)=ll2(i,1)+Pcenter(1,1);
    llF(5,i,2)=ll2(i,2)+Pcenter(1,2);
end


%SE-off-Boso%
Pcenter=[140.595,34.700];%[lon,lat] [digree]
a=0.200;%first radius (half radius)[m]
b=0.095;%second radius (half radius)[m]
digArfa=22.5;%amplitude, clockwise [digree]

radArfa=digArfa./180.*pi;%clockwise [radian]

theta=0-2.*pi./n;

%make ellips line
for i=1:n
    theta=theta+2.*pi./n;
    ll1(i,1)=a.*cos(theta);
    ll1(i,2)=b.*sin(theta);
    
    ll2(i,1)=ll1(i,1).*cos(radArfa)+ll1(i,2).*sin(radArfa);
    ll2(i,2)=-ll1(i,1).*sin(radArfa)+ll1(i,2).*cos(radArfa);
    
    llF(6,i,1)=ll2(i,1)+Pcenter(1,1);
    llF(6,i,2)=ll2(i,2)+Pcenter(1,2);
end


%S-Kujukuri%
Pcenter=[140.500,35.075];%[lon,lat] [digree]
a=0.230;%first radius (half radius)[m]
b=0.300;%second radius (half radius)[m]
digArfa=345.0;%amplitude, clockwise [digree]

radArfa=digArfa./180.*pi;%clockwise [radian]

theta=0-2.*pi./n;

%make ellips line
for i=1:n
    theta=theta+2.*pi./n;
    ll1(i,1)=a.*cos(theta);
    ll1(i,2)=b.*sin(theta);
    
    ll2(i,1)=ll1(i,1).*cos(radArfa)+ll1(i,2).*sin(radArfa);
    ll2(i,2)=-ll1(i,1).*sin(radArfa)+ll1(i,2).*cos(radArfa);
    
    llF(7,i,1)=ll2(i,1)+Pcenter(1,1);
    llF(7,i,2)=ll2(i,2)+Pcenter(1,2);
end


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
   define(i,1)=[1];
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
%% defineSlipBoso.m
function [Slip]=defineSlipBoso(triC,sitaS,llF,Velo,SVxy)
%A1=Rectangle asperity_1%
%p=minimum X and max Y, q=max X and minimum Y%
%y=+ wa south%


n=size(triC,1)
Slip=zeros(n,2);
define=zeros(n,1);

llFa(:,1:2)=llF(1,:,1:2);

for i=1:n
    %ID1=inpolygon(triC(i,1),triC(i,2),llF(1,:,1),llF(1,:,2));
    ID1=inpolygon(triC(i,1),triC(i,2),llFa(:,1),llFa(:,2));
    if ID1==1
        A1sSlip1=-cos(sitaS(i)).*Velo(1,1).*sin(SVxy(1,1))-sin(sitaS(i)).*Velo(1,1).*cos(SVxy(1,1));
        A1dSlip1=sin(sitaS(i)).*Velo(1,1).*sin(SVxy(1,1))+cos(sitaS(i)).*Velo(1,1).*cos(SVxy(1,1));
        Slip(i,:)=[A1sSlip1,A1dSlip1];
        define(i,1)=[1];%1=asperity,0=frictionless%
    end
    
    ID2=inpolygon(triC(i,1),triC(i,2),llF(2,:,1),llF(2,:,2));
    if ID2==1
        A1sSlip2=-cos(sitaS(i)).*Velo(2,1).*sin(SVxy(2,1))-sin(sitaS(i)).*Velo(2,1).*cos(SVxy(2,1));
        A1dSlip2=sin(sitaS(i)).*Velo(2,1).*sin(SVxy(2,1))+cos(sitaS(i)).*Velo(2,1).*cos(SVxy(2,1));
        Slip(i,:)=[A1sSlip2,A1dSlip2];
        define(i,1)=[1];
    end
    
    ID3=inpolygon(triC(i,1),triC(i,2),llF(3,:,1),llF(3,:,2));
    if ID3==1
        A1sSlip3=-cos(sitaS(i)).*Velo(3,1).*sin(SVxy(3,1))-sin(sitaS(i)).*Velo(3,1).*cos(SVxy(3,1));
        A1dSlip3=sin(sitaS(i)).*Velo(3,1).*sin(SVxy(3,1))+cos(sitaS(i)).*Velo(3,1).*cos(SVxy(3,1));
        Slip(i,:)=[A1sSlip3,A1dSlip3];
        define(i,1)=[1];
    end
    
    ID4=inpolygon(triC(i,1),triC(i,2),llF(4,:,1),llF(4,:,2));
    if ID4==1
        A1sSlip4=-cos(sitaS(i)).*Velo(4,1).*sin(SVxy(4,1))-sin(sitaS(i)).*Velo(4,1).*cos(SVxy(4,1));
        A1dSlip4=sin(sitaS(i)).*Velo(4,1).*sin(SVxy(4,1))+cos(sitaS(i)).*Velo(4,1).*cos(SVxy(4,1));
        Slip(i,:)=[A1sSlip4,A1dSlip4];
        define(i,1)=[1];
    end
    
    ID5=inpolygon(triC(i,1),triC(i,2),llF(5,:,1),llF(5,:,2));
    if ID5==1
        A1sSlip5=-cos(sitaS(i)).*Velo(5,1).*sin(SVxy(5,1))-sin(sitaS(i)).*Velo(5,1).*cos(SVxy(5,1));
        A1dSlip5=sin(sitaS(i)).*Velo(5,1).*sin(SVxy(5,1))+cos(sitaS(i)).*Velo(5,1).*cos(SVxy(5,1));
        Slip(i,:)=[A1sSlip5,A1dSlip5];
        define(i,1)=[1];
    end
    
    ID6=inpolygon(triC(i,1),triC(i,2),llF(6,:,1),llF(6,:,2));
    if ID6==1
        A1sSlip6=-cos(sitaS(i)).*Velo(6,1).*sin(SVxy(6,1))-sin(sitaS(i)).*Velo(6,1).*cos(SVxy(6,1));
        A1dSlip6=sin(sitaS(i)).*Velo(6,1).*sin(SVxy(6,1))+cos(sitaS(i)).*Velo(6,1).*cos(SVxy(6,1));
        Slip(i,:)=[A1sSlip6,A1dSlip6];
        define(i,1)=[1];
    end
    
    ID7=inpolygon(triC(i,1),triC(i,2),llF(7,:,1),llF(7,:,2));
    if ID7==1
        A1sSlip7=-cos(sitaS(i)).*Velo(7,1).*sin(SVxy(7,1))-sin(sitaS(i)).*Velo(7,1).*cos(SVxy(7,1));
        A1dSlip7=sin(sitaS(i)).*Velo(7,1).*sin(SVxy(7,1))+cos(sitaS(i)).*Velo(7,1).*cos(SVxy(7,1));
        Slip(i,:)=[A1sSlip7,A1dSlip7];
        define(i,1)=[1];
    end
end

end

%% diviSV.m
function diviSV

Fid3=fopen('/home/sasajima/Dropbox/yellow/addPAC_aft311.txt','r');
tmeca3=textscan(Fid3,'%f %f %f %f %f %f %f %f %f %f %f %f %f');
fclose(Fid3);
meca3=cell2mat(tmeca3);

Fid4=fopen('/home/sasajima/Dropbox/yellow/addPAC_before311.txt','r');
tmeca4=textscan(Fid4,'%f %f %f %f %f %f %f %f %f %f %f %f %f');
fclose(Fid4);
meca4=cell2mat(tmeca4);

Fid1=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_thrust97_11Mw9.txt','r');
tmeca1=textscan(Fid1,'%f %f %f %f %f %f %f %f %f %f %f %f %f');
fclose(Fid1);
meca1=cell2mat(tmeca1);

Fid2=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_thrust11Mw9_140215.txt','r');
tmeca2=textscan(Fid2,'%f %f %f %f %f %f %f %f %f %f %f %f %f');
fclose(Fid2);
meca2=cell2mat(tmeca2);

si1=size(meca1,1);
si2=size(meca2,1);
si3=size(meca3,1);
si4=size(meca4,1);

meca(1:si1,1:13)=meca1;
meca((si+1):(si1+si2),1:13)=meca2;
meca((si1+si2+1):(si1+si2+si3),1:13)=meca3;
meca((si1+si2+si3+1):(si1+si2+si3+si4),1:13)=meca4;

n=size(meca,1)

nj=27;

cumNM0=zeros(nj,1)
cumSV=zeros(nj,1)
NWcumSV=zeros(nj,1)
NWcount=zeros(nj,1)

jj=0;

for j=36:49:0.5
    
    jj=jj+1;
    
    for i=1:n
        
        if (meca(i,2)>j)&&(meca(i,2)<j+1)
            
            M0=meca(i,10)^(meca(i,11)-7);
            cumNM0(jj,1)=cumNM0(jj,1)+M0;
            cumSV(jj,1)=cumSV(jj,1)+(meca(i,4)+90)*M0;
            NWcumSV(jj,1)=NWcumSV(jj,1)+meca(i,4)+90;
            NWcount(jj,1)=NWconut(jj,1)+1;
        else
        end
        
    end
    
    avewSV(1:nj,1)=cumSV(1:nj,1)./cumNM0(1:nj,1);%clockwise from North, (toward trench-->);
    avenwSV(1:nj,1)=NWcumSV(1:nj,1)./NWcount(1:nj,1);%non weighted
    
    fprintf
    fprintf
    
end

end

%% euler2vel.m
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
%% grouping_CMT.m
function GroupingCMT
%for grouping CMT mechanism of earthquakes
%by Ryohei Sasajima
%final 2013/12/23


%[triC,tri3,tri,sll]=cmtmake_trill;
%save('/home_tmp/sasajima/DATA/PassiveSlip/PAC_cmt/tri','triC','tri3','tri','sll');
load('/home_tmp/sasajima/DATA/PassiveSlip/PAC_cmt/tri','triC','tri3','tri','sll');


%[trixyzC,trixyz3]=cmttrill2trixyz(triC,tri3,sll);
%save('/home_tmp/sasajima/DATA/PassiveSlip/PAC_cmt/trixyz','trixyzC','trixyz3');
load('/home_tmp/sasajima/DATA/PassiveSlip/PAC_cmt/trixyz','trixyzC','trixyz3');

%trixyzC(1:10,:)
%trixyz3(1:10,:,:)
%pause

%[sitaS,sitaD,normVec]=cmtstrike_dip(trixyzC,trixyz3);
%save('/home_tmp/sasajima/DATA/PassiveSlip/PAC_cmt/sita','sitaS','sitaD','normVec');
load('/home_tmp/sasajima/DATA/PassiveSlip/PAC_cmt/sita','sitaS','sitaD','normVec');


ll(:,1:2)=triC(:,1:2);

%sitaS(1:10,1)
%sitaD(20:30,1)


[vnorm,v,NARFA]=SasaEular2Velo(ll);

[n]=groupCMT(NARFA,sitaS,sitaD,triC);

end

%% least_pole.m
function least_pole
Fid1=fopen('/home/sasajima/DATA/PH-IA-lpdF.dat','r');
%Fid2=fopen('/home/sasajima/DATA/OUTPUT/PH-IH-leastpole.dat','w');
dF=textscan(Fid1,'%f %f %f %f %f %f %f %f %f');
dF=cell2mat(dF);
fclose(Fid1);

p0=7958.68798828;
B=[dF(:,3),dF(:,4),dF(:,5)];
A=[dF(:,1),dF(:,2)];
I=eye(3);
%P=[dF(:,9)/p0];
%P=diag(P);
%iP=inv(P);
%M=B*iP*B';
whos
M=B*I*B';
iM=inv(M);
w=[dF(:,6)];
K=A'*iM*A;
iK=inv(K);
whos
dlt=(-iK)*A'*iM*w
v=(I*B'*iM*A*iK*A'*iM*(-I)*B'*iM)*w;
whos
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
function [sSsn,sSdn,dSsn,dSdn]=loadMAT(trixyzC);

rootname='/home_tmp/sasajima/DATA/GreenF/PACTohoku2PACTohoku/';
ss='sSsn/sSsn';
sd='sSdn/sSdn';
ds='dSsn/dSsn';
dd='dSdn/dSdn';
extension='.dat';

n=length(trixyzC);

for i=1:n;
    
    loatMAT=i
    
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
function [sUxyz,dUxyz]=loadMAT2(xyz,trixyz3);

rootname='/home_tmp/sasajima/DATA/GreenF/PAC2test/';
sU='sU';
dU='dU';
extension='.dat';

m=size(xyz,1);
n=size(trixyz3,1);

for i=1:n;
    
    loadMAT=i
    
    w=num2str(i);
    
    filename1= [rootname,sU,w,extension];
    filename2= [rootname,dU,w,extension];
    
    load(filename1,'sUxyzi','-mat');
    load(filename2,'dUxyzi','-mat');
    
    sUxyz(:,i,1:3)=sSsni(:,1:3);
    dUxyz(:,i,1:3)=dSdni(:,1);
    
end
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

rootname='/home_tmp/sasajima/DATA/GreenF/PAC2test/'
dU='dU';
extension='.dat';

for i=1:n%numbers of fault%
    
    di=i
    w=num2str(i);
    
    x=[trixyz3(i,1,1),trixyz3(i,1,2),trixyz3(i,1,3)];
    y=[trixyz3(i,2,1),trixyz3(i,2,2),trixyz3(i,2,3)];
    z=[trixyz3(i,3,1),trixyz3(i,3,2),trixyz3(i,3,3)];
    
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

rootname='/home_tmp/sasajima/DATA/GreenF/PACtest2PACtest/'
dds='dSsn/dSsn';
ddn='dSdn/dSdn';
extension='.dat';

for i=1:n%numbers of fault%
    
    di=i
    w=num2str(i);
    
    x=[trixyz3(i,1,1),trixyz3(i,1,2),trixyz3(i,1,3)];
    y=[trixyz3(i,2,1),trixyz3(i,2,2),trixyz3(i,2,3)];
    z=[trixyz3(i,3,1),trixyz3(i,3,2),trixyz3(i,3,3)];
    
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
    dSsn(:,i)=-s2(i).*(c1(i).*s1(i).*(Syy(:,1)-Sxx(:,1))+((c1(i)).^2-(s1(i)).^2).*Sxy(:,1))+c2(i).*(c1(i).*Sxz(:,1)+s1(i).*Syz(:,1));
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

rootname='/home_tmp/sasajima/DATA/GreenF/PAC2test/'
sU='sU';
extension='.dat';

for i=1:n%numbers of fault%
    
    si=i
    w=num2str(i);
    
    x=[trixyz3(i,1,1),trixyz3(i,1,2),trixyz3(i,1,3)];
    y=[trixyz3(i,2,1),trixyz3(i,2,2),trixyz3(i,2,3)];
    z=[trixyz3(i,3,1),trixyz3(i,3,2),trixyz3(i,3,3)];
    
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


rootname='/home_tmp/sasajima/DATA/GreenF/PACtest2PACtest/'
ssn='sSsn/sSsn';
sdn='sSdn/sSdn';
extension='.dat';

for i=1:n%numbers of fault%
    
    si=i
    w=num2str(i);
    
    x=[trixyz3(i,1,1),trixyz3(i,1,2),trixyz3(i,1,3)];
    y=[trixyz3(i,2,1),trixyz3(i,2,2),trixyz3(i,2,3)];
    z=[trixyz3(i,3,1),trixyz3(i,3,2),trixyz3(i,3,3)];
    
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

for bl=1:blbl;
    
    Blue=inpolygon(mblack(bl,1),mblack(bl,2),PACblue(:,1),PACblue(:,2));
    
    if Blue==0
    else
        s=s+1;
        slon(s,1)=mblack(bl,1);
        slat(s,1)=mblack(bl,2);
    end
end

s

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
%nn=0;
%Fid=fopen('/home/sasajima/DATA/XYtriC.dat','w');
%Fid2=fopen('/home/sasajima/DATA/XYtri3.dat','w');
E=TriScatteredInterp(dep_main(:,1),dep_main(:,2),dep_main(:,3),'natural');
F=TriScatteredInterp(dep_sub(:,1),dep_sub(:,2),dep_sub(:,3),'natural');

%tri3=zeros(ntri,3,3);
%triC=zeros(ntri,3);
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
    %{
   if dep1<0.1
     dep1=0;
   else
   end
  
   if dep2<0.1
     dep2=0;
   else
   end

   if dep3<0.1
     dep3=0;
   else
   end
    %}
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

Fid6=fopen('/home_tmp/sasajima/DATA/PAC_tri3.txt','w');

for w=1:nn;
    for u=1:3;
        fprintf(Fid6,'%10.4f %9.4f %10.4f\n',tri3(w,1,u),tri3(w,2,u),tri3(w,3,u));
    end
end

fclose(Fid6);

%fprintf(Fid7,'%5d %11.5f %11.5f %11.5f\n',n,gx,gy,gdep);
%fprintf(Fid8,'%5d %11.5f %11.5f %11.5f\n',n,sx1,sy1,sz1,sx2,sy2,sz2,sx3,sy3,sz3);
%if ID==1
%  nn=nn+1;
%  s.tri(nn,:)=tri(n,:);
%end
%whos

%whos
%figure(47);
%plot(bound(:,1),bound(:,2),'r')
%hold on
%fclose(Fid);
%fclose(Fid2);
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

for r=1:rr;
    
    Red=inpolygon(mred(r,1),mred(r,2),PACred(:,1),PACred(:,2));
    
    if Red==0
    else
        s=s+1;
        slon(s,1)=mred(r,1);
        slat(s,1)=mred(r,2);
    end
end

sred=s


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

sll(:,1:2)=[slon(:,1),slat(:,1)];

%ll=[s.lon,s.lat];
%plot(s.lon,s.lat,'.')
%====================================================
tri = delaunay(slon,slat);
%====================================================
ntri=length(tri);
%nn=0;
%Fid=fopen('/home/sasajima/DATA/XYtriC.dat','w');
%Fid2=fopen('/home/sasajima/DATA/XYtri3.dat','w');
E=TriScatteredInterp(dep_main(:,1),dep_main(:,2),dep_main(:,3),'natural');%depth of interplate -km
F=TriScatteredInterp(dep_sub(:,1),dep_sub(:,2),dep_sub(:,3),'natural');%depth of seafloor -km

%tri3=zeros(ntri,3,3);
%triC=zeros(ntri,3);
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
    %{
   if dep1<0.1
     dep1=0;
   else
   end
  
   if dep2<0.1
     dep2=0;
   else
   end

   if dep3<0.1
     dep3=0;
   else
   end
    %}
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

Fid6=fopen('/home_tmp/sasajima/DATA/PAC_tri3.txt','w');

for w=1:nn;
    for u=1:3;
        fprintf(Fid6,'%10.4f %9.4f %10.4f\n',tri3(w,1,u),tri3(w,2,u),tri3(w,3,u));
    end
end

fclose(Fid6);

%fprintf(Fid7,'%5d %11.5f %11.5f %11.5f\n',n,gx,gy,gdep);
%fprintf(Fid8,'%5d %11.5f %11.5f %11.5f\n',n,sx1,sy1,sz1,sx2,sy2,sz2,sx3,sy3,sz3);
%if ID==1
%  nn=nn+1;
%  s.tri(nn,:)=tri(n,:);
%end
%whos

%whos
%figure(47);
%plot(bound(:,1),bound(:,2),'r')
%hold on
%fclose(Fid);
%fclose(Fid2);
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
xyz01=cell2mat(xyz01);

xyz=xyz01;

end
%====================================================
%% Slip2xyll.m
function [xySlip]=Slip2xyll(Slip,sitaS)
xySlip(:,1)=Slip(:,1).*cos(sitaS(:))+Slip(:,2).*sin(sitaS(:));
xySlip(:,2)=-Slip(:,1).*sin(sitaS(:))+Slip(:,2).*cos(sitaS(:));
end
%====================================================
%% smoothSV.m
function smoothSV(circle)

threshold=0;

n=size(circle,2);

for i=1:n
    
    IN==inpolygon(meca(:,1),meca(:,2),circle(:,1),circle(:,2));
    
    if IN==1;
    else
    end
    
end
end

%% strike_dip.m
function[sitaS,sitaD,normVec]=strike_dip(trixyzC,trixyz3)
%whos
%Define strike and dip of fault%
n=size(trixyzC,1);
ca=zeros(n,3);
cb=zeros(n,3);
InormVec=zeros(n,3);
normVec=zeros(n,3);
strikeVec=zeros(n,3);
dipVec=zeros(n,3);
RsitaS=zeros(n,1);
Rdip=zeros(n,1);
sitaS=zeros(n,1);
sitaD=zeros(n,1);

for i=1:n
    ca(i,:)=[trixyz3(i,1,1)-trixyz3(i,1,3),trixyz3(i,2,1)-trixyz3(i,2,3),trixyz3(i,3,1)-trixyz3(i,3,3)];
    cb(i,:)=[trixyz3(i,1,2)-trixyz3(i,1,3),trixyz3(i,2,2)-trixyz3(i,2,3),trixyz3(i,3,2)-trixyz3(i,3,3)];
    InormVec(i,:)              = cross(ca(i,:),cb(i,:),2);% direct to deeper %
    normVec(i,:)              = InormVec(i,:)./(sqrt((InormVec(i,1)).^2+(InormVec(i,2)).^2+(InormVec(i,3).^2)));
    if (normVec(i,3) < 0) % Enforce clockwise circulation
        normVec(i,:)               = -normVec(i,:);
        % [x(2) x(3)]               = swap(x(2), x(3));%???
        % [y(2) y(3)]               = swap(y(2), y(3));%???
        % [z(2) z(3)]               = swap(z(2), z(3));%???
    end
    strikeVec(i,:)              = [-sin(atan2(normVec(i,2),normVec(i,1))),cos(atan2(normVec(i,2),normVec(i,1))),0];%direct to left hand of who toward dip direction
    dipVec(i,:)                 = cross(normVec(i,:), strikeVec(i,:),2);
    
    %if normVec(1)==0 and normVec(2)==0
    % strikeVec=[1 0 0];
    % dipVec   =[0 1 0];
    % normVec  =[0 0 1];
    %end
    
    RsitaS(i)=atan2(strikeVec(i,2),strikeVec(i,1)); %from x-axis to y-axis rotation [rad]%
    %RsitaSY=2.5*pi-RsitaS; %from y-axis clockwise[rad]%
    %sitaSY=RsitaSY*180/pi;
    sitaS(i)=RsitaS(i); %from x-axis to y-axis rotation [raad]
    
    Rdip(i)=asin(dipVec(i,3));
    %dip=Rdip*180/pi;
    sitaD(i)=Rdip(i);
end
end
%====================================================
%% vecxy2ll.m
function [FSlipll]=vecxy2ll(trixyzC,xyFSlip)

%     PLTXY TRANSFORMS (X,Y) TO (ALAT,ALONG)
%  WHEN ICORD.NE.0  PLTXY MAKES NO CHANGE IN TRANSFORMATION  BETWEEN
%              (X,Y) AND (ALAT,ALONG).
%global LAT0 LON0

LAT0=38.3;
LON0=142.5;

whos

nn=length(trixyzC);

X(:,1)=trixyzC(:,1);%location%
Y(:,1)=trixyzC(:,2);

A=6.378160e3;
E2=6.6944541e-3;
E12=6.7395719e-3;
D=5.72958e1;
RD=1.0/D;
RLATO = LAT0.*RD;
SLATO = sin(LAT0);
R     = A.*(1-E2)./sqrt((1-E2.*SLATO.^2).^3);
AN    = A./sqrt(1.0-E2.*SLATO.^2);
V2    = 1 + E12.*cos(RLATO).^2;
C1    = D./R;
C2    = D./AN;

PH1   = LAT0+C1.*Y;
RPH1  = PH1.*RD;
TPHI1 = tan(RPH1);
CPHI1 = cos(RPH1);
ALAT0  = PH1-(C2.*X).^2.*V2.*TPHI1./(2.*D);
ALON0  = LON0+C2.*X./CPHI1-(C2.*X).^3.*(1.0+2.*TPHI1.^2)./(6.*D.^2.*CPHI1);

X2(:,1)=X(:,1)+xyFSlip(:,1);
Y2(:,1)=Y(:,1)+xyFSlip(:,2);

PH2   = LAT0+C1.*Y2;
RPH2  = PH2.*RD;
TPHI2 = tan(RPH2);
CPHI2 = cos(RPH2);
ALAT1  = PH2-(C2.*X2).^2.*V2.*TPHI2./(2.*D);
ALON1  = LON0+C2.*X2./CPHI2-(C2.*X).^3.*(1.0+2.*TPHI2.^2)./(6.*D.^2.*CPHI2);

NS=ALAT1-ALAT0;
EW=ALON1-ALON0;

Snorm(:,1)=(xyFSlip(:,1).^2+xyFSlip(:,2).^2).^0.5;%[m]

for j=1:nn;
    SV(j,1)=atan2(EW(j,1),NS(j,1));%Clock wise from North [radian]
    SNS(j,1)=Snorm(j,1).*cos(SV(j,1));%N=plus
    SEW(j,1)=Snorm(j,1).*sin(SV(j,1));%E=plus
end

%SVdig=SV.*180./pi;%[digree]

FSlipll(:,1)=SNS(:,1);
FSlipll(:,2)=SEW(:,1);

end
