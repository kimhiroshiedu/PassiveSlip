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
