function [SVxy]=SV_ll2xyz_boso(Vec)
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
%Change Value%
ALON0=142.5;
ALAT0=38.3;

PlateV(:,1:3)=Vec(:,1:3);%counter clock wise from North (degree)%

n=length(Vec);
%Vec =slip vector clockwise from N [digree]

for j=1:n

%calculation%
RPV=PlateV(j,3)/180*pi;
%!!mistake!! AftLON=PlateV(j,1)-cos(RPV);

dLAT=0.1./180.*pi;

 if cos(Vec(j,3))>=0 %NS component of slip vector is plus;

   if sin(Vec(j,3))>=0 %EW component of slip vector is plus;
    radVec(j,1)=Vec(j,3)./180.*pi;%radian
    dPSY=sin(dLAT).*cos(radVec(j,1));%radian
    dLON=
   else
    radVec(j,1)=(360-Vec(j,3))./180.*pi;%radian
    dPSY=sin(dLAT).*cos(radVec(j,1));%radian
   end


    AftLAT=Vec(j,2)+dLAT.*180./pi;%digree
    AftLON=Vec(j,1)+dLON.*180./pi;%digree

 else %NS component of slip vecotr is minus;

   if sin(Vec(j,3))>=0 %EW component of slip vector is plus;
    radVec(j,1)=(180-Vec(j,3))./180.*pi;%radian
   else
    radVec(j,1)=(Vec(j,3)-180)./180.*pi;%radian
   end

    dPSY=sin(dLAT).*cos(radVec(j,1));%radian

    AftLAT=Vec(j,2)+dLAT.*180./pi;%digree
    AftLON=Vec(j,1)+dLON.*180./pi;%digree


!!mistake!! AftLAT=PlateV(j,2)+sin(RPV);
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


