function SV_ll2xyz
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
EpiCLON=143;
EpiCLAT=38.2;
PlateV=68.0;%counter clock wise from North (degree)%

%calculation%
RPV=PlateV/180*pi;
AftLON=EpiCLON-cos(RPV);
AftLAT=EpiCLAT+sin(RPV);
Vlon=[EpiCLON;AftLON];
Vlat=[EpiCLAT;AftLAT];

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
SVxyz=atan2(VY(2)-VY(1),VX(2)-VX(1))%counter clock wise from Y-axis (radian)%
SVxyzdig=SVxyz.*180./pi
end
