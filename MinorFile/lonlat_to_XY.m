function lonlat_to_XY
%-------------------
%  PLTXY TRANSFORMS (ALAT,ALON) TO (X,Y)
%  WHEN ICORD.NE.0  PLTXY MAKES NO CHANGE IN 
%  TRANSFORMATION BETWEEN (X,Y) AND (ALAT,ALON).
%  D=5.72958e1[degree]=1.0[radian],
%  RLAT=ALAT[radian]
%-------------------

Fid=fopen('/home_tmp/sasajima/DATA/Tohoku_xy.dat','w');

%Change value%
ALAT=32.99;
ALAT0=38.5;
ALON0=143.0;
%Constant%
A=6.378160e3;
E2=6.6944541e-3;
E12=6.7395719e-3;
D=5.72958e1;
RD=1.0/D;

 for n=1:1201
  ALAT=ALAT+0.01;
  ALON=138.99;
  
   for m=1:901
    ALON=ALON+0.01;

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
   
   fprintf(Fid,'%9.3f %9.3f %15.7f %15.7f\n',ALAT,ALON,X,Y);

   end
 end
 fclose(Fid);
end
