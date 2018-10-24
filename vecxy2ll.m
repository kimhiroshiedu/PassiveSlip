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

