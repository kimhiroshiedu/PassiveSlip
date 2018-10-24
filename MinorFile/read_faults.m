function [FLOC,DA,STR,DIP,NV,AL,AW,FLOC_GMT]=read_faults(list)
Fid=fopen(list,'r');
NN=0;
while 1
  tline = fgetl(Fid);
  if ~ischar(tline); break; end
  NN=NN+1;
  loc_f=fscanf(Fid,'%f %f %f \n', [3 4])';
  tline = fgetl(Fid);
  [FLOC(NN,:),DA(NN,1),STR(NN,1),DIP(NN,1),NV(NN,:),AL(NN,1),AW(NN,1)]=est_fault_p(loc_f);
  FLOC_GMT(NN,:,:)=loc_f;
end
fclose(Fid);
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
