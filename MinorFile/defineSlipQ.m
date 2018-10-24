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

