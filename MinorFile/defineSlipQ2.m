function [Slip]=defineSlipQ(triC,sitaS)
%A1=Rectangle asperity_1%
%p=minimum X and max Y, q=max X and minimum Y%
%y=+ wa south%

Fid1=fopen('/home_tmp/sasajima/DATA/VxP4.dat','r');
Fid2=fopne('/home_tmp/sasajima/DATA/PAC-NA-SVxyz.dat','r');
Vx=textscan(Fid1,'%f %f %f');
Vx=cell2mat(Vx); 
SV=textscan(Fid2,'%f %f %f');
SV=cell2mat(SV);

n=size(triC,1);
Slip=zeros(n,2);
define=zeros(n,1);

InterpSV=TriScatteredInterp(Vx(:,1),Vx(:,2),Vx(:,3));
InterpVx=TriScatteredInterp(SV(:,1),SV(:,2),SV(:,3));

%2011 Tohoku-oki%

A1plon=142.8;
A1plat=38.5;
A1qlon=143.5;
A1qlat=37.5;

 for i=1:n
  if triC(i,1)<A1plon,Slip(i,:)=[0,0];
  elseif triC(i,1)>A1qlon,Slip(i,:)=[0,0];
  elseif triC(i,2)<A1qlat,Slip(i,:)=[0,0];
  elseif triC(i,2)>A1plat,Slip(i,:)=[0,0];
  else
   A1sSlip=-cos(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*sin(InterpSV(triC(i,1),triC(i,2)))-sin(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*cos(InterpSV(triC(i,1),triC(i,2)));
   A1dSlip=sin(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*sin(InterpSV(triC(i,1),triC(i,2)))+cos(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*cos(InterpSV(triC(i,1),triC(i,2)));
   Slip(i,:)=[A1sSlip,A1dSlip];
   define(i,1)=[1];
  end
 end

%Hokkaido 500year%

A2plon=145.0;
A2plat=42.2;
A2qlon=146.5;
A2qlat=41.25;

 for i=1:n
  if     define(i,1)==1;
  elseif triC(i,1)<A2plon,Slip(i,:)=[0,0];
  elseif triC(i,1)>A2qlon,Slip(i,:)=[0,0];
  elseif triC(i,2)<A2qlat,Slip(i,:)=[0,0];
  elseif triC(i,2)>A2plat,Slip(i,:)=[0,0];
  else
   A2sSlip=-cos(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*sin(InterpSV(triC(i,1),triC(i,2)))-sin(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*cos(InterpSV(triC(i,1),triC(i,2)));
   A2dSlip=sin(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*sin(InterpSV(triC(i,1),triC(i,2)))+cos(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*cos(InterpSV(triC(i,1),triC(i,2)));
   Slip(i,:)=[A2sSlip,A2dSlip];
   define(i,1)=[1];
  end
 end

%1896 Meiji Sanriku-oki%

A3plon=143.7;
A3plat=39.6;
A3qlon=144.15;
A3qlat=38.8;

 for i=1:n
  if     define(i,1)==1;
  elseif triC(i,1)<A3plon,Slip(i,:)=[0,0];
  elseif triC(i,1)>A3qlon,Slip(i,:)=[0,0];
  elseif triC(i,2)<A3qlat,Slip(i,:)=[0,0];
  elseif triC(i,2)>A3plat,Slip(i,:)=[0,0];
  else
   A3sSlip=-cos(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*sin(InterpSV(triC(i,1),triC(i,2)))-sin(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*cos(InterpSV(triC(i,1),triC(i,2)));
   A3dSlip=sin(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*sin(InterpSV(triC(i,1),triC(i,2)))+cos(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*cos(InterpSV(triC(i,1),triC(i,2)));
   Slip(i,:)=[A3sSlip,A3dSlip];
   define(i,1)=[1];
  end
 end

%2003 Tokachi-oki%
A4plon=143.65;
A4plat=42.2;
A4qlon=144.2;
A4qlat=41.7;

 for i=1:n
  if     define(i,1)==1;
  elseif triC(i,1)<A4plon,Slip(i,:)=[0,0];
  elseif triC(i,1)>A4qlon,Slip(i,:)=[0,0];
  elseif triC(i,2)<A4qlat,Slip(i,:)=[0,0];
  elseif triC(i,2)>A4plat,Slip(i,:)=[0,0];
  else
   A4sSlip=-cos(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*sin(InterpSV(triC(i,1),triC(i,2)))-sin(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*cos(InterpSV(triC(i,1),triC(i,2)));
   A4dSlip=sin(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*sin(InterpSV(triC(i,1),triC(i,2)))+cos(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*cos(InterpSV(triC(i,1),triC(i,2)));
   Slip(i,:)=[A4sSlip,A4dSlip];
   define(i,1)=[1];
  end
 end 

%2011 3.11 15:15 Ibaraki-oki Mw7.9%
A5plon=140.1;
A5plat=36.35;
A5qlon=140.5;
A5qlat=35.9;

 for i=1:n
  if     define(i,1)==1;
  elseif triC(i,1)<A5plon,Slip(i,:)=[0,0];
  elseif triC(i,1)>A5qlon,Slip(i,:)=[0,0];
  elseif triC(i,2)<A5qlat,Slip(i,:)=[0,0];
  elseif triC(i,2)>A5plat,Slip(i,:)=[0,0];
  else
   A5sSlip=-cos(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*sin(InterpSV(triC(i,1),triC(i,2)))-sin(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*cos(InterpSV(triC(i,1),triC(i,2)));
   A5dSlip=sin(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*sin(InterpSV(triC(i,1),triC(i,2)))+cos(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*cos(InterpSV(triC(i,1),triC(i,2)));
   Slip(i,:)=[A5sSlip,A5dSlip];
   define(i,1)=[1];
  end
 end

%1994 Sanriku-harukaoki Mw7.8%
A6plon=142.475;
A6plat=40.475;
A6qlon=142.85;
A6qlat=40.0;

 for i=1:n
  if     define(i,1)==1;
  elseif triC(i,1)<A6plon,Slip(i,:)=[0,0];
  elseif triC(i,1)>A6qlon,Slip(i,:)=[0,0];
  elseif triC(i,2)<A6qlat,Slip(i,:)=[0,0];
  elseif triC(i,2)>A6plat,Slip(i,:)=[0,0];
  else
   A6sSlip=-cos(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*sin(InterpSV(triC(i,1),triC(i,2)))-sin(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*cos(InterpSV(triC(i,1),triC(i,2)));
   A6dSlip=sin(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*sin(InterpSV(triC(i,1),triC(i,2)))+cos(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*cos(InterpSV(triC(i,1),triC(i,2)));
   Slip(i,:)=[A6sSlip,A6dSlip];
   define(i,1)=[1];
  end
 end

%1968 Aomori-oki Mw8.2; northern asperity of 2 asperities%
A7plon=142.1;
A7plat=41.2;
A7qlon=142.4;
A7qlat=40.75;

 for i=1:n
  if     define(i,1)==1;
  elseif triC(i,1)<A7plon,Slip(i,:)=[0,0];
  elseif triC(i,1)>A7qlon,Slip(i,:)=[0,0];
  elseif triC(i,2)<A7qlat,Slip(i,:)=[0,0];
  elseif triC(i,2)>A7plat,Slip(i,:)=[0,0];
  else
   A7sSlip=-cos(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*sin(InterpSV(triC(i,1),triC(i,2)))-sin(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*cos(InterpSV(triC(i,1),triC(i,2)));
   A7dSlip=sin(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*sin(InterpSV(triC(i,1),triC(i,2)))+cos(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*cos(InterpSV(triC(i,1),triC(i,2)));
   Slip(i,:)=[A7sSlip,A7dSlip];
   define(i,1)=[1];
  end
 end

%1973 Nemuro-oki Mw7.8%
A8plon=145.5;
A8plat=42.95;
A8qlon=145.9;
A8qlat=42.6;

 for i=1:n
  if     define(i,1)==1;
  elseif triC(i,1)<A8plon,Slip(i,:)=[0,0];
  elseif triC(i,1)>A8qlon,Slip(i,:)=[0,0];
  elseif triC(i,2)<A8qlat,Slip(i,:)=[0,0];
  elseif triC(i,2)>A8plat,Slip(i,:)=[0,0];
  else
   A8sSlip=-cos(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*sin(InterpSV(triC(i,1),triC(i,2)))-sin(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*cos(InterpSV(triC(i,1),triC(i,2)));
   A8dSlip=sin(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*sin(InterpSV(triC(i,1),triC(i,2)))+cos(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*cos(InterpSV(triC(i,1),triC(i,2)));
   Slip(i,:)=[A8sSlip,A8dSlip];
   define(i,1)=[1];
  end
 end

%1938 Fukushima-oki%
A9plon=141.65;
A9plat=37.5;
A9qlon=141.95;
A9qlat=37.05;

 for i=1:n
  if     define(i,1)==1;
  elseif triC(i,1)<A9plon,Slip(i,:)=[0,0];
  elseif triC(i,1)>A9qlon,Slip(i,:)=[0,0];
  elseif triC(i,2)<A9qlat,Slip(i,:)=[0,0];
  elseif triC(i,2)>A9plat,Slip(i,:)=[0,0];
  else
   A9sSlip=-cos(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*sin(InterpSV(triC(i,1),triC(i,2)))-sin(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*cos(InterpSV(triC(i,1),triC(i,2)));
   A9dSlip=sin(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*sin(InterpSV(triC(i,1),triC(i,2)))+cos(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*cos(InterpSV(triC(i,1),triC(i,2)));
   Slip(i,:)=[A9sSlip,A9dSlip];
   define(i,1)=[1];
  end
 end

%1978 Miyagi-oki%
A10plon=141.8;
A10plat=38.35;
A10qlon=142.15;
A10qlat=38.025;

 for i=1:n
  if     define(i,1)==1;
  elseif triC(i,1)<A10plon,Slip(i,:)=[0,0];
  elseif triC(i,1)>A10qlon,Slip(i,:)=[0,0];
  elseif triC(i,2)<A10qlat,Slip(i,:)=[0,0];
  elseif triC(i,2)>A10plat,Slip(i,:)=[0,0];
  else
   A10sSlip=-cos(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*sin(InterpSV(triC(i,1),triC(i,2)))-sin(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*cos(InterpSV(triC(i,1),triC(i,2)));
   A10dSlip=sin(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*sin(InterpSV(triC(i,1),triC(i,2)))+cos(sitaS(i)).*InterpVx(triC(i,1),triC(i,2)).*cos(InterpSV(triC(i,1),triC(i,2)));
   Slip(i,:)=[A10sSlip,A10dSlip];
   define(i,1)=[1];
  end
 end

end

