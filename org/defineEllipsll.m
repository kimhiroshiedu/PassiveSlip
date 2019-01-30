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

