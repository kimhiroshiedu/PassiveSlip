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

