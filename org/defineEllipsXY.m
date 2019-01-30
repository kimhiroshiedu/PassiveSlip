function [xyF]=defineEllipsXY

%change value
Pcenter=[,];%[x,y] [m]
a=;%first radius (half radius)[m]
b=;%second radius (half radius)[m]
digArfa=;%amplitude, clockwise [digree]
n=90;
%===================================
radArfa=digArfa./180.*pi;%clockwise [radian]

theta=0-2.*pi./n;
 %make ellips line
 for i=1:n
  theta=theta+2.*pi./n;
  xy1(i,1)=a.*cos(theta);
  xy1(i,2)=b.*cos(theta);

  xy2(i,1)=xy1(i,1).*cos(radArfa)+xy1(i,2).*sin(radArfa);
  xy2(i,2)=-xy1(i,1).*sin(radArfa)+xy1(i,2).*cos(radArfa);

  xyF(i,1)=xy2(i,1)+Pcenter(1,1);
  xyF(i,2)=xy2(i,2)+Pcenter(1,2);
 end

end
