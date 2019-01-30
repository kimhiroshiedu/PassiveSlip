function [AnormaryVecxy]=AnormaryVecCalc(vnorm,SVxy,xyFSlip)

n=length(vnorm);

xyplateV(:,1)=vnorm(:,1).*sin(SVxy(:,1));%x-component of Velocity [m/y]
xyplateV(:,2)=vnorm(:,1).*cos(SVxy(:,1));%y-component of Velocity [m/y]
xyCreep(:,1:2)=xyFSlip(:,1:2)+xyplateV(:,1:2);%[m/y]

 for i=1:n
  XinstSVxy(i,1)=atan2(xyCreep(i,2),xyCreep(i,1));%countour clock wise from X [radian]
  instSVxy(i,1)=2.*pi-XinstSVxy(i,1)+pi.*0.5;%clock wise from Y [radian]
 
  if instSVxy(i,1)>=(2.*pi)
   instSVxy(i,1)=instSVxy(i,1)-2.*pi;
  else
   continue
  end
  
 end

  AnormaryVecxy(:,1)=(SVxy(:,1)-instSVxy(:,1)).*180./pi;%[digree]

end

