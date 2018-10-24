function[sitaS,sitaD,normVec]=cmtstrike_dip(trixyzC,trixyz3)
%whos
%Define strike and dip of fault%
n=size(trixyzC,1);
ca=zeros(n,3);
cb=zeros(n,3);
InormVec=zeros(n,3);
normVec=zeros(n,3);
strikeVec=zeros(n,3);
dipVec=zeros(n,3);
RsitaS=zeros(n,1);
Rdip=zeros(n,1);
sitaS=zeros(n,1);
sitaD=zeros(n,1);

for i=1:n
 ca(i,:)=[trixyz3(i,1,1)-trixyz3(i,1,3),trixyz3(i,2,1)-trixyz3(i,2,3),trixyz3(i,3,1)-trixyz3(i,3,3)];
 cb(i,:)=[trixyz3(i,1,2)-trixyz3(i,1,3),trixyz3(i,2,2)-trixyz3(i,2,3),trixyz3(i,3,2)-trixyz3(i,3,3)];
 InormVec(i,:)              = cross(ca(i,:),cb(i,:),2);% direct to deeper %
 normVec(i,:)              = InormVec(i,:)./(sqrt((InormVec(i,1)).^2+(InormVec(i,2)).^2+(InormVec(i,3).^2)));
   if (normVec(i,3) < 0) % Enforce clockwise circulation
    normVec(i,:)               = -normVec(i,:);
   % [x(2) x(3)]               = swap(x(2), x(3));%???
   % [y(2) y(3)]               = swap(y(2), y(3));%???
   % [z(2) z(3)]               = swap(z(2), z(3));%???
   end
 strikeVec(i,:)              = [-sin(atan2(normVec(i,2),normVec(i,1))),cos(atan2(normVec(i,2),normVec(i,1))),0];%direct to left hand of who toward dip direction
 dipVec(i,:)                 = cross(normVec(i,:), strikeVec(i,:),2);

   %if normVec(1)==0 and normVec(2)==0
   % strikeVec=[1 0 0];
   % dipVec   =[0 1 0];
   % normVec  =[0 0 1];
   %end

  RsitaS(i)=atan2(strikeVec(i,2),strikeVec(i,1)); %from x-axis to y-axis rotation [rad]%
  %RsitaSY=2.5*pi-RsitaS; %from y-axis clockwise[rad]%
  %sitaSY=RsitaSY*180/pi;
  sitaS(i)=RsitaS(i)-0.5*pi; %from Y-axis(=North) to clockwise [rad]

  Rdip(i)=asin(dipVec(i,3));
  %dip=Rdip*180/pi;
  sitaD(i)=Rdip(i);%dip angle[rad]
end
end
