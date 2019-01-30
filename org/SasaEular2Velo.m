function [vnorm,v,NARFA]=SasaEular2Velo(ll)

Fid=fopen('/home_tmp/sasajima/DATA/Pole_data/PAC-OKH_MORVEL56.dat','r');
Epole=textscan(Fid,'%f %f %f');
fclose(Fid);
Epole=cell2mat(Epole);

%Fid2=fopen('/home_tmp/sasajima/DATA/boso/IA-KA_vll.txt','w');

R(1,1)=6371000;%[m]
%a=6378137.0;
%b=6356752.3;
%aabb=(a.^2)./(b.^2);

diglonP(1,1)=Epole(1,1);%longutitude of Eular pole [digree]%
diglatP(1,1)=Epole(1,2);%latitude of Eular pole [digree]%
digwMa(1,1)=Epole(1,3);%Rotation velocity [digree/Ma]%
radlonP=diglonP./180.*pi;
radlatP(1,1)=diglatP(1,1)./180.*pi;
radwMa=digwMa./180.*pi;

radlonX(:,1)=ll(:,1)./180.*pi;%lon of observed point[radian]
radlatX(:,1)=ll(:,2)./180.*pi;%lat of observed point[radian]
diglonX(:,1)=ll(:,2);%lon of observed poiint [digree]

NARFA(:,1)=radlonX(:,1);
NARFA(:,2)=radlatX(:,1);

CE(1,1)=0.5.*pi-radlatP(1,1);%[radian],scolar
CB(:,1)=0.5.*pi-radlatX(:,1);%[radian],vector
RMD(:,1)=radlonP-radlonX(:,1);%[radian],vector

cCE(1,1)=cos(CE);
sCE(1,1)=sin(CE);

n=length(ll);
v=zeros(n,2);
vnorm=zeros(n,1);

 for i=1:n;
  cCB(1,1)=cos(CB(i,1));
  sCB(1,1)=sin(CB(i,1));
  cRMD(1,1)=cos(RMD(i,1));

  cDLT(1,1)=cCB(1,1).*cCE(1,1)+sCB(1,1).*sCE(1,1).*cRMD(1,1);
  DLT(1,1)=acos(cDLT(1,1));

  sDLT(1,1)=sin(DLT(1,1));
  L(1,1)=sDLT(1,1).*R(1,1);%distance between obs point and Eular Vector [m]

  vmMa(1,1)=L(1,1).*radwMa(1,1);%velocity[m/Ma]
  vnorm(i,1)=vmMa(1,1)./10000;%velocity[cm/year]

  cARFA(1,1)=(cCE(1,1)-(cCB(1,1).*cDLT(1,1)))./(sCB(1,1).*sDLT(1,1));
  ARFA(1,1)=acos(cARFA(1,1));%[radian]

  dltdiglon(1,1)=diglonP(1,1)-diglonX(i,1);%[digree]

   if dltdiglon>0;
    NARFA(i,3)=0.5.*pi-ARFA(1,1);%countour? clock wise from North [radian]
   else
    NARFA(i,3)=0.5.*pi+ARFA(1,1);%countour? clock wise from North[radian]
   end

  v(i,1)=vnorm(i,1).*cos(NARFA(i,3));%vlat;N=plus[cm/year]
  v(i,2)=vnorm(i,1).*sin(NARFA(i,3));%vlon;E=plus[cm/year]

  %digNARFA(i,1)=NARFA(i,3).*180./pi;%clock wise from North[digree]
   
  %fprintf(Fid2,'%15.8f %15.8f %15.8f %15.8f\n',v(i,1),v(i,2),vnorm(i,1),digNARFA(i,1));

 end

 %fclose(Fid2);
end

