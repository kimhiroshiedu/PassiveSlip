function [Area]=AreaCalc(tri3)
n=length(tri3);
ab=zeros(n,3);
bc=zeros(n,3);
ac=zeros(n,3);
Nab=zeros(n,1);
Nac=zeros(n,1);
Nbc=zeros(n,1);
s=zeros(n,1);
Area=zeros(n,1);
for i=1:n
ab(i,1:3)=[tri3(i,1,2)-tri3(i,1,1),tri3(i,2,2)-tri3(i,2,1),tri3(i,3,2)-tri3(i,3,1)];
ac(i,1:3)=[tri3(i,1,3)-tri3(i,1,1),tri3(i,2,3)-tri3(i,2,1),tri3(i,3,3)-tri3(i,3,1)];
bc(i,1:3)=[tri3(i,1,3)-tri3(i,1,2),tri3(i,2,3)-tri3(i,2,2),tri3(i,3,3)-tri3(i,3,2)];
Nab(i,1)=norm(ab(i,1:3));
Nac(i,1)=norm(ac(i,1:3));
Nbc(i,1)=norm(bc(i,1:3));
s(i,1)=(Nab(i)+Nac(i)+Nbc(i))./2;
Area(i,1)=(s(i).*(s(i)-Nab(i)).*(s(i)-Nac(i)).*(s(i)-Nbc(i))).^0.5;
end
end
