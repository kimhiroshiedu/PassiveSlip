function [Mw,M0]=MwCalc(Area,MSlip)
n=size(Area,1);
G=40.*(10.^9);
LSlip=zeros(n,1);
for i=1:n
LSlip(i,1)=norm(MSlip(i,1:2));
end
ILSlip=LSlip';
AD=ILSlip*Area;
M0=G.*AD;
LM0=log10(M0);
Mw=(LM0-9.1)./1.5;
end