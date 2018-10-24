function [I_Ssn,I_Sdn]=I2_CalcStrain(I_Slip,sSsn,sSdn,dSsn,dSdn)
n=size(I_Slip,1);
%
slip=reshape(I_Slip,2*n,1);
sdn=[sSsn sSdn; dSsn dSdn];
%whos
I_S=sdn*slip;
I_Ssn=I_S(1:n);
I_Sdn=I_S(n+1:2*n);
%{
TSlip=I_Slip';

TsSsn=sSsn';
TsSdn=sSdn';
TdSsn=dSsn';
TdSdn=dSdn';
Sum_Ssn=zeros(n,1);
Sum_Sdn=zeros(n,1);

 for j=1:n
  Sum_Ssn(j)=TSlip(1,:)*TsSsn(:,j)+TSlip(2,:)*TdSsn(:,j);
  Sum_Sdn(j)=TSlip(1,:)*TsSdn(:,j)+TSlip(2,:)*TdSdn(:,j);
 end
  I_Ssn=Sum_Ssn;
  I_Sdn=Sum_Sdn;
%}
end

