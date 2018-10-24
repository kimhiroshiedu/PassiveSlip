function [R_Ssn,R_Sdn,SUM_Ssn,SUM_Sdn]=I2_CalcStrain12(R_Slip,SUM_Ssn,SUM_Sdn,sSsn,sSdn,dSsn,dSdn)
n=size(R_Slip,1);
%
slip=reshape(R_Slip,2*n,1);
sdn=[sSsn dSsn; sSdn dSdn];
%whos
R_S=sdn*slip;
R_Ssn=R_S(1:n);
R_Sdn=R_S(n+1:2*n);
SUM_Ssn=SUM_Ssn+R_Ssn;
SUM_Sdn=SUM_Sdn+R_Sdn;
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

