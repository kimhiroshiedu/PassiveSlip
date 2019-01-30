function [A1_Ssn,A1_Sdn]=I_CalcStrain(Slip,sSsn,sSdn,dSsn,dSdn)
n=size(Slip,1);
SSlip=Slip(1:n,1:2);
slip=reshape(SSlip,2*n,1);
sdn=[sSsn dSsn; sSdn dSdn];

whos

A1_S=sdn*slip;
A1_Ssn=A1_S(1:n);
A1_Sdn=A1_S(n+1:2*n);
end
