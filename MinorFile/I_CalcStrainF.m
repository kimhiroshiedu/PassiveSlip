function [F_Ssn,F_Sdn]=I_CalcStrainF(F_Slip,sSsn,sSdn,dSsn,dSdn)

n=size(F_Slip,1);
FF_Slip=F_Slip(1:n,1:2);
slip=reshape(FF_Slip,2*n,1);
sdn=[sSsn dSsn; sSdn dSdn];
F_S=sdn*slip;
F_Ssn=F_S(1:n);
F_Sdn=F_S(n+1:2*n);
end
