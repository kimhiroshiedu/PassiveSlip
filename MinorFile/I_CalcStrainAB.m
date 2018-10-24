function [B1_Ssn,B1_Sdn,AAA1_Slip]=I_CalcStrainAB(A1_Slip,sSsn,sSdn,dSsn,dSdn)
n=size(A1_Slip,1);
index_out_barri=(A1_Slip(:,3))<1.5;
AA1_Slip=-A1_Slip;
AA1_Slip(index_out_barri,1)=0;
AA1_Slip(index_out_barri,2)=0;
AAA1_Slip=AA1_Slip(1:n,1:2);
ASlip=reshape(AAA1_Slip,2*n,1);
sdn=[sSsn dSsn; sSdn dSdn];
B1_S=sdn*ASlip;
B1_Ssn=B1_S(1:n);
B1_Sdn=B1_S(n+1:2*n);
end
