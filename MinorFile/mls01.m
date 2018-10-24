function [F_Slip]=mls01(O_Ssn,O_Sdn,O_sSsn,O_sSdn,O_dSsn,O_dSdn,Slip)
index_asp=(Slip(:,1).^2+Slip(:,2).^2)~=0;
k=size(O_Ssn,1);
d=[O_Ssn;O_Sdn];
G=[O_sSsn,O_sSdn;O_dSsn,O_dSdn];
m=G\d;
F_Slip=Slip;
F_Slip(~index_asp,1:2)=-[m(1:k) m(k+1:end)];

%{
n=length(Slip);
k=length(O_Ssn);
d(1:k,1)=O_Ssn;
d(k+1:2*k,1)=O_Sdn;
G=[O_sSsn,O_sSdn;O_dSsn,O_dSdn];
whos
m=G\d;
j=0;
F_sSlip=zeros(n,1);
F_dSlip=zeros(n,1);
for i=1:n;
normSlip=norm(Slip(i,:));

 if normSlip==0
  j=j+1;
  F_sSlip(i)=-m(j);
  F_dSlip(i)=-m(j+k);
 else
  F_sSlip(i)=Slip(i,1);
  F_dSlip(i)=Slip(i,2);
 end
end
F_Slip(:,1)=F_sSlip(:);
F_Slip(:,2)=F_dSlip(:);
%}