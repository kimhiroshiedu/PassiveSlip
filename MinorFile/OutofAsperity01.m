function [A1_Slip]=OutofAsperity01(A1_Ssn,A1_Sdn,Slip,sSsn,sSdn,dSsn,dSdn)
index_out_asp=(Slip(:,1).^2+Slip(:,2).^2)==0;
%O_sSsn=sSsn(~index_asp,~index_asp);
%O_sSdn=sSdn(~index_asp,~index_asp);
%O_dSsn=dSsn(~index_asp,~index_asp);
%O_dSdn=dSdn(~index_asp,~index_asp);
%O_Ssn=A_Ssn(~index_asp);
%O_Sdn=A_Sdn(~index_asp);
%
k=sum(index_out_asp)
d=[A1_Ssn(index_out_asp);A1_Sdn(index_out_asp)];
G=[sSsn(index_out_asp,index_out_asp),dSsn(index_out_asp,index_out_asp);...
   sSdn(index_out_asp,index_out_asp),dSdn(index_out_asp,index_out_asp)];
m=G\d;
A1_Slip=Slip;
A1_Slip(index_out_asp,1:2)=-[m(1:k) m(k+1:end)];

%{
k=0;
for i=1:n
normSlip=norm(Slip(i,:));
 if normSlip==0
  k=k+1;
  k_sSsn(k,:)=sSsn(i,:);
  k_sSdn(k,:)=sSdn(i,:);
  k_dSsn(k,:)=dSsn(i,:);
  k_dSdn(k,:)=dSdn(i,:);
  O_Ssn(k,1)=A_Ssn(i);
  O_Sdn(k,1)=A_Sdn(i);
 else
 end
end

m=0;

for j=1:n
normSlip=norm(Slip(j,:));
 if normSlip==0
  m=m+1;
  O_sSsn(:,m)=k_sSsn(:,j);
  O_sSdn(:,m)=k_sSdn(:,j);
  O_dSsn(:,m)=k_dSsn(:,j);
  O_dSdn(:,m)=k_dSdn(:,j);
 else
 end
end
%}
end
