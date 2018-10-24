function [FO_Ssn,FO_Sdn]=OutofAsperity11(F_Ssn,F_Sdn,Slip)
n=length(Slip);
k=0;
for i=1:n
normSlip=norm(Slip(i,:));
 if normSlip==0
  k=k+1;
 
  FO_Ssn(k,1)=F_Ssn(i);
  FO_Sdn(k,1)=F_Sdn(i);
 else
 end
end
end

