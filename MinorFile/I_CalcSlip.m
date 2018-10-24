function [R_Slip]=I_CalcSlip(A_Ssn,A_Sdn,sSsn,dSdn,Slip)
n=length(Slip);
I_Slip=zeros(n,2);
sumsSsn=sum(sSsn);    
sumdSdn=sum(dSdn);

 for i=1:n
 
   R_Slip(i,:)=[-A_Ssn(i)/-sumsSsn(i),-A_Sdn(i)/sumdSdn(i)];%(i,1)=strike slip, (i,2)=dip slip%%Decrease the over slip to max norm of initial slip%  
  
   if abs(Slip(i,1))>0.1
    I_Slip(i,:)=Slip(i,:);
   elseif abs(Slip(i,2))>0.1
    I_Slip(i,:)=Slip(i,:);
   else
    I_Slip(i,:)=I_Slip(i,:);
   end
 end
end

