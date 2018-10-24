function [I_Slip]=I2_CalcSlip(DistC,I_Slip,A_Ssn,A_Sdn,I_Ssn,I_Sdn,sSsn,dSdn,Slip,define)
n=length(I_Ssn);
maxslip=[define(5),define(6)];
C_Slip=zeros(n,2);
atm=atan((abs(maxslip(1,2))/(abs(maxslip(1,1)))));

 for i=1:n
   C_Slip(i,:)=[-I_Ssn(i)/sSsn(i,i),-I_Sdn(i)/dSdn(i,i)];%(i,1)=strike slip, (i,2)=dip slip%%Decrease the over slip to max norm of initial slip%  
  
  if (C_Slip(i,1)-I_Slip(i,1))>(abs(C_Slip(i,1))*0.5),C_Slip(i,1)=I_Slip(i,1)*0.5+C_Slip(i,1)*0.5;
   else C_Slip(i,1)=C_Slip(i,1);
  end
  
  if (C_Slip(i,2)-I_Slip(i,2))>(abs(C_Slip(i,2))*0.8),C_Slip(i,2)=I_Slip(i,2)*0.5+C_Slip(i,2)*0.5;
   else C_Slip(i,2)=C_Slip(i,2);
  end
   
  if A_Ssn(i)>0 && C_Slip(i,1)<0, C_Slip(i,1)=-C_Slip(i,1)*0.5;
  elseif A_Ssn(i)<0 && C_Slip(i,1)>0, C_Slip(i,1)=-C_Slip(i,1)*0.5;
  else C_Slip(i,1)=C_Slip(i,1);
  end
  
  if A_Sdn(i)>0 && C_Slip(i,2)>0, C_Slip(i,2)=-C_Slip(i,2)*0.5;
  elseif A_Sdn(i)<0 && C_Slip(i,2)<0, C_Slip(i,2)=-C_Slip(i,2)*0.5;
  else C_Slip(i,2)=C_Slip(i,2);
  end
  
    atc=atan(abs(C_Slip(i,2))/abs(C_Slip(i,1)));
  
  if abs(atm-atc)>(0.15*pi),C_Slip(i,1)=C_Slip(i,1)/abs(C_Slip(i,1))*(abs(C_Slip(i,2))*tan(0.15*pi));
   else C_Slip(i,1)=C_Slip(i,1);
  end
  
  if norm(C_Slip(i,:))>norm(maxslip),C_Slip(i,:)=C_Slip(i,:)*(norm(maxslip)/norm(C_Slip(i,:)));
   else C_Slip(i,:)=C_Slip(i,:);
  end
  

      
  if abs(Slip(i,1))>0.1,C_Slip(i,:)=Slip(i,:);
   elseif abs(Slip(i,2))>0.1,C_Slip(i,:)=Slip(i,:);
   else C_Slip(i,:)=C_Slip(i,:);
  end
  
  
  
  
  
  
  %if abs(I_Slip(i,1))>abs(maxslip(1,1)),I_Slip(i,1)=maxslip(1,1);
  % else I_Slip(i,:)=I_Slip(i,:);
  %end
  
  %if abs(I_Slip(i,2))>abs(maxslip(1,2)),I_Slip(i,2)=maxslip(1,2);
  % else I_Slip(i,:)=I_Slip(i,:);
  %end
  
 end
 I_Slip=C_Slip;
 
end

