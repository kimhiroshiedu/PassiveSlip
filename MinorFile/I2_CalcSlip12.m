function [R_Slip,K_Slip,C_Strain,SUM_Slip,sumsSsn,sumdSdn]=I2_CalcSlip12(DistC,R_Slip,SUM_Slip,R_Ssn,R_Sdn,sSsn,sSdn,dSsn,dSdn,Slip,define)
n=length(R_Ssn);
maxslip=[define(5),define(6)];
C_Slip=zeros(n,2);
I_Slip=zeros(2.*n,1);
%atm=atan((abs(maxslip(1,2))/(abs(maxslip(1,1)))))
sumsSsn=sum(sSsn,1);    
sumdSdn=sum(dSdn,1);
ss=diag(sSsn);
dd=diag(dSdn);
whos
%I_Slip(1:n)= R_Ssn./(sumsSsn');
%I_Slip(n+1:2*n)=-R_Sdn./(sumdSdn');%(i,1)=strike slip, (i,2)=dip slip%%Decrease the over slip to max norm of initial slip%  
I_Slip(1:n)=R_Ssn./ss;
I_Slip(n+1:2.*n)=-R_Sdn./dd;
I2_Strain=[sSsn dSsn; sSdn dSdn]*I_Slip;
C_Strain=[I2_Strain(1:n) I2_Strain(n+1:2*n)];
K_Slip=[I_Slip(1:n) I_Slip(n+1:2*n)];
whos
for i=1:n
 
   C_Slip(i,:)=[-R_Ssn(i)./(-C_Strain(i,1)./I_Slip(i)),-R_Sdn(i)./(C_Strain(i,2)./I_Slip(n+i))];%(i,1)=strike slip, (i,2)=dip slip%%Decrease the over slip to max norm of initial slip%  
  
  %for carm down the shaking slip% 
  %{
  if abs(C_Slip(i,1))>abs(I_Slip(i,1))*0.5
   C_Slip(i,1)=I_Slip(i,1).*0.5+C_Slip(i,1).*0.5;
  else
   C_Slip(i,1)=C_Slip(i,1);
  end
  
  if abs(C_Slip(i,2))>abs(I_Slip(i,2))*0.5
   C_Slip(i,2)=I_Slip(i,2).*0.5+C_Slip(i,2).*0.5;
  else
   C_Slip(i,2)=C_Slip(i,2);
  end
  %}
  %{
  %for carm down the shaking slip%
  if (C_Slip(i,1)-I_Slip(i,1))>(abs(C_Slip(i,1))*0.8),C_Slip(i,1)=I_Slip(i,1)*0.5+C_Slip(i,1)*0.5;
   else C_Slip(i,1)=C_Slip(i,1);
  end
  if (C_Slip(i,2)-I_Slip(i,2))>(abs(C_Slip(i,2))*0.8),C_Slip(i,2)=I_Slip(i,2)*0.5+C_Slip(i,2)*0.5;
   else C_Slip(i,2)=C_Slip(i,2);
  end
  %
  %for forbidden the inverse slip to strain direction
  %
  if A_Ssn(i)>0 && C_Slip(i,1)<0, C_Slip(i,1)=-C_Slip(i,1)*0.5;
  elseif A_Ssn(i)<0 && C_Slip(i,1)>0, C_Slip(i,1)=-C_Slip(i,1)*0.5;
  else C_Slip(i,1)=C_Slip(i,1);
  end
  if A_Sdn(i)>0 && C_Slip(i,2)>0, C_Slip(i,2)=-C_Slip(i,2)*0.5;
  elseif A_Sdn(i)<0 && C_Slip(i,2)<0, C_Slip(i,2)=-C_Slip(i,2)*0.5;
  else C_Slip(i,2)=C_Slip(i,2);
  end
  %}
  
  %  atc=atan(abs(C_Slip(i,2))/abs(C_Slip(i,1)));
  % for constraint the direction of slip
  %if abs(atm-atc)>(0.1*pi),C_Slip(i,:)=abs(C_Slip(i,:)).*maxslip./(abs(maxslip)).*0.5+C_Slip(i,:).*0.5;
  % else C_Slip(i,:)=C_Slip(i,:);
  %end
  %
  
  %for forbidden larger slip than one on the asperity
   
   if norm(C_Slip(i,:))>norm(maxslip)
    R_Slip(i,:)=C_Slip(i,:).*(norm(maxslip)./norm(C_Slip(i,:)));
   else
    R_Slip(i,:)=C_Slip(i,:);
   end

   %for fix the slip on the asperity%
   if abs(Slip(i,1))>0.1
      R_Slip(i,:)=[0,0];
   elseif abs(Slip(i,2))>0.1
      R_Slip(i,:)=[0,0];
   else
      R_Slip(i,:)=R_Slip(i,:);
   end
   
  %Sum up the Slip history%
  SUM_Slip(i,:)=SUM_Slip(i,:)+R_Slip(i,:);
      
  %for fix again the slip on the asperity%
  %if abs(Slip(i,1))>0.1
  %  SUM_Slip(i,:)=Slip(i,:);
  %elseif abs(Slip(i,2))>0.1
  %  SUM_Slip(i,:)=Slip(i,:);
  %else
  %    
  %  SUM_Slip(i,:)=SUM_Slip(i,:);
  %end
  
  
  %if abs(I_Slip(i,1))>abs(maxslip(1,1)),I_Slip(i,1)=maxslip(1,1);
  % else I_Slip(i,:)=I_Slip(i,:);
  %end
 
  %if abs(I_Slip(i,2))>abs(maxslip(1,2)),I_Slip(i,2)=maxslip(1,2);
  % else I_Slip(i,:)=I_Slip(i,:);
  %end
 end
 
 %{
  %for the smoothing of slip%
 for i=1:n
  for j=1:n
    if 100000<=DistC(i,j)<200000
        I_Slip(i,:)=I_Slip(j,:).*0.1+I_Slip(i,:).*0.9;
     elseif 50000<=DistC(i,j)<100000
        I_Slip(i,:)=I_Slip(j,:).*0.2+I_Slip(i,:).*0.8;
     elseif 20000<=DistC(i,j)<50000
        I_Slip(i,:)=I_Slip(j,:).*0.3+I_Slip(i,:).*0.7;
     elseif 10000<=DistC(i,j)<20000
        I_Slip(i,:)=I_Slip(j,:).*0.4+I_Slip(i,:).*0.6;
     elseif 5000<=DistC(i,j)<10000
        I_Slip(i,:)=I_Slip(j,:).*0.5+I_Slip(i,:).*0.5;
     elseif 1000<=DistC(i,j)<5000
        I_Slip(i,:)=I_Slip(j,:).*0.5+I_Slip(i,:).*0.5;
     elseif 0<=DistC(i,j)<1000
        I_Slip(i,:)=I_Slip(j,:).*0.5+I_Slip(i,:).*0.5;
     else C_Slip(i,:)=C_Slip(i,:);
    end
  end
  %
   if abs(Slip(i,1))>0.1
    I_Slip(i,:)=Slip(i,:);
   elseif abs(Slip(i,2))>0.1
    I_Slip(i,:)=Slip(i,:);
   else
    I_Slip(i,:)=C_Slip(i,:);
   end
 end

 
  %}
 
end
