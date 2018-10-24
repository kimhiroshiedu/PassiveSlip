function [I_Slip]=I2_CalcSlip2(DistC,I_Slip,A_Ssn,A_Sdn,I_Ssn,I_Sdn,sSsn,dSdn,Slip,define)
n=length(I_Ssn);
maxslip=[define(5),define(6)];
C_Slip=zeros(n,2);
%atm=atan((abs(maxslip(1,2))/(abs(maxslip(1,1)))));
%sumsSsn=sum(sSsn,2);    
%sumdSdn=sum(dSdn,2);
whos
 for i=1:n
 
   %C_Slip(i,:)=[-I_Ssn(i)/-sumsSsn(i,1),-I_Sdn(i)/sumdSdn(i,1)];%(i,1)=strike slip, (i,2)=dip slip%%Decrease the over slip to max norm of initial slip%  
  
   C_Slip(i,:)=[1.*I_Ssn(i)./sSsn(i,i),-1.*I_Sdn(i)./dSdn(i,i)];
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
  %
  %
  %for carm down the shaking slip%
  if (C_Slip(i,1)-I_Slip(i,1))>(abs(C_Slip(i,1))*0.8),C_Slip(i,1)=I_Slip(i,1)*0.5+C_Slip(i,1)*0.5;
   else C_Slip(i,1)=C_Slip(i,1);
  end
  if (C_Slip(i,2)-I_Slip(i,2))>(abs(C_Slip(i,2))*0.8),C_Slip(i,2)=I_Slip(i,2)*0.5+C_Slip(i,2)*0.5;
   else C_Slip(i,2)=C_Slip(i,2);
  end
  %}
 
  
  
  %}
  
    %atc=atan(abs(C_Slip(i,2))/abs(C_Slip(i,1)));
  % for constraint the direction of slip
  %if abs(atm-atc)>(0.1*pi),C_Slip(i,:)=abs(C_Slip(i,:)).*maxslip./(abs(maxslip)).*0.5+C_Slip(i,:).*0.5;
  % else C_Slip(i,:)=C_Slip(i,:);
  %end
  %
  %for forbidden larger slip than one on the asperity
  %if norm(C_Slip(i,:))>norm(maxslip),C_Slip(i,:)=C_Slip(i,:)*(norm(maxslip)/norm(C_Slip(i,:)));
  % else C_Slip(i,:)=C_Slip(i,:);
  %end
  %}
  %for fix the slip on the asperity%
  %if abs(Slip(i,1))>0.1,C_Slip(i,:)=Slip(i,:);
  % elseif abs(Slip(i,2))>0.1,C_Slip(i,:)=Slip(i,:);
  % else C_Slip(i,:)=C_Slip(i,:);
  %end
  

  %Sum up the Slip history%
  I_Slip(i,:)=I_Slip(i,:)+C_Slip(i,:);
  
   %for forbidden the inverse slip to strain direction
  if I_Slip(i)*maxslip <0
    I_Slip(i,:)=[0,0];
  else
    I_Slip(i,:)=I_Slip(i,:);
  end
  %{
   %for the smoothing of slip%
 %for i=1:n
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
     else I_Slip(i,:)=I_Slip(i,:);
    end
  end
  %}
  
  
  %for forbidden larger slip than one on the asperity
  if norm(I_Slip(i,:))>norm(maxslip)
    I_Slip(i,:)=I_Slip(i,:).*(norm(maxslip)/norm(I_Slip(i,:)));
  else
    I_Slip(i,:)=I_Slip(i,:);
  end
  
  %for fix again the slip on the asperity%
  if abs(Slip(i,1))>0.1
    I_Slip(i,:)=Slip(i,:);
  elseif abs(Slip(i,2))>0.1
    I_Slip(i,:)=Slip(i,:);
  else
    I_Slip(i,:)=I_Slip(i,:);
  end
  
  
  %if abs(I_Slip(i,1))>abs(maxslip(1,1)),I_Slip(i,1)=maxslip(1,1);
  % else I_Slip(i,:)=I_Slip(i,:);
  %end
  
  %if abs(I_Slip(i,2))>abs(maxslip(1,2)),I_Slip(i,2)=maxslip(1,2);
  % else I_Slip(i,:)=I_Slip(i,:);
  %end
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
 

