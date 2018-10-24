function [Slip,define]=defineSlipAB(triC)
%A1=Rectangle asperity_1%
%B1=Rectangle barrier_1%
%p=minimum X and max Y, q=max X and minimum Y%
A1px=-70000;
A1py=280000;
A1qx=20000;
A1qy=240000;
B1px=-100000;
B1py=170000;
B1qx=100000;
B1qy=50000;
A1sSlip=0;
A1dSlip=-10;
n=length(triC);
define=[A1px,A1py,A1qx,A1qy,A1sSlip,A1dSlip];
Slip=zeros(n,3);
 %Asperity1%
 for i=1:n
  if triC(i,1)<A1px
     Slip(i,:)=[0,0,0];
   elseif triC(i,1)>A1qx
     Slip(i,:)=[0,0,0];
   elseif triC(i,2)<A1qy
     Slip(i,:)=[0,0,0];
   elseif triC(i,2)>A1py
     Slip(i,:)=[0,0,0];
   else 
     Slip(i,:)=[A1sSlip,A1dSlip,1];
  end
 end 
 %Barrier1%
 for i=1:n
    if triC(i,1)<B1px
     elseif triC(i,1)>B1qx
     elseif triC(i,2)<B1qy
     elseif triC(i,2)>B1py
     else
     Slip(i,:)=[0,0,2];
    end
 end
end

