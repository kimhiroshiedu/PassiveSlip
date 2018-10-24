function [Slip,define]=defineSlip2(triC)
%A1=Rectangle asperity_1%
%p=minimum X and max Y, q=max X and minimum Y%
A1px=-100000;
A1py=170000;
A1qx=100000;
A1qy=50000;
A1sSlip=0;
A1dSlip=0.1;
A2px=-70000;
A2py=280000;
A2qx=20000;
A2qy=240000;
A2sSlip=0;
A2dSlip=0.1;
n=length(triC);
define=[A1px,A1py,A1qx,A1qy,A1sSlip,A1dSlip];
Slip=zeros(n,2);
 for i=1:n
  if triC(i,1)<A1px,Slip(i,:)=[0,0];
  elseif triC(i,1)>A1qx,Slip(i,:)=[0,0];
  elseif triC(i,2)<A1qy,Slip(i,:)=[0,0];
  elseif triC(i,2)>A1py,Slip(i,:)=[0,0];
  else Slip(i,:)=[A1sSlip,A1dSlip];
  end
 end
 
 for j=1:n
  if triC(j,1)<A2px
  elseif triC(j,1)>A2qx
  elseif triC(j,2)<A2qy
  elseif triC(j,2)>A2py
  else Slip(j,:)=[A2sSlip,A2dSlip];
  end
 end
end
