function [I_Slip]=make_I_Slip(Slip)
slipX=-50;
slipY=0;
n=length(Slip);
I_Slip=zeros(n,2);
 for i=1:n
   I_Slip(i,:)=[slipX,slipY]; 
 end
end