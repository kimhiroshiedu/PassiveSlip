 function [DistC]=distance_xyz(triC)
n=length(triC);
DistC=zeros(n,n);

 for i=1:n
  for j=1:n
    DistC(i,j)=((triC(i,1)-triC(j,1)).^2+(triC(i,2)-triC(j,2)).^2+(triC(i,3)-triC(j,3)).^2).^0.5;
  end
 end
end
