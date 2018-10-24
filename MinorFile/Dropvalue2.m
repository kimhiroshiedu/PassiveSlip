function Dropvalue2
 load('/home/sasajima/DATA/MAT/triC_2','triC');
 load('/home/sasajima/DATA/MAT/tri3_2','tri3');
 load('/home/sasajima/DATA/MAT/tri_2','tri');
 load('/home/sasajima/DATA/MAT/xy_2','xy');
 load('/home/sasajima/DATA/MAT/sSsn_2','sSsn');
 load('/home/sasajima/DATA/MAT/sSdn_2','sSdn');
 load('/home/sasajima/DATA/MAT/sUxyz_2','sUxyz');
 load('/home/sasajima/DATA/MAT/dSsn_2','dSsn');
 load('/home/sasajima/DATA/MAT/dSdn_2','dSdn');
 load('/home/sasajima/DATA/MAT/sUxyz_2','sUxyz');
 load('/home/sasajima/DATA/MAT/dUxyz_2','dUxyz');
 load('/home/sasajima/DATA/MAT/Slip_2','Slip');
 load('/home/sasajima/DATA/MAT/define_2','define');
 load('/home/sasajima/DATA/MAT/A_Ssn_2','A_Ssn');
 load('/home/sasajima/DATA/MAT/A_Sdn_2','A_Sdn');
 load('/home/sasajima/DATA/MAT/DistC','DistC');
 %load('/home/sasajima/DATA/MAT/I_Slip_2','I1_Slip');
 load('/home/sasajima/DATA/MAT/IF_Ssn_2','I_Ssn');
 load('/home/sasajima/DATA/MAT/IF_Sdn_2','I_Sdn');
 load('/home/sasajima/DATA/MAT/IF_Slip_2','I_Slip');
 %load('/home/sasajima/DATA/MAT/Sum_Slip','Sum_Slip');
sx=triC(:,1);
sy=triC(:,2);
sz=triC(:,3);
%Slip(:,1)
%Slip(:,2)
%I1_Ssn
%I1_Sdn
%dSsn
%sSdn
%dSsn
%sSdn
%Slip(:,2)
%I1_Slip(:,2)
%I2_Slip(:,2)
%sumSlip=Slip+I1_Slip+I2_Slip;
whos
%DistC
Dx(:,1)=Slip(:,1);
Dy(:,1)=Slip(:,2);
%Dx1(:,1)=Dx(:,1)+I1_Slip(:,1);
%Dy1(:,1)=Dy(:,1)+I1_Slip(:,2);
%Dx2(:,1)=Dx1(:,1)+I2_Slip(:,1);
%Dy2(:,1)=Dy1(:,1)+I2_Slip(:,2);
Dx1(:,1)=dSsn(:,1);
Dy1(:,1)=dSdn(:,1);
Ex=I_Ssn;
Ey=I_Sdn;
%AEx=A_Ssn;
%AEy=A_Sdn;


Data=dSdn(:,1);
dUx=dUxyz(:,2,1);
dUy=dUxyz(:,2,2);
dUz=-dUxyz(:,2,3);
%whos
min_s=0;
max_s=500000;
index_max=Data>max_s;
index_min=Data<min_s;
index_nan=isnan(Data);
color_index=jet(128);
RRR=fix((Data-min_s)./((max_s-min_s)./128)+1);
RRR(index_max)=128;
RRR(index_min)=1;
RRR(index_nan)=1;
RR(:,1)=color_index(RRR,1);
RR(:,2)=color_index(RRR,2);
RR(:,3)=color_index(RRR,3);
%figure(1500)
%scatter(sx,sy,20,RR,'filled')
%scatter(sx,sy,RSyz,'filled')
%fclose(Fid2);
%figure(101); quiver(sx,sy,Ex,Ey,5,'r')
figure(1); quiver(sx,sy,Dx1,Dy1,5,'g')
hold on
 triplot(tri,xy(:,1),xy(:,2))
%figure(200); quiver(sx,sy,Dx,Dy,'b')
%figure(3); quiver3(sx,sy,sz,dUx,dUy,dUz,'b')
%figure(4); quiver(sx,sy,dUx,dUy,50,'r'); grid on;
end
