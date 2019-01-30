function DropValueQ
 load('/home_tmp/sasajima/DATA/PassiveSlip/Tohoku/tri_Tohoku','triC','tri3','tri','sll');
 load('/home_tmp/sasajima/DATA/PassiveSlip/Tohoku/trixyz_Tohoku','trixyzC','trixyz3','sxyz');
 load('/home_tmp/sasajima/DATA/PassiveSlip/Tohoku/sita_Tohoku','sitaS','sitaD','normVec');

 %[sSsn,sSdn,dSsn,dSdn]=loatMAT(trixyzC);

 %load('/home_tmp/sasajima/DATA/MAT/sSsn_Boso','sSsn');
 %load('/home_tmp/sasajima/DATA/MAT/sSdn_Q07','sSdn');
 %load('/home_tmp/sasajima/DATA/MAT/dSsn_Q07','dSsn');
 %load('/home_tmp/sasajima/DATA/MAT/dSdn_Boso','dSdn');

 %load('/home_tmp/sasajima/DATA/MAT/sxyz_Q07','sUxyz');
 %load('/home_tmp/sasajima/DATA/MAT/dUxyz_Q07','dUxyz');

sx=triC(:,1);
sy=triC(:,2);
sz=triC(:,3);
%Slip(:,1)
%Slip(:,2)
%I1_R2Syz
%I1_R2Syz
%dR2Syz
%sR2Syz
%dR2Sxz
%sR2Sxz
%Slip(:,2)
%I1_Slip(:,2)
%I2_Slip(:,2)
%sumSlip=Slip+I1_Slip+I2_Slip;
whos
%A_Ssn(:);
%dSdn(101,:)
%Value=triC(:,3);
Value=sitaS(:,1);

min_s=pi./2
max_s=pi
index_max=Value>max_s;
index_min=Value<min_s;
index_nan=isnan(Value);
color_index=jet(128);
RRR=fix((Value-min_s)./((max_s-min_s)./128)+1);
RRR(index_max)=128;
RRR(index_min)=1;
RRR(index_nan)=1;
RR(:,1)=color_index(RRR,1);
RR(:,2)=color_index(RRR,2);
RR(:,3)=color_index(RRR,3);
figure(1);
scatter(sx,sy,20,RR,'filled')
%hold on
%figure(2)
%triplot(tri,1000.*sll(:,1),1000.*sll(:,2));
%fclose(Fid2);
%figure(15); quiver3(sx,sy,sz,dUxyz(:,1,1),dUxyz(:,1,2),dUxyz(:,1,3),'r')
%figure(16); quiver(sx,sy,normVec(:,1),-normVec(:,2),'r');
%grid on;
end

