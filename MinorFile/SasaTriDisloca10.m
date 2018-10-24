function SasaTriDisloca10
 %[triC,tri3,xy,tri]=Sasa_make_trixyz02;
 %save('/home/sasajima/DATA/MAT/triC_2','triC');
 %save('/home/sasajima/DATA/MAT/tri3_2','tri3');
 %save('/home/sasajima/DATA/MAT/tri_2','tri');
 %save('/home/sasajima/DATA/MAT/xy_2','xy');
 load('/home/sasajima/DATA/MAT/tri_2','tri');
 load('/home/sasajima/DATA/MAT/xy_2','xy');
 load('/home/sasajima/DATA/MAT/triC_2','triC');
 load('/home/sasajima/DATA/MAT/tri3_2','tri3');
 %[DistC]=distance_xyz(triC);
 %save('/home/sasajima/DATA/MAT/DistC','DistC');
 load('/home/sasajima/DATA/MAT/DistC','DistC');
 %[sSsn,sSdn,sUxyz]=makeG_s_03(triC,tri3);
 %save('/home/sasajima/DATA/MAT/sSsn_2','sSsn');
 %save('/home/sasajima/DATA/MAT/sSdn_2','sSdn');
 %save('/home/sasajima/DATA/MAT/sUxyz_2','sUxyz');
 %[dSsn,dSdn,dUxyz]=makeG_d_03(triC,tri3);
 %save('/home/sasajima/DATA/MAT/dSsn_2','dSsn');
 %save('/home/sasajima/DATA/MAT/dSdn_2','dSdn');
 %save('/home/sasajima/DATA/MAT/dUxyz_2','dUxyz');
 load('/home/sasajima/DATA/MAT/sSsn_2','sSsn');
 load('/home/sasajima/DATA/MAT/sSdn_2','sSdn');
 load('/home/sasajima/DATA/MAT/dSsn_2','dSsn');
 load('/home/sasajima/DATA/MAT/dSdn_2','dSdn');
 %ss=diag(sSsn);
 %dd=diag(dSdn);
 %figure(5); quiver(triC(:,1),-triC(:,2),ss,dd,1,'r')
 %figure(6); quiver(triC(:,1),-triC(:,2),sSsn(:,3),dSdn(:,3),10,'r')
 %RsSsn=sSsn';
 %RdSdn=dSdn';
 %figure(7); quiver(triC(:,1),-triC(:,2),RsSsn(:,1),RdSdn(:,1),10,'r')
 %figure(107); quiver(triC(:,1),-triC(:,2),RsSsn(:,3),RdSdn(:,3),10,'r')
 %[Slip,define]=defineSlip1(triC);
 %save('/home/sasajima/DATA/MAT/Slip_2','Slip');
 %save('/home/sasajima/DATA/MAT/define_2','define');
 load('/home/sasajima/DATA/MAT/Slip_2','Slip');
 load('/home/sasajima/DATA/MAT/define_2','define');
 %[A_Ssn,A_Sdn]=I_CalcStrain(Slip,sSsn,sSdn,dSsn,dSdn);
 %save('/home/sasajima/DATA/MAT/A_Ssn_2','A_Ssn');
 %save('/home/sasajima/DATA/MAT/A_Sdn_2','A_Sdn');
 load('/home/sasajima/DATA/MAT/A_Ssn_2','A_Ssn');
 load('/home/sasajima/DATA/MAT/A_Sdn_2','A_Sdn');
 %[I_Slip]=I_CalcSlip(A_Ssn,A_Sdn,sSsn,dSdn,Slip);
 %save('/home/sasajima/DATA/MAT/I_Slip_2','I_Slip');
 %load('/home/sasajima/DATA/MAT/I_Slip_2','I_Slip');
 
 %%n=length(Slip);
 %%[I_Slip]=make_I_Slip(Slip);
 %%save('/home/sasajima/DATA/MAT/I_Slip','I_Slip')
 %%load('/home/sasajima/DATA/MAT/I_Slip','I_Slip')
 I_Slip=Slip;
 for i=1:5                                                                                   
 %figure(105); quiver(triC(:,1),triC(:,2),Slip(:,1),Slip(:,2),5,'b')
 %figure(105); quiver(triC(:,1),-triC(:,2),diag(sSsn),diag(dSdn),5,'b')
 %drawnow
 [I_Ssn,I_Sdn]=I2_CalcStrain(I_Slip,sSsn,sSdn,dSsn,dSdn);
 figure(103); quiver(triC(:,1),-triC(:,2),I_Ssn,I_Sdn,3,'r')
 %hold on
 %triplot(tri,xy(:,1),-xy(:,2))
 drawnow

 [I_Slip]=I2_CalcSlip2(DistC,I_Slip,A_Ssn,A_Sdn,I_Ssn,I_Sdn,sSsn,dSdn,Slip,define); 
 figure(100); quiver(triC(:,1),-triC(:,2),I_Slip(:,1),I_Slip(:,2),1,'g')
 %S=sumsSsn;
 %D=sumdSdn;
 %figure(4); quiver(triC(:,1),-triC(:,2),S,D,5,'r')
 %hold on
 %triplot(tri,xy(:,1),-xy(:,2))
 drawnow
 i
 pause(1)
 end
 save('/home/sasajima/DATA/MAT/IF_Ssn_2','I_Ssn');
 save('/home/sasajima/DATA/MAT/IF_Sdn_2','I_Sdn');
 %load('/home/sasajima/DATA/MAT/IF_Ssn_2','I_Ssn');
 %load('/home/sasajima/DATA/MAT/IF_Sdn_2','I_Sdn');
 save('/home/sasajima/DATA/MAT/IF_Slip_2','I_Slip');
 %%save('/home/sasajima/DATA/MAT/Sum_Slip','Sum_Slip');
 %load('/home/sasajima/DATA/MAT/IF_Slip_2','I_Slip');
 %%load('/home/sasajima/DATA/MAT/Sum_Slip','Sum_Slip');
end