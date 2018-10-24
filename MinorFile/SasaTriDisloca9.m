function SasaTriDisloca9

 [triC,tri3]=Sasa_make_trill;
 save('/home/sasajima/DATA/MAT/triC_Q','triC');
 save('/home/sasajima/DATA/MAT/tri3_Q','tri3');
 %load('/home/sasajima/DATA/MAT/triC_Q','triC');
 %load('/home/sasajima/DATA/MAT/tri3_Q','tri3');

 [trixyzC,trixyz3]=Sasa_ll2xyz;
 save('/home/sasajima/DATA/MAT/triC_Q','triC');
 save('/home/sasajima/DATA/MAT/tri3_Q','tri3');
 %load('/home/sasajima/DATA/MAT/triC_Q','triC');
 %load('/home/sasajima/DATA/MAT/tri3_Q','tri3');

 [sitaS,sitaD]=strike_dip(trixyzC,trixyz3);
 save('/home/sasajima/DATA/MAT/sitaS','sitaS');
 save('/home/sasajima/DATA/MAT/sitaD','sitaD');
 %load('/home/sasajima/DATA/MAT/sitaS','sitaS');
 %load('/home/sasajima/DATA/MAT/sitaD','sitaD');

 [sR2Sxz,sR2Syz]=makeG_s_Q(trixyzC,trixyz3,sitaS,sitaD);
 save('/home/sasajima/DATA/MAT/sRSxz_Q','sR2Sxz');
 save('/home/sasajima/DATA/MAT/sRSyz_Q','sR2Syz');

 [dR2Sxz,dR2Syz]=makeG_d_Q(trixyzC,trixyz3,sitaS,sitaD);
 save('/home/sasajima/DATA/MAT/dRSxz_Q','dR2Sxz');
 save('/home/sasajima/DATA/MAT/dRSyz_Q','dR2Syz');
 %load('/home/sasajima/DATA/MAT/sRSxz_Q','sR2Sxz');
 %load('/home/sasajima/DATA/MAT/sRSyz_Q','sR2Syz');
 %load('/home/sasajima/DATA/MAT/dRSxz_Q','dR2Sxz');
 %load('/home/sasajima/DATA/MAT/dRSyz_Q','dR2Syz');

 [Slip,define]=defineSlipQ(triC,tri3,trixyzC,trixyz3);
 save('/home/sasajima/DATA/MAT/Slip_Q1','Slip');
 save('/home/sasajima/DATA/MAT/define_Q1','define');
 %load('/home/sasajima/DATA/MAT/Slip_Q1','Slip');
 %load('/home/sasajima/DATA/MAT/define_Q1','define');

 [I1_R2Sxz,I1_R2Syz]=I_CalcStrain(Slip,sR2Sxz,sR2Syz,dR2Sxz,dR2Syz);
 save('/home/sasajima/DATA/MAT/I1_R2Sxz_Q1','I1_R2Sxz');
 save('/home/sasajima/DATA/MAT/I1_R2Syz_Q1','I1_R2Syz');
 %load('/home/sasajima/DATA/MAT/I1_R2Sxz_Q1','I1_R2Sxz');
 %load('/home/sasajima/DATA/MAT/I1_R2Syz_Q1','I1_R2Syz');

 [I1_Slip]=I_CalcSlip(I1_R2Sxz,I1_R2Syz,sR2Sxz,sR2Syz,dR2Sxz,dR2Syz,Slip,define);
 save('/home/sasajima/DATA/MAT/I1_Slip_Q1','I1_Slip');
 %load('/home/sasajima/DATA/MAT/I1_Slip_Q1','I1_Slip');

 [I2_R2Sxz,I2_R2Syz]=I2_CalcStrain(I1_Slip,sR2Sxz,sR2Syz,dR2Sxz,dR2Syz);
 save('/home/sasajima/DATA/MAT/I2_R2Sxz_Q1','I2_R2Sxz');
 save('/home/sasajima/DATA/MAT/I2_R2Syz_Q1','I2_R2Syz');
 %load('/home/sasajima/DATA/MAT/I2_R2Sxz_Q1','I2_R2Sxz');
 %load('/home/sasajima/DATA/MAT/I2_R2Syz_Q1','I2_R2Syz');

 [I2_Slip]=I2_CalcSlip(I1_R2Sxz,I1_R2Syz,I2_R2Sxz,I2_R2Syz,sR2Sxz,sR2Syz,dR2Sxz,dR2Syz,Slip,define);
 save('/home/sasajima/DATA/MAT/I2_Slip_Q1','I2_Slip');
 %load('/home/sasajima/DATA/MAT/I2_Slip_Q1','I2_Slip');

end
