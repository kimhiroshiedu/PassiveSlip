function [n]=groupCMT(NARFA,sitaD,sitaS,triC)


 Fid0=fopen('/home/sasajima/Dropbox/yellow/PACdepth201312.txt','r');
 dep_main=textscan(Fid0,'%f %f %f');
 fclose(Fid0);
 dep_main=cell2mat(dep_main);

 E=TriScatteredInterp(dep_main(:,1),dep_main(:,2),dep_main(:,3),'linear');%depth of interplate -km 
 F=TriScatteredInterp(triC(:,1),triC(:,2),sitaS(:,1),'linear');
 G=TriScatteredInterp(triC(:,1),triC(:,2),sitaD(:,1),'linear');
 H=TriScatteredInterp(triC(:,1),triC(:,2),NARFA(:,3),'linear');
 
 Fid1=fopen('/home/sasajima/Dropbox/yellow/NEaftMw9.txt','r');
 tmeca=textscan(Fid1,'%f %f %f %f %f %f %f %f %f %f %f %f %f');
 fclose(Fid1);
 meca=cell2mat(tmeca); 

 %{
 Fid2=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_th_befoMw9.txt','w');
 Fid3=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEhang_befoMw9.txt','w');
 Fid4=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_intra_upp_befoMw9.txt','w');
 Fid5=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_intra_bot_befoMw9.txt','w');
 %}
 %
 Fid2=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_th_aftMw9.txt','w');
 Fid3=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEhang_aftMw9.txt','w');
 Fid4=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_intra_upp_aftMw9.txt','w');
 Fid5=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_intra_bot_aftMw9.txt','w');
 %

 n=size(meca,1);
 tt=0;  
 RD=0;

 for i=1:n  
  i 
  tt=tt+1;

  RD=RD+1;

  lon=meca(i,1);
  lat=meca(i,2);
  focaldep=meca(i,3);
  focalstrike=meca(i,7);
  focaldip=meca(i,8);
  focalSV=meca(i,4)+90;

  thrustdep=-E(lon,lat);
  thruststrike=180+F(lon,lat)*180/pi;
  thrustdip=G(lon,lat)*180/pi;
  thrustSV=H(lon,lat)*180/pi+180;

 if lon<147 

   if (focaldep>(thrustdep-10))&&(focaldep<(thrustdep+10))&&(focalstrike>(thruststrike-45))&&(focalstrike<(thruststrike+45))&&(focaldip>(thrustdip-12.5))&&(focaldip<(thrustdip+12.5))&&(focalSV>(thrustSV-17.5))&&(focalSV<(thrustSV+20));

     fprintf(Fid2,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));

   elseif (focaldep<thrustdep)&&(focaldep<70);

     fprintf(Fid3,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));

   elseif focaldep>thrustdep+25;

     fprintf(Fid5,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));
   
   else 

     fprintf(Fid4,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));

   end

 %
 else

   if (focaldep>(thrustdep-15))&&(focaldep<(thrustdep+15))&&(focalstrike>(thruststrike-45))&&(focalstrike<(thruststrike+45))&&(focaldip>(thrustdip-15))&&(focaldip<(thrustdip+15))&&(focalSV>(thrustSV-20))&&(focalSV<(thrustSV+22.5));

     fprintf(Fid2,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));

   elseif (focaldep<thrustdep)&&(focaldep<70);

     fprintf(Fid3,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));

   elseif focaldep>thrustdep+25;

     fprintf(Fid5,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));

   else

     fprintf(Fid4,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));

   end
 
 end
 %

 end

fclose(Fid2);
fclose(Fid3);
fclose(Fid4);
fclose(Fid5);

end
