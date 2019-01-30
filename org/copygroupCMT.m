function [n]=groupCMT(NARFA,sitaD,sitaS,triC)


 Fid0=fopen('/home/sasajima/Dropbox/yellow/PACdepth201312.txt','r');
 dep_main=textscan(Fid0,'%f %f %f');
 fclose(Fid0);
 dep_main=cell2mat(dep_main);

 E=TriScatteredInterp(dep_main(:,1),dep_main(:,2),dep_main(:,3),'linear');%depth of interplate -km 
 F=TriScatteredInterp(triC(:,1),triC(:,2),sitaS(:,1),'linear');
 G=TriScatteredInterp(triC(:,1),triC(:,2),sitaD(:,1),'linear');
 H=TriScatteredInterp(triC(:,1),triC(:,2),NARFA(:,3),'linear');
 
 Fid1=fopen('/home/sasajima/Dropbox/yellow/NEall.txt','r');
 tmeca=textscan(Fid1,'%f %f %f %f %f %f %f %f %f %f %f %f %f');
 fclose(Fid1);
 meca=cell2mat(tmeca); 

 Fid00=fopen('/home/sasajima/Dropbox/yellow/NEallDate.txt','r');
 
 Fid2=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_thrust97_11Mw9.txt','w');
 Fid3=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEhang97_11Mw9.txt','w');
 Fid4=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_intra_upp97_11Mw9.txt','w');
 Fid5=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_intra_bot97_11Mw9.txt','w');

 Fid7=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_thrust11Mw9_140215.txt','w');
 Fid8=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEhang11Mw9_140215.txt','w');
 Fid9=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_intra_upp11Mw9_140215.txt','w');
 Fid10=fopen('/home_tmp/sasajima/DATA/JPmecaALL/NEPAC_intra_bot11Mw9_140215.txt','w');

 n=size(meca,1);
 tt=0;  
 RD=0;

 for i=1:n  
  i 
  whos
  tline=fgetl(Fid00);
  tt=tt+1;

  RD=RD+1;
  serialTime(1,1)=datenum(tline(1:16),'yyyy/mm/dd,HH:MM');

  d311='2011 03 11 14:46:00';
  t311=datenum(d311,'yyyy mm dd HH:MM:SS');

  lon=meca(i,1);
  lat=meca(i,2);
  focaldep=meca(i,6);
  focalstrike=meca(i,7);
  focaldip=meca(i,8);
  focalSV=meca(i,4)+90;

  thrustdep=-E(lon,lat);
  thruststrike=F(lon,lat);
  thrustdip=G(lon,lat);
  thrustSV=H(lon,lat)*180/pi;

 if serialTime(1,1)<t311;

   if (focaldep>(thrustdep-10))&&(focaldep<(thrustdep+10))&&(focalstrike>(thruststrike-45))&&(focalstrike<(thruststrike+45))&&(focaldip>(thrustdip-12.5))&&(focaldip<(thrustdip+12.5))&&(focalSV>(thrustSV-22.5))&&(focalSV<(thrustSV+22.5));

     fprintf(Fid2,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));

   elseif (focaldep<thrustdep);

     fprintf(Fid3,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));

   elseif (focaldep>thrustdep)&&(focaldep<thrustdep+25);

     fprintf(Fid4,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));
   
   else 

     fprintf(Fid5,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));

   end



 else

   if (focaldep>(thrustdep-10))&&(focaldep<(thrustdep+10))&&(focalstrike>(thruststrike-45))&&(focalstrike<(thruststrike+45))&&(focaldip>(thrustdip-12.5))&&(focaldip<(thrustdip+12.5))&&(focalSV>(thrustSV-22.5))&&(focalSV<(thrustSV+22.5));

     fprintf(Fid7,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));

   elseif (focaldep<thrustdep);

     fprintf(Fid8,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));

   elseif (focaldep>thrustdep)&&(focaldep<thrustdep+25);

     fprintf(Fid9,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));

   else

     fprintf(Fid10,'%10.4f %9.4f %5.0f %5.0f %5.0f %5.0f %5.0f %5.0f %5.2f %4.0f %2.0f %2.0f\n',meca(i,1),meca(i,2),meca(i,3),meca(i,4),meca(i,5),meca(i,6),meca(i,7),meca(i,8),meca(i,9),meca(i,10),meca(i,11),meca(i,12),meca(i,13));

   end

 end

fclose(Fid2);
fclose(Fid3);
fclose(Fid4);
fclose(Fid5);
fclose(Fid7);
fclose(Fid8);
fclose(Fid9);
fclose(Fid10);
fclose(Fid00);

end
