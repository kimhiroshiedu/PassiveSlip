function Expand_Allsample(folder,burnin)
% burnin: enter for percent scale
folder  = fullfile(pwd,folder);
[input] = read_chafile(folder);
[tcha]  = cal_avestdbin(input,burnin);

% save
outfile=[folder,'/tcha.mat'];
save(outfile,'tcha','-v7.3');
end

%% read all cha mat file
function [input] = read_chafile(folder)
ext  = 'cha_test*.mat';
file = dir([folder,'/',ext]);
[nit,~] = size(file);
for ii = 1:nit
  input(ii).fname = fullfile(file(ii).folder,file(ii).name);
end
end

%% calculate average, std, bin
function [tcha] = cal_avestdbin(input,burnin)
% load('./result/test_06/cha_test.mat'); % test
nit = size(input,2);
% sfactor = 2^8;
sfactor = 2^16;
mpbin = [-1e-7:1e-10:1e-7];
mcbin = [-1:1e-3:1];
mibin = [-1e-7:1e-10:1e-7];
mabin = [0:0.1:100];
smpint = 50; % sampling interval
burnin = floor(burnin*nit/100)+1;
acctotal = 0;
for ii = 1:nit
  load(input(ii).fname);
  accflag = isfield(cha,'ajr');
  if ii==1
    nch = size(cha.mpcompress.smpmp,2);
    npol   = size(cha.mpcompress.npol,2);
    nflt   = size(cha.mccompress.nflt,2);
    nine   = size(cha.micompress.nine,2);
    nasp   = size(cha.macompress.nasp,2);
    naspid = size(cha.maidcompress.smpmaid,1);
    sumpol   = zeros(npol,1);
    sumflt   = zeros(nflt,1);
    sumine   = zeros(nine,1);
    sumasp   = zeros(nasp,1);
    sumaspid = zeros(naspid,1);
    sumpolpair   = zeros(npol,npol);
    sumfltpair   = zeros(nflt,nflt);
    suminepair   = zeros(nine,nine);
    sumasppair   = zeros(nasp,nasp);
    sumaspidpair = zeros(naspid,naspid);
    ndatapol   = 0;
    ndataflt   = 0;
    ndataine   = 0;
    ndataasp   = 0;
    ndataaspid = 0;
    mphist = zeros(npol,size(mpbin,2)-1);
    mchist = zeros(nflt,size(mcbin,2)-1);
    mihist = zeros(nine,size(mibin,2)-1);
    mahist = zeros(nasp,size(mabin,2)-1);
    smpid = [1:smpint:nch];
    medpol   = zeros(npol,1);
    medflt   = zeros(nflt,1);
    medine   = zeros(nine,1);
    medasp   = zeros(nasp,1);
    medaspid = zeros(naspid,1);
    asmppol   = [];
    asmpflt   = [];
    asmpine   = [];
    asmpasp   = [];
    asmpaspid = [];
  end
  if ii>burnin
    %
    for np = 1:npol
      infid = cha.mpcompress.npol(np).mpscale==inf;
      if ~infid
        smppol(np,:) = (double(cha.mpcompress.smpmp(np,:))+(sfactor/2))./((sfactor-1).*cha.mpcompress.npol(np).mpscale)+cha.mpcompress.npol(np).mpmin;
      else
        smppol(np,:) = ones(1,nch).*cha.mpcompress.npol(np).mpmax;
      end
      mphist(np,:) = mphist(np,:)+histcounts(smppol(np,:),mpbin);
    end
    %
    for nf = 1:nflt
      infid = cha.mccompress.nflt(nf).mcscale==inf;
      if ~infid
        smpflt(nf,:) = (double(cha.mccompress.smpmc(nf,:))+(sfactor/2))./((sfactor-1).*cha.mccompress.nflt(nf).mcscale)+cha.mccompress.nflt(nf).mcmin;
      else
        smpflt(nf,:) = ones(1,nch).*cha.mccompress.nflt(nf).mcmax;
      end
      mchist(nf,:) = mchist(nf,:)+histcounts(smpflt(nf,:),mcbin);
    end
    %
    for ni = 1:nine
      infid = cha.micompress.nine(ni).miscale==inf;
      if ~infid
        smpine(ni,:) = (double(cha.micompress.smpmi(ni,:))+(sfactor/2))./((sfactor-1).*cha.micompress.nine(ni).miscale)+cha.micompress.nine(ni).mimin;
      else
        smpine(ni,:) = ones(1,nch).*cha.micompress.nine(ni).mimax;
      end
      mihist(ni,:) = mihist(ni,:)+histcounts(smpine(ni,:),mibin);
    end
    %
    for na = 1:nasp
      infid = cha.macompress.nasp(na).mascale==inf;
      if ~infid
        smpasp(na,:) = (double(cha.macompress.smpma(na,:))+(sfactor/2))./((sfactor-1).*cha.macompress.nasp(na).mascale)+cha.macompress.nasp(na).mamin;
      else
        smpasp(na,:) = ones(1,nch).*cha.macompress.nasp(na).mamax;
      end
      mahist(na,:) = mahist(na,:)+histcounts(smpasp(na,:),mabin);
    end
    %
    smpaspid = cha.maidcompress.smpmaid;
    
    ndatapol   = ndatapol   + nch;
    ndataflt   = ndataflt   + nch;
    ndataine   = ndataine   + nch;
    ndataasp   = ndataasp   + nch;
    ndataaspid = ndataaspid + nch;
    sumpol   = sumpol   + sum(smppol,2);
    sumflt   = sumflt   + sum(smpflt,2);
    sumine   = sumine   + sum(smpine,2);
    sumasp   = sumasp   + sum(smpasp,2);
    sumaspid = sumaspid + sum(smpaspid,2);
    sumpolpair   = sumpolpair   + smppol*smppol';
    sumfltpair   = sumfltpair   + smpflt*smpflt';
    suminepair   = suminepair   + smpine*smpine';
    sumasppair   = sumasppair   + smpasp*smpasp';
    sumaspidpair = sumaspidpair + double(smpaspid)*double(smpaspid');
  else
    %
    for np = 1:npol
      infid = cha.mpcompress.npol(np).mpscale==inf;
      if ~infid
        smppol(np,:) = (double(cha.mpcompress.smpmp(np,:))+(sfactor/2))./((sfactor-1).*cha.mpcompress.npol(np).mpscale)+cha.mpcompress.npol(np).mpmin;
      else
        smppol(np,:) = ones(1,nch).*cha.mpcompress.npol(np).mpmax;
      end
    end
    %
    for nf = 1:nflt
      infid = cha.mccompress.nflt(nf).mcscale==inf;
      if ~infid
        smpflt(nf,:) = (double(cha.mccompress.smpmc(nf,:))+(sfactor/2))./((sfactor-1).*cha.mccompress.nflt(nf).mcscale)+cha.mccompress.nflt(nf).mcmin;
      else
        smpflt(nf,:) = ones(1,nch).*cha.mccompress.nflt(nf).mcmax;
      end
    end
    %
    for ni = 1:nine
      infid = cha.micompress.nine(ni).miscale==inf;
      if ~infid
        smpine(ni,:) = (double(cha.micompress.smpmi(ni,:))+(sfactor/2))./((sfactor-1).*cha.micompress.nine(ni).miscale)+cha.micompress.nine(ni).mimin;
      else
        smpine(ni,:) = ones(1,nch).*cha.micompress.nine(ni).mimax;
      end
    end
    %
    for na = 1:nasp
      infid = cha.macompress.nasp(na).mascale==inf;
      if ~infid
        smpasp(na,:) = (double(cha.macompress.smpma(na,:))+(sfactor/2))./((sfactor-1).*cha.macompress.nasp(na).mascale)+cha.macompress.nasp(na).mamin;
      else
        smpasp(na,:) = ones(1,nch).*cha.macompress.nasp(na).mamax;
      end
    end
    %
    smpaspid = cha.maidcompress.smpmaid;
  end
  if accflag
    acctotal = acctotal + cha.ajr;
  end
  asmppol   = [asmppol smppol(:,smpid)];
  asmpflt   = [asmpflt smpflt(:,smpid)];
  asmpine   = [asmpine smpine(:,smpid)];
  asmpasp   = [asmpasp smpasp(:,smpid)];
  asmpaspid = [asmpaspid smpaspid(:,smpid)];
  clear cha
  fprintf('now finised at %i/%i\n',ii,nit)
end
avepol   = sumpol./ndatapol;
aveflt   = sumflt./ndataflt;
aveine   = sumine./ndataine;
aveasp   = sumasp./ndataine;
aveaspid = sumaspid./ndataaspid;
[~,medpolid] = max(mphist,[],2);
[~,medfltid] = max(mchist,[],2);
[~,medineid] = max(mihist,[],2);
[~,medaspid] = max(mahist,[],2);
for np = 1:npol
  medpol(np) = 0.5*(mpbin(medpolid(np))+mpbin(medpolid(np)+1));
end
for nf = 1:nflt
  medflt(nf) = 0.5*(mcbin(medfltid(nf))+mcbin(medfltid(nf)+1));
end
for ni = 1:nine
  medine(ni) = 0.5*(mibin(medineid(ni))+mibin(medineid(ni)+1));
end
for na = 1:nasp
  medasp(na) = 0.5*(mabin(medaspid(na))+mabin(medaspid(na)+1));
end
covpol   = sumpolpair./ndatapol - avepol*avepol';
covflt   = sumfltpair./ndataflt - aveflt*aveflt';
covine   = suminepair./ndataine - aveine*aveine';
covasp   = sumasppair./ndataasp - aveasp*aveasp';
covaspid = sumaspidpair./ndataaspid - aveaspid*aveaspid';
stdpol   = sqrt(diag(covpol));
stdflt   = sqrt(diag(covflt));
stdine   = sqrt(diag(covine));
stdasp   = sqrt(diag(covasp));
stdaspid = sqrt(diag(covaspid));
corpol   = covpol./(stdpol*stdpol');
corflt   = covflt./(stdflt*stdflt');
corine   = covine./(stdine*stdine');
corasp   = covasp./(stdasp*stdasp');
coraspid = covaspid./(stdaspid*stdaspid');

% Save
if accflag
  tcha.acctotal = acctotal;
end
tcha.burnin   = burnin;
tcha.smpint   = smpint;
tcha.mpbin    = mpbin;
tcha.mcbin    = mcbin;
tcha.mibin    = mibin;
tcha.mabin    = mabin;
tcha.avepol   = single(avepol);
tcha.aveflt   = single(aveflt);
tcha.aveine   = single(aveine);
tcha.aveasp   = single(aveasp);
tcha.aveaspid = single(aveaspid);
tcha.medpol   = single(medpol);
tcha.medflt   = single(medflt);
tcha.medine   = single(medine);
tcha.medasp   = single(medasp);
tcha.stdpol   = single(stdpol);
tcha.stdflt   = single(stdflt);
tcha.stdine   = single(stdine);
tcha.stdasp   = single(stdasp);
tcha.stdaspid = single(stdaspid);
tcha.covpol   = single(covpol);
tcha.covflt   = single(covflt);
tcha.covine   = single(covine);
tcha.covasp   = single(covasp);
tcha.covaspid = single(covaspid);
tcha.corpol   = single(corpol);
tcha.corflt   = single(corflt);
tcha.corine   = single(corine);
tcha.corasp   = single(corasp);
tcha.coraspid = single(coraspid);
tcha.histpol  = single(mphist);
tcha.histflt  = single(mchist);
tcha.histine  = single(mihist);
tcha.histasp  = single(mahist);
tcha.smppol   = single(asmppol);
tcha.smpflt   = single(asmpflt);
tcha.smpine   = single(asmpine);
tcha.smpasp   = single(asmpasp);
tcha.smpaspid = single(asmpaspid);
tcha.ndatpol   = single(ndatapol);
tcha.ndatflt   = single(ndataflt);
tcha.ndatine   = single(ndataine);
tcha.ndatasp   = single(ndataasp);
tcha.ndataspid = single(ndataaspid);
end
