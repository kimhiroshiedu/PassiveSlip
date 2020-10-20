function Expand_Allsample(folder,burnin)
% burnin: enter for percent scale
folder  = fullfile(pwd,folder);
load(fullfile(folder,'prm.mat'));
[input] = read_chafile(folder);
[tcha]  = cal_avestdbin(input,burnin,prm);

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
function [tcha] = cal_avestdbin(input,burnin,prm)
nit = size(input,2);
% sfactor = 2^8;
sfactor = 2^16;
smpint = 50; % resampling interval
burnin = floor(burnin*nit/100) + 1;
acctotal = 0;
for ii = 1:nit
  load(input(ii).fname);
  accflag = isfield(cha,'ajr');
  if ii==1
    nrep = prm.nrep;
    nch = size(cha.mpcompress.smpmp,2);
    npol = size(cha.mpcompress.npol,2);
    nflt = size(cha.mccompress.nflt,2);
    nine = size(cha.micompress.nine,2);
    nasp = size(cha.macompress.nasp,2);
    naid = size(cha.iacompress.smpia,1);
    sumpol = zeros(npol,1,nrep);
    sumflt = zeros(nflt,1,nrep);
    sumine = zeros(nine,1,nrep);
    sumasp = zeros(nasp,1,nrep);
    sumaid = zeros(naid,1,nrep);
    sumpolpair = zeros(npol,npol,nrep);
    sumfltpair = zeros(nflt,nflt,nrep);
    suminepair = zeros(nine,nine,nrep);
    sumasppair = zeros(nasp,nasp,nrep);
    sumaidpair = zeros(naid,naid,nrep);
    ndatapol = 0;
    ndataflt = 0;
    ndataine = 0;
    ndataasp = 0;
    ndataaid = 0;
    smpid = 1:smpint:nch;
    asmppol = [];
    asmpflt = [];
    asmpine = [];
    asmpasp = [];
    asmpaid = [];
  end
  %
  smppol = zeros(npol,nch,nrep);
  smpflt = zeros(nflt,nch,nrep);
  smpine = zeros(nine,nch,nrep);
  smpasp = zeros(nasp,nch,nrep);
  % Euler vector
  for np = 1:npol
      infid = cha.mpcompress.npol(np).mpscale==inf;
      smppol(np,:,~infid) = (double(cha.mpcompress.smpmp(np,:,~infid))+(sfactor/2))./((sfactor-1).*cha.mpcompress.npol(np).mpscale(:,:,~infid))+cha.mpcompress.npol(np).mpmin(:,:,~infid);
      smppol(np,:, infid) = ones(1,nch,sum(infid)).*cha.mpcompress.npol(np).mpmax(:,:,infid);
  end
  % Coupling ratio
  for nf = 1:nflt
      infid = cha.mccompress.nflt(nf).mcscale==inf;
      smpflt(nf,:,~infid) = (double(cha.mccompress.smpmc(nf,:,~infid))+(sfactor/2))./((sfactor-1).*cha.mccompress.nflt(nf).mcscale(:,:,~infid))+cha.mccompress.nflt(nf).mcmin(:,:,~infid);
      smpflt(nf,:, infid) = ones(1,nch,sum(infid)).*cha.mccompress.nflt(nf).mcmax(:,:,infid);
  end
  % Internal strain rate
  for ni = 1:nine
      infid = cha.micompress.nine(ni).miscale==inf;
      smpine(ni,:,~infid) = (double(cha.micompress.smpmi(ni,:,~infid))+(sfactor/2))./((sfactor-1).*cha.micompress.nine(ni).miscale(:,:,~infid))+cha.micompress.nine(ni).mimin(:,:,~infid);
      smpine(ni,:, infid) = ones(1,nch,sum(infid)).*cha.micompress.nine(ni).mimax(:,:,infid);
  end
  % Depth limit of asperities
  for na = 1:nasp
      infid = cha.macompress.nasp(na).mascale==inf;
      smpasp(na,:,~infid) = (double(cha.macompress.smpma(na,:,~infid))+(sfactor/2))./((sfactor-1).*cha.macompress.nasp(na).mascale(:,:,~infid))+cha.macompress.nasp(na).mamin(:,:,~infid);
      smpasp(na,:, infid) = ones(1,nch,sum(infid)).*cha.macompress.nasp(na).mamax(:,:,infid);
  end
  % Asperity binaries
  smpaid = cha.iacompress.smpia;
  %
  if ii>burnin
    ndatapol = ndatapol + nch;
    ndataflt = ndataflt + nch;
    ndataine = ndataine + nch;
    ndataasp = ndataasp + nch;
    ndataaid = ndataaid + nch;
    sumpol = sumpol + sum(smppol,2);
    sumflt = sumflt + sum(smpflt,2);
    sumine = sumine + sum(smpine,2);
    sumasp = sumasp + sum(smpasp,2);
    sumaid = sumaid + sum(smpaid,2);
    for rep = 1:nrep
      sumpolpair(:,:,rep) = sumpolpair(:,:,rep) + smppol(:,:,rep)*smppol(:,:,rep)';
      sumfltpair(:,:,rep) = sumfltpair(:,:,rep) + smpflt(:,:,rep)*smpflt(:,:,rep)';
      suminepair(:,:,rep) = suminepair(:,:,rep) + smpine(:,:,rep)*smpine(:,:,rep)';
      sumasppair(:,:,rep) = sumasppair(:,:,rep) + smpasp(:,:,rep)*smpasp(:,:,rep)';
      sumaidpair(:,:,rep) = sumaidpair(:,:,rep) + double(smpaid(:,:,rep))*double(smpaid(:,:,rep)');
    end
  end
  %
  if accflag
    acctotal = acctotal + cha.ajr;
  end
  asmppol = [asmppol smppol(:,smpid,:)];
  asmpflt = [asmpflt smpflt(:,smpid,:)];
  asmpine = [asmpine smpine(:,smpid,:)];
  asmpasp = [asmpasp smpasp(:,smpid,:)];
  asmpaid = [asmpaid smpaid(:,smpid,:)];
  clear cha
  fprintf('now finised at %i/%i\n',ii,nit)
end
avepol = sumpol./ndatapol;
aveflt = sumflt./ndataflt;
aveine = sumine./ndataine;
aveasp = sumasp./ndataine;
aveaid = sumaid./ndataaid;
medpol = median(asmppol,2);
medflt = median(asmpflt,2);
medine = median(asmpine,2);
medasp = median(asmpasp,2);
medaid = median(asmpaid,2);
for rep = 1:nrep
  covpol(:,:,rep) = sumpolpair(:,:,rep)./ndatapol - avepol(:,:,rep)*avepol(:,:,rep)';
  covflt(:,:,rep) = sumfltpair(:,:,rep)./ndataflt - aveflt(:,:,rep)*aveflt(:,:,rep)';
  covine(:,:,rep) = suminepair(:,:,rep)./ndataine - aveine(:,:,rep)*aveine(:,:,rep)';
  covasp(:,:,rep) = sumasppair(:,:,rep)./ndataasp - aveasp(:,:,rep)*aveasp(:,:,rep)';
  covaid(:,:,rep) = sumaidpair(:,:,rep)./ndataaid - aveaid(:,:,rep)*aveaid(:,:,rep)';
  stdpol(:,:,rep) = sqrt(diag(covpol(:,:,rep)));
  stdflt(:,:,rep) = sqrt(diag(covflt(:,:,rep)));
  stdine(:,:,rep) = sqrt(diag(covine(:,:,rep)));
  stdasp(:,:,rep) = sqrt(diag(covasp(:,:,rep)));
  stdaid(:,:,rep) = sqrt(diag(covaid(:,:,rep)));
  corpol(:,:,rep) = covpol(:,:,rep)./(stdpol(:,:,rep)*stdpol(:,:,rep)');
  corflt(:,:,rep) = covflt(:,:,rep)./(stdflt(:,:,rep)*stdflt(:,:,rep)');
  corine(:,:,rep) = covine(:,:,rep)./(stdine(:,:,rep)*stdine(:,:,rep)');
  corasp(:,:,rep) = covasp(:,:,rep)./(stdasp(:,:,rep)*stdasp(:,:,rep)');
  coraid(:,:,rep) = covaid(:,:,rep)./(stdaid(:,:,rep)*stdaid(:,:,rep)');
end

% Save
if accflag
  tcha.acctotal = acctotal;
end
tcha.burnin = burnin;
tcha.nit    = nit;
tcha.smpint = smpint;
tcha.avepol = single(avepol);
tcha.aveflt = single(aveflt);
tcha.aveine = single(aveine);
tcha.aveasp = single(aveasp);
tcha.aveaid = single(aveaid);
tcha.medpol = single(medpol);
tcha.medflt = single(medflt);
tcha.medine = single(medine);
tcha.medasp = single(medasp);
tcha.medaid = single(medaid);
tcha.stdpol = single(stdpol);
tcha.stdflt = single(stdflt);
tcha.stdine = single(stdine);
tcha.stdasp = single(stdasp);
tcha.stdaid = single(stdaid);
tcha.covpol = single(covpol);
tcha.covflt = single(covflt);
tcha.covine = single(covine);
tcha.covasp = single(covasp);
tcha.covaid = single(covaid);
tcha.corpol = single(corpol);
tcha.corflt = single(corflt);
tcha.corine = single(corine);
tcha.corasp = single(corasp);
tcha.coraid = single(coraid);
tcha.smppol = single(asmppol);
tcha.smpflt = single(asmpflt);
tcha.smpine = single(asmpine);
tcha.smpasp = single(asmpasp);
tcha.smpaid = single(asmpaid);
tcha.ndatpol = single(ndatapol);
tcha.ndatflt = single(ndataflt);
tcha.ndatine = single(ndataine);
tcha.ndatasp = single(ndataasp);
tcha.ndataid = single(ndataaid);
end
