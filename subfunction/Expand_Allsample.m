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
smpint = 50; % sampling interval
burnin = floor(burnin*nit/100)+1;
acctotal = 0;
for ii = 1:nit
  load(input(ii).fname);
  accflag = isfield(cha,'ajr');
  if ii==1
    nch = size(cha.mpcompress.smpmp,2);
    npol = length(cha.mpcompress.npol);
    nflt = length(cha.mccompress.nflt);
    nine = length(cha.micompress.nine);
    sumpol = zeros(npol,1);
    sumflt = zeros(nflt,1);    
    sumine = zeros(nine,1);    
    sumpolpair = zeros(npol,npol);
    sumfltpair = zeros(nflt,nflt);
    suminepair = zeros(nine,nine);
    ndatapol = 0;
    ndataflt = 0;
    ndataine = 0;
    mphist = zeros(npol,size(mpbin,2)-1);
    mchist = zeros(nflt,size(mcbin,2)-1);
    mihist = zeros(nine,size(mibin,2)-1);
    smpid = [1:smpint:nch];
    smppol = zeros(npol,nch);
    smpflt = zeros(nflt,nch);
    smpine = zeros(nine,nch);
    medpol = zeros(npol,1);
    medflt = zeros(nflt,1);
    medine = zeros(nine,1);
    smppol = [];
    smpflt = [];
    smpine = [];
  end
  if ii>burnin
    for jj = 1:npol
      infid = cha.mpcompress.npol(jj).mpscale==inf;
      if ~infid
        smppol(jj,:) = (double(cha.mpcompress.smpmp(jj,:))+(sfactor/2))./((sfactor-1).*cha.mpcompress.npol(jj).mpscale)+cha.mpcompress.npol(jj).mpmin;
      else
        smppol(jj,:) = ones(1,nch).*cha.mpcompress.npol(jj).mpmax;
      end
      mphist(jj,:) = mphist(jj,:)+histcounts(smppol(jj,:),mpbin);
    end
    for kk = 1:nflt
      infid = cha.mccompress.nflt(kk).mcscale==inf;
      if ~infid
        smpflt(kk,:) = (double(cha.mccompress.smpmc(kk,:))+(sfactor/2))./((sfactor-1).*cha.mccompress.nflt(kk).mcscale)+cha.mccompress.nflt(kk).mcmin;
      else
        smpflt(kk,:) = ones(1,nch).*cha.mccompress.nflt(kk).mcmax;
      end
      mchist(kk,:) = mchist(kk,:)+histcounts(smpflt(kk,:),mcbin);
    end
    for ll = 1:nine
      infid = cha.micompress.nine(ll).miscale==inf;
      if ~infid
        smpine(ll,:) = (double(cha.micompress.smpmi(ll,:))+(sfactor/2))./((sfactor-1).*cha.micompress.nine(ll).miscale)+cha.micompress.nine(ll).mimin;
      else
        smpine(ll,:) = ones(1,nch).*cha.micompress.nine(ll).mimax;
      end
      mihist(ll,:) = mihist(ll,:)+histcounts(smpine(ll,:),mibin);
    end
    ndatapol = ndatapol+nch;
    ndataflt = ndataflt+nch;
    ndataine = ndataine+nch;
    sumpol = sumpol+sum(smppol,2);
    sumflt = sumflt+sum(smpflt,2);
    sumine = sumine+sum(smpine,2);
    sumpolpair = sumpolpair+smppol*smppol';
    sumfltpair = sumfltpair+smpflt*smpflt';
    suminepair = suminepair+smpine*smpine';
  else
    for jj = 1:npol
      infid = cha.mpcompress.npol(jj).mpscale==inf;
      if ~infid
        smppol(jj,:) = (double(cha.mpcompress.smpmp(jj,:))+(sfactor/2))./((sfactor-1).*cha.mpcompress.npol(jj).mpscale)+cha.mpcompress.npol(jj).mpmin;
      else
        smppol(jj,:) = ones(1,nch).*cha.mpcompress.npol(jj).mpmax;
      end
    end
    for kk = 1:nflt
      infid = cha.mccompress.nflt(kk).mcscale==inf;
      if ~infid
        smpflt(kk,:) = (double(cha.mccompress.smpmc(kk,:))+(sfactor/2))./((sfactor-1).*cha.mccompress.nflt(kk).mcscale)+cha.mccompress.nflt(kk).mcmin;
      else
        smpflt(kk,:) = ones(1,nch).*cha.mccompress.nflt(kk).mcmax;
      end
    end
    for ll = 1:nine
      infid = cha.micompress.nine(ll).miscale==inf;
      if ~infid
        smpine(ll,:) = (double(cha.micompress.smpmi(ll,:))+(sfactor/2))./((sfactor-1).*cha.micompress.nine(ll).miscale)+cha.micompress.nine(ll).mimin;
      else
        smpine(ll,:) = ones(1,nch).*cha.micompress.nine(ll).mimax;
      end
    end
  end
  if accflag
    acctotal = acctotal+cha.ajr;
  end
  smppol = [smppol smppol(:,smpid)];
  smpflt = [smpflt smpflt(:,smpid)];
  smpine = [smpine smpine(:,smpid)];
  clear cha
  fprintf('now finised at %i/%i\n',ii,nit)
end
avepol = sumpol./ndatapol;
aveflt = sumflt./ndataflt;
aveine = sumine./ndataine;
[~,medpolid] = max(mphist,[],2);
[~,medfltid] = max(mchist,[],2);
[~,medineid] = max(mihist,[],2);
for np = 1:npol
  medpol(np) = 0.5*(mpbin(medpolid(np))+mpbin(medpolid(np)+1));
end
for nf = 1:nflt
  medflt(nf) = 0.5*(mcbin(medfltid(nf))+mcbin(medfltid(nf)+1));
end
for ni = 1:nine
  medine(ni) = 0.5*(mibin(medineid(ni))+mibin(medineid(ni)+1));
end
covpol = sumpolpair./ndatapol-avepol*avepol';
covflt = sumfltpair./ndataflt-aveflt*aveflt';
covine = suminepair./ndataine-aveine*aveine';
stdpol = sqrt(diag(covpol));
stdflt = sqrt(diag(covflt));
stdine = sqrt(diag(covine));
corpol = covpol./(stdpol*stdpol');
corflt = covflt./(stdflt*stdflt');
corine = covine./(stdine*stdine');
% output
if accflag
  tcha.acctotal = acctotal;
end
tcha.burnin  = burnin;
tcha.smpint  = smpint;
tcha.mpbin   = mpbin;
tcha.mcbin   = mcbin;
tcha.mibin   = mibin;
tcha.avepol  = single(avepol);
tcha.aveflt  = single(aveflt);
tcha.aveine  = single(aveine);
tcha.medpol  = single(medpol);
tcha.medflt  = single(medflt);
tcha.medine  = single(medine);
tcha.stdpol  = single(stdpol);
tcha.stdflt  = single(stdflt);
tcha.stdine  = single(stdine);
tcha.covpol  = single(covpol);
tcha.covflt  = single(covflt);
tcha.covine  = single(covine);
tcha.corpol  = single(corpol);
tcha.corflt  = single(corflt);
tcha.corine  = single(corine);
tcha.histpol = single(mphist);
tcha.histflt = single(mchist);
tcha.histine = single(mihist);
tcha.ndatpol = single(ndatapol);
tcha.ndatflt = single(ndataflt);
tcha.ndatine = single(ndataine);
tcha.smppol  = single(smppol);
tcha.smpflt  = single(smpflt);
tcha.smpine  = single(smpine);
end
