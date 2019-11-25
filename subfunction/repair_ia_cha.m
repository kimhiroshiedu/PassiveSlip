function repair_ia_cha(folder)
savedir = fullfile(pwd,folder);
fprintf('Now loading %s ...',fullfile(savedir,'/prm.mat'))
load(fullfile(savedir,'/prm.mat')); fprintf('load\n')
fprintf('Now loading %s ...',fullfile(savedir,'/grn.mat'))
load(fullfile(savedir,'/grn.mat')); fprintf('load\n')
G(1).tb_kin = full(G(1).tb_kin);
G(1).tb_mec = full(G(1).tb_mec);

[input] = read_chafile(savedir);
repair_ia(input,prm,G);

end

function y = Heaviside(x)
% Calculate Heveaside step function
%   y = Heaviside(x)
%   y and x are possible to be either scaler or vector.
y = 1 .* (x >= 0);
end

function repair_ia(input,prm,G)
nit = size(input,2);
% sfactor = 2^8;
sfactor = 2^16;

Hu = repmat(Heaviside(G(1).zulim - G(1).zc),1,prm.nrep);
Hd = repmat(Heaviside(G(1).zdlim - G(1).zc),1,prm.nrep);
Hlim = Hd - Hu;

for nfile = 1:size(input,2)
  load(input(nfile).fname);
  nch = prm.cha;
  nasp = size(cha.macompress.nasp,2);
  smpasp = zeros(nasp,nch);
  for na = 1:nasp
      infid = cha.macompress.nasp(na).mascale==inf;
      if ~infid
          smpasp(na,:) = (double(cha.macompress.smpma(na,:))+(sfactor/2))./((sfactor-1).*cha.macompress.nasp(na).mascale)+cha.macompress.nasp(na).mamin;
      else
          smpasp(na,:) = ones(1,nch).*cha.macompress.nasp(na).mamax;
      end
  end
  ma.smp = smpasp;
  ia.smp = (Heaviside(G(1).zd*ma.smp-G(1).zc) - Heaviside(G(1).zu*ma.smp-G(1).zc)) .* Hlim;
  cha.iacompress.smpia = ia.smp;
  cha.iacompress.meania = mean(int8(ia.smp),2);
  cha.iacompress.covia  = cov(ia.smp');
  fprintf('now finised at %i/%i\n',nfile,nit)
  save(input(nfile).fname,'cha','-v7.3');
end

end

function [input] = read_chafile(folder)
ext  = 'cha_test*.mat';
file = dir([folder,'/',ext]);
[nit,~] = size(file);
for ii = 1:nit
  input(ii).fname = fullfile(file(ii).folder,file(ii).name);
end
end
