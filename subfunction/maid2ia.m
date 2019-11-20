function maid2ia(varargin)
if nargin == 1
  mode = 'both';
  folder  = fullfile(pwd,char(varargin{1}));
elseif nargin == 2
  mode = char(varargin{2});
  folder  = fullfile(pwd,char(varargin{1}));
elseif nargin == 0
  folder = input('Enter the result file:  ','s');
else
  error('Too much variables')
end
fin = 0;
while fin ~= 1
  if strcmpi(mode,'both')
    [datin] = read_chafile(folder);
    rename_maid2ia_cha(datin);
    rename_maid2ia_tcha(folder)
    fin = 1;
  elseif strcmpi(mode,'cha')
    [datin] = read_chafile(folder);
    rename_maid2ia_cha(datin);
    fin = 1;
  elseif strcmpi(mode,'tcha')
    rename_maid2ia_tcha(folder);
    fin = 1;
  else
    mode = input('Error. Choose correct option, both/cha/tcha: ','s');
  end
end

end

function rename_maid2ia_tcha(folder)
file = fullfile(folder,'tcha.mat');
load(file);
tcha.aveaid  = tcha.aveaspid ;
tcha.stdaid  = tcha.stdaspid ;
tcha.covaid  = tcha.covaspid ;
tcha.coraid  = tcha.coraspid ;
tcha.smpaid  = tcha.smpaspid ;
tcha.ndataid = tcha.ndataspid;
tcha = rmfield(tcha,{'aveaspid','stdaspid','covaspid','coraspid','smpaspid','ndataspid'});
save(file,'tcha','-v7.3');
end

function rename_maid2ia_cha(datin)
nit = size(datin,2);
for ii = 1:nit
  load(datin(ii).fname); fprintf('load %s ...',datin(ii).fname);
  cha.iacompress.covia = cha.maidcompress.covmaid;
  cha.iacompress.meania = cha.maidcompress.meanmaid;
  cha.iacompress.smpia = cha.maidcompress.smpmaid;
  cha = rmfield(cha,'maidcompress');
  save(datin(ii).fname); fprintf('save %s\n',datin(ii).fname);
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
