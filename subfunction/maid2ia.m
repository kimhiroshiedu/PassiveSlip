function maid2ia(folder)
folder  = fullfile(pwd,folder);
[input] = read_chafile(folder);
rename_maid2ia(input);
end

function [input] = read_chafile(folder)
ext  = 'cha_test*.mat';
file = dir([folder,'/',ext]);
[nit,~] = size(file);
for ii = 1:nit
  input(ii).fname = fullfile(file(ii).folder,file(ii).name);
end
end

function rename_maid2ia(input)
nit = size(input,2);
for ii = 1:nit
  load(input(ii).fname); fprintf('load %s ...',input(ii).fname);
  cha.iacompress.covia = cha.maidcompress.covmaid;
  cha.iacompress.meania = cha.maidcompress.meanmaid;
  cha.iacompress.smpia = cha.maidcompress.smpmaid;
  cha = rmfield(cha,'maidcompress');
  save(input(ii).fname); fprintf('save %s\n',input(ii).fname);
end
end