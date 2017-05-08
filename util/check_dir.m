function check_dir(folder)

if ~exist(folder, 'dir')
  mkdir(folder);
end