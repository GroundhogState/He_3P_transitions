data_dir = 'E:\Data\Spectroscopy\spectroscopy_data';
fname = '18hr_CS_lock.txt';
fprintf('Loading...\n')
fid= fopen(fullfile(data_dir,fname),'r');
rawline=fread(fid,Inf,'*char')';
fclose(fid);
ai_dat=jsondecode(rawline);
fprintf('Done.\n')