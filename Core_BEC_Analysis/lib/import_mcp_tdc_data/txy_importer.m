function  data=txy_importer(filepath,filenum)
%this function imports txy files faster by using strongly a matched import
%double precision
%simply pass the filepath and file number
% 2019-10-31 DKS - filenum will accept integer and int string - error raised
%   for  non-int/str filenum
%% ERROR CHECK
% Check filenum variable type
if isa(filenum,'char')
    if isempty(str2double(filenum))
        error('filenum: %s is not an integer string',filenum);
    else
        filenum=str2double(filenum);
    end
elseif ~isa(filenum,'double')
    % error on data type
    error('filenum is not an int or int str');
end

% Check filenum is an integer
if ~(filenum==floor(filenum))
    error('filenum: %u is not an integer - returning null data',filenum);
end

fileID = fopen([filepath,'_txy_forc',num2str(filenum),'.txt'],'r');
data = textscan(fileID,'%f64%f64%f64',...%f32
    'Delimiter', ',', ...
    'MultipleDelimsAsOne',0,...
    'ReturnOnError', true,...
    'NumCharactersToSkip',0,... 
    'CollectOutput', true,...
    'EndOfLine','\n');
fclose(fileID);
data=data{1};
end
