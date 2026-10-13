function createOutputDirectory(path)
if isfolder(path) || isfile(path), error('SQ1:Exists','Refusing to overwrite %s',path); end
mkdir(path); mkdir(fullfile(path,'logs')); mkdir(fullfile(path,'failure_logs'));
end
