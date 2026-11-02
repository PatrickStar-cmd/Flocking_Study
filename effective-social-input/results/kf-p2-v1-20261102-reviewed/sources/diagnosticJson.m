function diagnosticJson(path,value)
% For metadata owned by the current newly-created run only.
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0,'Cannot write %s',path);
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
