function writeJson(path,value)
f=fopen(path,'w','n','UTF-8'); assert(f~=-1); cleanup=onCleanup(@()fclose(f));
fprintf(f,'%s\n',jsonencode(value,'PrettyPrint',true));
end
