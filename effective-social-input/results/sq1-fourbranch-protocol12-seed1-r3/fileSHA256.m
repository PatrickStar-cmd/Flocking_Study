function value=fileSHA256(path)
md=java.security.MessageDigest.getInstance('SHA-256');
f=fopen(path,'rb'); assert(f~=-1); cleanup=onCleanup(@()fclose(f));
while ~feof(f), bytes=fread(f,1024*1024,'*uint8'); md.update(bytes); end
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
