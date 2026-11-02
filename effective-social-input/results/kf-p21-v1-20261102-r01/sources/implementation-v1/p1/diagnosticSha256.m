function value=diagnosticSha256(path)
fid=fopen(path,'rb'); assert(fid>=0,'Cannot read %s',path);
cleanup=onCleanup(@()fclose(fid));
digest=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid)
    bytes=fread(fid,1024*1024,'*uint8'); digest.update(bytes);
end
value=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
