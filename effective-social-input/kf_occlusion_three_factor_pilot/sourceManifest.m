function source=sourceManifest(root)
% 实际实现源码哈希，冻结运行前版本；不假报Git提交。
files=dir(fullfile(root,'*.m')); source=struct([]);
for i=1:numel(files)
    filename=fullfile(root,files(i).name); fid=fopen(filename,'r'); bytes=fread(fid,inf,'*uint8'); fclose(fid);
    md=java.security.MessageDigest.getInstance('SHA-256'); md.update(bytes);
    hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
    row.file=files(i).name; row.sha256=hash; row.bytes=files(i).bytes; row.loadedPath=which(files(i).name(1:end-2));
    if i==1, source=row; else, source(i)=row; end
end
end
