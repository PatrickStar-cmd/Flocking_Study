function manifest=buildSourceManifest(root,outDir)
files=[dir(fullfile(root,'*.m'));dir(fullfile(root,'tests','*.m'));dir(fullfile(root,'sources','**','*'))];
files=files(~[files.isdir]); rows=struct([]);
for q=1:numel(files)
    path=fullfile(files(q).folder,files(q).name);
    rows(q).relativePath=strrep(path(numel(root)+2:end),'\','/');
    rows(q).actualPath=path; rows(q).SHA256=fileSHA256(path); rows(q).bytes=files(q).bytes;
end
manifest.sourceFiles=rows; manifest.protocolCommit='1e11d64 (protocol freeze); snapshot YCwork@b869447b29a99fbcbafb610ad69b36a8122739f5';
manifest.ACommit='992e71aba3b26896335b52a3182cf01cfd542a9d';
manifest.protocolSHA=fileSHA256(fullfile(root,'sources','protocol','PROTOCOL.md'));
manifest.configSHA=fileSHA256(fullfile(root,'defaultFourBranchConfig.m'));
manifest.implementationVersion='r3 source SHA256 manifest; local files, no Git commit claimed';
for name={'simulateProtocolSegment','controllerInputFromObservation','buildRingConnectivity','localSocialInputMetrics'}
    k=name{1}; manifest.loadedDependencies.(k).path=which(k);
    manifest.loadedDependencies.(k).SHA256=fileSHA256(which(k));
end
writeJson(fullfile(outDir,'source_manifest.json'),manifest);
end
