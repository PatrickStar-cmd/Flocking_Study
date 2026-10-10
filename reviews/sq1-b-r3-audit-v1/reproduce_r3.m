% Independent new software output; never edit B's delivered results.
auditRoot=fileparts(mfilename('fullpath'));src=fullfile(auditRoot,'source');
restoredefaultpath;addpath(src);
out=run_protocol12_seed1_r3('C_recheck_20261009');diary('off');
for name={'tests.json','audit.json','metadata.json','resume_check.csv','T0_maxdiff.csv','weight_contract.csv','summary.csv'}
 copyfile(fullfile(out,name{1}),fullfile(auditRoot,['C_' name{1}]));
end
fprintf('C independent software recheck generated in a new directory: %s\n',out);
