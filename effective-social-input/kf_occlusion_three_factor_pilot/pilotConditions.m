function conditions=pilotConditions()
% 三条单因素扫描共享中心条件，共7种唯一条件。
h=[2 4 6 4 4 4 4]; A=[.228 .228 .228 .228 .228 .114 .456];
label={'random','random','random','concentrated','dispersed','random','random'};
rule={'random','random','random','max_hidden_R1','min_hidden_R1','random','random'};
for i=1:7
    conditions(i)=struct('id',sprintf('h%d_%s_A%d',h(i),label{i},round(A(i)*1000)),...
        'h',h(i),'A0',A(i),'label',label{i},'maskRule',rule{i}); %#ok<AGROW>
end
end
