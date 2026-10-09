# R01–R05 与导出修复说明

| 项目 | 实现与验证入口 |
|---|---|
| R01 checkpoint | simulateProtocolSegment及checkpointToStart逐字段赋值，maskDiagnostics作为完整cell字段，不再由struct构造器展开；G4断言两个checkpoint为标量且20个诊断完整，对四模式连续6600步及两级续跑逐步比较 |
| R02 phase | 从occlusionStartStep及occlusionDurationSteps计算；G7覆盖d0/1/10/34/100及d10的6211恢复；d0不生成mask |
| R03 tie | makeProtocolMasks先排序实际ID，再枚举ID元组；1e-12内取字典序第一；G9覆盖观察者17同向tie→1,2,3,4、极值、±pi、刚体旋转、独立随机流及均匀候选抽样 |
| R04 空输入 | 只激活有限方位及正权重；空场/矩数值0并配no_input，方向/熵/gap为NaN；G6/G11覆盖k0/k1/h0/hk/中间/零权重/无缓存hold |
| R05 熵 | weightEntropy恢复为-sum(p log p)，旧值另存weightEntropyNormalized；G11断言ln4/ln8及同向/均匀角等熵 |

补齐逐边lastBearing/lastWeight/observedWeight、神经活动总量/集中度/原始及分析角、隐藏/剩余/全体结构、逐步u、两类oracle、四窗口及配对摘要、真实manifest状态、日志、来源和哈希。每个模式MAT仅保存自身run。

另恢复原神经更新中的globalInhibition和原航向保留阈值1e-9。主配置globalInhibition=0；航向阈值的影响通过G1原核心回归及r2/r3实测对照记录，不以分析方向标记改写运动。

重要版本差异：r2随机mask候选按循环目标索引排序，r3按协议的实际ID字典序排序。两者相同mask随机seed下可能抽到不同ID集合，因此非fresh主轨迹可能改变。本包重新从初始化运行，不复用r2 checkpoint、不拼接或编造缺字段；保留r2/r3对照计数、哈希和差异。此变化不表示控制规则总体收益。
