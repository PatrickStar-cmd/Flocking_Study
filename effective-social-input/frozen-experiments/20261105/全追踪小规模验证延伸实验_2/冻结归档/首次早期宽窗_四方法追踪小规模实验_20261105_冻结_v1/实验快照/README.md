# 首次早期宽窗四方法追踪实验

已完成65条闭环轨迹。优先阅读`实验结果说明.md`。
条件：12架、1领航者、5个新种子、起点20/40、单窗宽30/40、120步；原目标轨迹与KF参数不变。
主条件S40/D40，n=5探索；四个条件遮挡均值的角度与距离均为KF优于Hold和Drop，距离收益较小。

## 运行与复核

PowerShell中输入：

```powershell
cd 'C:\Users\WHQ\Studying\大三上学期\无人机群\11.2周\全追踪小规模验证延伸实验_2\首次早期宽窗_四方法追踪小规模实验_20261105'
& 'C:\Users\WHQ\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' run_all.py
```

已有结果会校验哈希并跳过已完成仿真；如果改动科学参数或引擎，必须建立新版本。
脚本使用`requirements.lock.txt`的依赖；bootstrap会寻找本目录.runtime_deps或`11.2周/追踪验证/.runtime_deps`。
仿真代码不依赖旧目录，旧目录存在时会附加核对来源和旧种子回归；旧种子不纳入实验样本。
`run_all.py`重新生成的图形必须与已经目视检查的图形哈希一致，否则最终验收会要求重新查看。

只读完整性核验：

```powershell
& 'C:\Users\WHQ\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' verify_delivery.py
```

保留`results/checks/VISUAL_QA.json`，其中记录本轮已目视检查的图形指纹；若重新作图有变化，应重新检查。
`run.log`为原始仿真日志，`independent_audit.log`为独立控制与指标重放日志。
角度和距离是跟随者逐机指标，模型长度/时间没有换算为米/秒。
