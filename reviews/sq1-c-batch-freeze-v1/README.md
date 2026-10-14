# 四分支开发集执行配置冻结包

日期：2026-10-09；负责人：杨骋。实现负责人：王浩权。

**杨骋当前负责的配置与清单已在此冻结；王浩权可据此开发和检查批处理接口。** r3基础软件已通过独立验收，批处理适配尚未交付，因此配置中允许适配检查，完整560条批次仍待适配验收放行。本次没有启动MATLAB或生产正式实验数据。

## 使用顺序

1. 王浩权先读 [B_HANDOFF.md](B_HANDOFF.md)，实现批处理接口。
2. 运行前执行 `python reviews/sq1-c-batch-freeze-v1/validate_freeze.py` 校验本包；脚本需Python 3标准库及本仓库中对应的Git对象。
3. 读取 [execution_config.json](execution_config.json)、[SOURCE_LOCK.json](SOURCE_LOCK.json) 和两份bindings；按 [ADAPTER_CHECK.json](ADAPTER_CHECK.json) 运行seed 1全部条件。
4. 杨骋按 [BATCH_ACCEPTANCE.md](BATCH_ACCEPTANCE.md) 核查结果，生成新放行记录，固定实际批处理提交SHA后再开始完整批次。

## 文件与含义

| 文件 | 用途 |
|---|---|
| `execution_config.json` | 科学参数、统计方法、调度策略、当前放行状态 |
| `SOURCE_LOCK.json` | 已验收r3提交、43个源文件哈希、15个冻结核心文件与输入来源 |
| `RUN_BINDINGS.json` | 560个计划运行键逐项绑定MATLAB配置、窗口和输出相对路径 |
| `T0_BINDINGS.json` | 80个软件测试计划键，全部仍标PLANNED |
| `ADAPTER_CHECK.json` | 固定seed 1的28个主键和4个T=0键 |
| `references/protocol/` | `YCwork@1e11d64`协议、CSV清单与接口契约的原字节快照 |
| `references/c_acceptance/` | r3独立验收报告、15项测试与原MAT复算证据快照 |
| `SHA256SUMS.json` | 本包文件哈希，排除清单自身与验证记录，避免循环引用 |
| `VALIDATION.json` | 校验实测结果；不代表仿真完成 |

旧协议快照中的`execution.allowed=false`、B历史metadata中的`formalBatchAllowed=false`均保留原样，描述当时状态。当前状态以本包`execution_config.json.authorization`为准。未来完整运行须另存放行记录，不回写历史快照。

## 已冻结的设计

协议科学版本仍为1.2，未因软件验收改变实验设计：N=20、M=100、dt=0.3、周期边界L=1000、k=8、h=4、总权重0.228、单边权重0.0285、固定循环目标ID、无噪声方位、等权输入、无预测、hold不按年龄丢弃。

seed 1–20均为已使用开发种子。掩码代码值为`random / max_hidden_R1 / min_hidden_R1`，后两者展示名分别为concentrated/dispersed。缺失步数为10/34/100；四模式为fresh/drop/drop-normalized/hold。不改变现有CSV的名称和运行键。

每seed 1条共享fresh + 27条非fresh，共560条；另80条T=0软件轨迹。数据仍为探索性开发集，不采用二元成功阈值，不据中途效果删条件。固定k不能单独识别邻居数的因果效应。

## 配对窗口和计数

以下为绝对步号，端点均包含：

| d | 缺失期 | 主要恢复100步 | 延长恢复300步 | 固定事件后400步 |
|---|---|---|---|---|
| 10 | 6201–6210 | 6211–6310 | 6211–6510 | 6201–6600 |
| 34 | 6201–6234 | 6235–6334 | 6235–6534 | 6201–6600 |
| 100 | 6201–6300 | 6301–6400 | 6301–6600 | 6201–6600 |

fresh运行本身无缺失，清单duration=0、mask=not_applicable。比较时用对应非fresh条件的窗口读取同一条fresh参考，不能按fresh的duration=0错位取窗。参考轨迹的phase不作为每个比较条件的事件边界依据。

必须交付540条非fresh运行×4窗口=2160条配对摘要，其中主要终点540条。可另保留720条fresh条件窗口摘要（合计2880行），但须标记比较上下文与共享reference_id，不增加20个独立seed或560条独立运行。两种模式规则差值从同seed、同mask、同d、同window的配对指标计算。

同状态回放须覆盖20seed×3mask×3d=180个逻辑上下文，每个上下文含四候选；可合并存储或共享相同前缀，但不能只回放random/d34。另有60组冻结掩码、40个逻辑检查点；不要求与物理文件数一一对应。

## 验证与版本边界

`prepare_freeze.py`仅生成本配置包，不运行实验；新文件已存在时拒绝覆盖，不应对本已冻结目录再次生成。`validate_freeze.py`读取现有文件，独立枚举运行键、核查逐项配置/窗口/配对来源、Git快照和哈希，并通过5个故障注入检查验证拒绝逻辑。`--record`只用于首次写入新的验证记录，记录存在时拒绝覆盖。

本次工作分支为YCwork，HEAD=`b869447b29a99fbcbafb610ad69b36a8122739f5`；固定基础实现为`54fb85d91151e0b28e8e570318b4a9937d6a0425`。**基础实现SHA不冒充未来批处理实现SHA**；后者目前为null，必须在新放行记录补齐实际提交、差异说明及测试证据。运行器/汇总/回放/审计中的硬编码均属于B适配范围。

只新增本目录，既有验收目录、协议、历史结果和原始代码未修改；没有自动改main，也未提交或推送。main保持`92b8b78cd0752dc91fc96d40f18e5af0a09c953c`。建议将本包作为YCwork协作材料供B使用，当前不建议合并或启动未验收的批处理实现。若撤回，保留历史证据，仅撤销本次新增材料；公开提交后使用`git revert`。
