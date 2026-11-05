# 两项冻结实验归档（2026-11-05）

完整保留两个本地“冻结归档”目录，包括实验条件、源码、原始结果、汇总、图表、校验和与ZIP。原文件未修改。

- `全追踪小规模验证延伸实验_1/冻结归档/强领航随机初向_FreshDrop验证_20261104_冻结_v1/`：强领航随机初向Fresh/Drop实验，详细条件和结果在归档ZIP及冻结说明中，ZIP包含592个载荷文件。
- `全追踪小规模验证延伸实验_2/冻结归档/首次早期宽窗_四方法追踪小规模实验_20261105_冻结_v1/`：早期单次宽窗四方法实验。12架、1领航者，种子971001—971005，起点20/40、窗宽30/40、120步。全部65条轨迹和4组结果保留。主条件S40/D40下KF角度误差较Hold低2.53°、较Drop低9.89°；距离收益较小。n=5探索性结果，不声称长期稳定或普遍有效。

`UPLOAD_INVENTORY.csv`记录所有214个本地归档文件的仓库路径、大小和SHA256，可核对下载内容。实验自身的冻结清单和校验脚本原样保留。

## 下载与校验

第一个实验的ZIP超过GitHub普通Git单文件100 MiB限制，仅该文件使用Git LFS存储；其余归档文件使用普通Git。

推荐安装Git LFS后：

```bash
git lfs install
git clone --branch WHQ_work --single-branch https://github.com/PatrickStar-cmd/Flocking_Study.git
cd Flocking_Study
git lfs pull
```

网页单独下载大ZIP时，请使用文件页的原始下载按钮；若得到几行`version https://git-lfs.github.com/spec/v1`文本，说明下载的是指针，应通过Git LFS获取真实ZIP。

第一个实验运行`校验归档.py`；第二个实验运行`verify_frozen.py`，可用`--zip`指定同级ZIP进一步校验。Python运行环境/第三方安装目录未随归档打包，版本和依赖要求在快照内。

归档内保留了历史绝对路径用于来源记录；复现实验应先复制到新的工作目录并按说明配置依赖，不覆盖冻结版本。
