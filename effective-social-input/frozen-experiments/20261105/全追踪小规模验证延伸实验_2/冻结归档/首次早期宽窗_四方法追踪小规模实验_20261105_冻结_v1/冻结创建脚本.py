"""Create a new read-only archive without modifying the verified working package."""
from pathlib import Path
import csv
import hashlib
import io
import json
import shutil
import stat
import subprocess
import sys
import zipfile

PARENT=Path(r'C:\Users\WHQ\Studying\大三上学期\无人机群\11.2周\全追踪小规模验证延伸实验_2')
SOURCE=PARENT/'首次早期宽窗_四方法追踪小规模实验_20261105'
ARCHIVES=PARENT/'冻结归档'
VERSION=SOURCE.name+'_冻结_v1'
DESTINATION=ARCHIVES/VERSION
FREEZE_ID='EARLY_SINGLE_WINDOW_TRACKING_20261105_v1'


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def files(root):
    return sorted(p for p in root.rglob('*') if p.is_file() and '__pycache__' not in p.relative_to(root).parts)


def read_csv(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))


VERIFY_SCRIPT = '''"""Read-only frozen-package check; does not run simulations or alter the snapshot."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import stat
import subprocess
import sys
import zipfile

ROOT=Path(__file__).resolve().parent

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--zip',type=Path,help='Optionally verify the adjacent ZIP and its .sha256 file')
    args=parser.parse_args()
    with (ROOT/'FROZEN_MANIFEST.csv').open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
    assert len({r['path'] for r in rows})==len(rows)
    actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file()
            and '__pycache__' not in p.relative_to(ROOT).parts and p!=ROOT/'FROZEN_MANIFEST.csv'}
    assert actual=={r['path'] for r in rows},'Frozen file-set changed'
    for r in rows:
        p=ROOT/r['path']
        assert p.stat().st_size==int(r['bytes']) and sha(p)==r['sha256'],r['path']
    meta=json.loads((ROOT/'冻结索引.json').read_text(encoding='utf-8'))
    snapshot=ROOT/'实验快照'
    original=meta['source_file_sha256']
    snapshot_files={p.relative_to(snapshot).as_posix():sha(p) for p in snapshot.rglob('*')
                    if p.is_file() and '__pycache__' not in p.relative_to(snapshot).parts}
    assert snapshot_files==original,'Snapshot differs from the verified source at freeze time'
    assert sha(snapshot/'experiment_config.json')==meta['config_sha256']
    assert sha(snapshot/'results/analysis/PRIMARY_RESULTS.json')==meta['primary_results_sha256']
    subprocess.run([sys.executable,str(snapshot/'verify_delivery.py')],check=True)
    readonly=sum(not bool(p.stat().st_mode & stat.S_IWRITE) for p in ROOT.rglob('*') if p.is_file())
    zip_checked=False
    if args.zip is not None:
        archive=args.zip.resolve()
        expected=archive.with_name(archive.name+'.sha256').read_text(encoding='utf-8').split()[0]
        assert sha(archive)==expected,'ZIP SHA256 mismatch'
        paths=[p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.relative_to(ROOT).parts]
        with zipfile.ZipFile(archive) as z:
            assert z.testzip() is None,'ZIP CRC check failed'
            expected_names={ROOT.name+'/'+p.relative_to(ROOT).as_posix() for p in paths}
            assert set(z.namelist())==expected_names and len(z.namelist())==len(expected_names)
            for p in paths:
                name=ROOT.name+'/'+p.relative_to(ROOT).as_posix()
                h=hashlib.sha256()
                with z.open(name) as f:
                    for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
                assert h.hexdigest()==sha(p),name
        zip_checked=True
    print('FROZEN PASS:',meta['freeze_id'],';',len(snapshot_files),'source files identical;',len(rows)+1,
          'frozen files;',readonly,'read-only files; ZIP checked:',zip_checked)

if __name__=='__main__':main()
'''


def main():
    assert SOURCE.is_dir()
    archive=ARCHIVES/(VERSION+'.zip')
    assert not DESTINATION.exists() and not archive.exists(),'Existing frozen version will never be overwritten.'
    subprocess.run([sys.executable,str(SOURCE/'verify_delivery.py')],check=True)
    source_hashes={p.relative_to(SOURCE).as_posix():digest(p) for p in files(SOURCE)}
    source_bytes=sum(p.stat().st_size for p in files(SOURCE))
    cfg=json.loads((SOURCE/'experiment_config.json').read_text(encoding='utf-8-sig'))
    primary=json.loads((SOURCE/'results/analysis/PRIMARY_RESULTS.json').read_text(encoding='utf-8'))
    performance=read_csv(SOURCE/'results/analysis/performance_summary.csv')
    paired=read_csv(SOURCE/'results/analysis/paired_summary.csv')
    assert cfg['pilot_seeds']==[971001,971002,971003,971004,971005]
    assert cfg['unique_runs']['total']==65
    assert all(int(r['improved_seeds'])==5 for r in paired if r['scope']=='occlusion_mean'
               and r['contrast'] in ('kf-mean - hold','kf-mean - drop'))
    ARCHIVES.mkdir(exist_ok=True)
    DESTINATION.mkdir()
    snapshot=DESTINATION/'实验快照'
    shutil.copytree(SOURCE,snapshot,ignore=shutil.ignore_patterns('__pycache__'))
    copied={p.relative_to(snapshot).as_posix():digest(p) for p in files(snapshot)}
    assert copied==source_hashes,'Snapshot copy mismatch'
    meta=dict(freeze_id=FREEZE_ID,client_date='2026-11-05',status='FROZEN',version=1,
        source_folder=str(SOURCE),snapshot_folder='实验快照',source_file_count=len(source_hashes),source_bytes=source_bytes,
        source_file_sha256=source_hashes,config_sha256=digest(snapshot/'experiment_config.json'),
        primary_results_sha256=digest(snapshot/'results/analysis/PRIMARY_RESULTS.json'),
        original_delivery_manifest_sha256=digest(snapshot/'DELIVERY_MANIFEST.csv'),
        exact_experiment_config=cfg,primary_results=primary,
        all_fixed_condition_occlusion_performance=[r for r in performance if r['scope']=='occlusion_mean'],
        all_fixed_condition_KF_paired_comparisons=[r for r in paired if r['scope']=='occlusion_mean'
            and r['contrast'] in ('kf-mean - hold','kf-mean - drop')],
        simulation_rerun_during_freeze=False,source_modified_during_freeze=False,
        protection='Independent snapshot plus read-only file attributes and SHA256 manifests; future edits require a new working copy/version.',
        dependencies='Requirements and observed environment are included; external Python runtime and installed dependency directory are not bundled.')
    (DESTINATION/'冻结索引.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    lines=['# 首次早期宽窗四方法追踪实验：冻结 v1','',
        f'冻结编号：`{FREEZE_ID}`；冻结日期：2026年11月5日。','',
        '该目录是本次已完成实验的独立基准存档。冻结时没有重新仿真，没有修改实验条件或结果；实验快照与原工作包逐文件SHA256一致。','',
        '## 冻结的具体实验条件','',
        '- 12架飞机：1领航者、11跟随者。领航者独立追当前目标，跟随者初始航向随机，且不直接接收目标输入。',
        '- 5个种子：971001、971002、971003、971004、971005。总长120步，dt=0.3，为未标定的模型时间单位。',
        '- 仅一次早期遮挡：0基起点20或40，窗宽30或40步，全部4组组合保留；主条件固定为S40/D40。按第1次更新计起点为第21或41次更新。',
        '- 目标轨迹x(t)=235+0.6t，y(t)=200+4sin(2πt/24)；周期性保持原样，没有单独研究或排除其影响。',
        '- 每个跟随者接收1条领航边、6条普通边，幅值分别0.24、0.012；遮挡领航边和3条普通边，共44/77条有向边缺失。',
        '- Fresh、Drop、Hold、KF-mean四方法；KF状态为相对方位和方位变化率，q=R=1e-4，观测噪声标准差0.01rad。',
        '- 遮挡前状态、噪声及目标按种子匹配；未遮挡时直接使用新鲜观测。共65条独立闭环轨迹：5条Fresh、60条处理轨迹。',
        '- 主指标是遮挡期间的跟随者逐机绝对角度误差和到当前目标的欧氏距离；主评分按11跟随者、时间段、种子依次汇总。',
        '- 完整参数以`实验快照/experiment_config.json`为准，运行前方案在`实验快照/EXPERIMENT_PROTOCOL.md`。','',
        '## 冻结的实验结果','',
        '4个预先固定条件中，KF的遮挡期平均角度和距离误差相对Hold、Drop全部为5/5种子改善。角度收益较明确，距离收益较小。','',
        '主条件S40/D40，以下为5种子均值，越小越好：','',
        '| 方法 | 平均角度误差/° | 平均距离误差/模型单位 |','|---|---:|---:|']
    for method in cfg['methods']:
        vals=[next(float(r['mean']) for r in performance if r['start']=='40' and r['duration']=='40'
                   and r['method']==method and r['scope']=='occlusion_mean' and r['metric']==metric)
              for metric in ('angle_deg','distance')]
        lines.append(f'| {method} | {vals[0]:.6f} | {vals[1]:.6f} |')
    lines += ['', '全部条件中KF相对参照降低的遮挡期平均误差（正数为改善）：','',
        '| 条件 | 角度较Hold降低/° | 角度较Drop降低/° | 距离较Hold降低 | 距离较Drop降低 |',
        '|---|---:|---:|---:|---:|']
    for a in (20,40):
        for w in (30,40):
            vals=[-float(next(r['mean'] for r in paired if r['start']==str(a) and r['duration']==str(w)
                and r['scope']=='occlusion_mean' and r['contrast']==contrast and r['metric']==metric))
                for metric,contrast in [('angle_deg','kf-mean - hold'),('angle_deg','kf-mean - drop'),
                                       ('distance','kf-mean - hold'),('distance','kf-mean - drop')]]
            lines.append(f'| S{a}/D{w} | '+' | '.join(f'{v:.6f}' for v in vals)+' |')
    lines += ['', '完整逐种子结果、所有评价时段、样本标准差和改善种子数在实验快照的`results/analysis/`中。',
        '这仍是n=5探索性结果；没有证明长期稳定追踪、任意场景有效或初始随机方向的因果作用，也不能把追踪分数较低解释成方位重构必然更真实。','',
        '## 存档内容与使用方式','',
        '- `实验快照/`：原工作包的完整副本，包括所有代码、配置、协议、65个原始NPZ、65个逐步CSV、汇总结果、4张图、日志和验收记录。',
        '- `冻结索引.json`：完整具体条件、主结果、全部4组结果及原包逐文件指纹。',
        '- `FROZEN_MANIFEST.csv`：整个冻结目录的SHA256清单，清单自身除外。',
        '- `verify_frozen.py`：只读核验冻结目录；加`--zip 压缩包路径`可一并核验ZIP、外部SHA256和压缩包每个条目。',
        '- 同级同名`.zip`和`.zip.sha256`：可整体备份、发送或迁移的冻结包和校验码。',
        '- 代码依赖版本与原运行环境已记录；不包含外部Python安装及依赖安装目录。','',
        '冻结目录及压缩包文件设置为只读，避免误覆盖。后续实验请另建工作副本和版本；复现时可把“实验快照”复制到原实验同级的新目录，以沿用依赖搜索路径。完整归档内的数据不用于继续写入。','',
        '只读核验示例（在冻结目录执行）：','',
        '```powershell','python verify_frozen.py','```','',
        '原工作目录保留原样，本冻结副本作为以后对比、引用与复核的固定基准。','']
    (DESTINATION/'冻结说明_条件与结果.md').write_text('\n'.join(lines),encoding='utf-8')
    (DESTINATION/'verify_frozen.py').write_text(VERIFY_SCRIPT,encoding='utf-8')
    shutil.copyfile(Path(__file__),DESTINATION/'冻结创建脚本.py')
    rows=[dict(path=p.relative_to(DESTINATION).as_posix(),bytes=p.stat().st_size,sha256=digest(p)) for p in files(DESTINATION)]
    with (DESTINATION/'FROZEN_MANIFEST.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['path','bytes','sha256']);writer.writeheader();writer.writerows(rows)
    assert {p.relative_to(SOURCE).as_posix():digest(p) for p in files(SOURCE)}==source_hashes,'Original changed during freeze'
    for p in files(DESTINATION):p.chmod(stat.S_IREAD)
    subprocess.run([sys.executable,str(DESTINATION/'verify_frozen.py')],check=True)
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
        for p in files(DESTINATION):z.write(p,arcname=VERSION+'/'+p.relative_to(DESTINATION).as_posix())
    checksum=archive.with_name(archive.name+'.sha256')
    checksum.write_text(digest(archive)+'  '+archive.name+'\n',encoding='utf-8')
    archive.chmod(stat.S_IREAD);checksum.chmod(stat.S_IREAD)
    subprocess.run([sys.executable,str(DESTINATION/'verify_frozen.py'),'--zip',str(archive)],check=True)
    subprocess.run([sys.executable,str(SOURCE/'verify_delivery.py')],check=True)
    print('FREEZE COMPLETE',json.dumps(dict(freeze_id=FREEZE_ID,folder=str(DESTINATION),zip=str(archive),
        source_files=len(source_hashes),frozen_files=len(rows)+1,zip_bytes=archive.stat().st_size,
        zip_sha256=digest(archive),source_unchanged=True),ensure_ascii=True),flush=True)


if __name__=='__main__':main()
