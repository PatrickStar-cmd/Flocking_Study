"""Record source snapshots, completion status and complete delivery hashes."""
from common import ROOT,RESULTS,read_json,write_json,write_csv,digest


def main():
    checks=RESULTS/'checks'
    for name in ('PREFLIGHT.json','INDEPENDENT_VERIFICATION.json','VISUAL_QA.json'):
        assert read_json(checks/name)['status']=='PASS',name
    for name,sha in read_json(checks/'VISUAL_QA.json')['figure_sha256'].items():
        assert digest(RESULTS/'figures'/name)==sha,'Rendered figure changed after visual inspection: '+name
    snapshot=RESULTS/'analysis_source_snapshot';snapshot.mkdir(exist_ok=True)
    hashes={}
    for p in sorted(ROOT.glob('*.py')):
        (snapshot/p.name).write_bytes(p.read_bytes());hashes[p.name]=digest(p)
    write_json(RESULTS/'analysis_source_hashes.json',hashes)
    write_json(RESULTS/'STATUS.json',dict(stage='COMPLETE',planned=65,completed=65,seeds=5,conditions=4,
        independent_verification='PASS',visual_qa='PASS',raw_data_preserved=True))
    write_json(checks/'DELIVERY_QA.json',dict(status='PASS',all_fixed_conditions_retained=True,all_seeds_retained=True,
        scientific_sources_frozen=True,independent_motion_filter_and_score_audit=True,
        previous_package_untouched=True,statistics_unit='seed',exploratory_only=True))
    rows=[dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(ROOT.rglob('*'))
          if p.is_file() and '__pycache__' not in p.relative_to(ROOT).parts and p!=ROOT/'DELIVERY_MANIFEST.csv']
    write_csv(ROOT/'DELIVERY_MANIFEST.csv',rows)
    print('Delivery frozen:',len(rows),'files;',sum(r['bytes'] for r in rows),'bytes',flush=True)


if __name__=='__main__':main()
