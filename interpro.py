import time
import json
import requests
from Bio.PDB import PDBParser, PPBuilder
import re
from typing import List, Dict, Any
from goatools.obo_parser import GODag
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from collections import Counter
from typing import List, Dict, Any

# =====================================================
# 函数1️⃣：从PDB文件中提取蛋白质序列
# =====================================================
def extract_sequence_from_pdb(pdb_file, chain_id=None):
    """
    从 PDB 文件中提取蛋白质氨基酸序列
    参数:
        pdb_file (str): PDB文件路径
        chain_id (str, 可选): 指定链ID（默认提取第一个模型的所有链）
    返回:
        str: 蛋白质序列
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    ppb = PPBuilder()

    sequence = ""
    for model in structure:
        for chain in model:
            if chain_id and chain.id != chain_id:
                continue
            for pp in ppb.build_peptides(chain):
                sequence += str(pp.get_sequence())
        break  # 只取第一个 model
    return sequence

# =====================================================
# 函数2️⃣：基于 InterProScan API 的序列功能预测
# =====================================================
def interpro_predict(sequence, email="your_email@example.com", wait_interval=20):
    """
    使用 EBI InterProScan API 对蛋白质序列进行功能预测

    参数:
        sequence (str): 蛋白质氨基酸序列
        email (str): 你的邮箱（EBI要求用于任务标识）
        wait_interval (int): 查询任务状态的间隔秒数

    返回:
        dict: InterProScan 的预测结果(JSON对象)
    """
    # Step 1: 提交任务
    run_url = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run/"
    payload = {
        "email": email,
        "title": "InterProScan_from_PDB",
        "sequence": sequence
    }

    print("[🚀] 提交序列到 InterProScan ...")
    response = requests.post(run_url, data=payload)
    if response.status_code != 200:
        raise Exception(f"提交失败: {response.text}")
    job_id = response.text.strip()
    print(f"[✅] 任务已提交，Job ID: {job_id}")

    # Step 2: 等待任务完成
    status_url = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status/{job_id}"
    while True:
        status = requests.get(status_url).text.strip()
        print(f"[⌛] 当前状态: {status}")
        if status == "FINISHED":
            print("[✅] 任务完成！")
            break
        elif status in ["ERROR", "FAILURE"]:
            raise Exception("InterProScan 任务失败！")
        time.sleep(wait_interval)

    # Step 3: 获取结果（JSON 格式）
    result_url = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/{job_id}/json"
    result_response = requests.get(result_url)
    if result_response.status_code != 200:
        raise Exception(f"结果获取失败: {result_response.text}")

    result = result_response.json()
    print("[💾] 获取预测结果成功！")
    return result

# =====================================================
# 函数3️⃣：整理数据
# =====================================================
def extract_main_data(result: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    从 interproscan-style 的 result 中提取数据，只包含 locations 中 representative=True 的记录。
    每条记录包含：
      - accession (from signature)
      - name (from signature)
      - signature_description (from signature.description)
      - entry_accession (from signature.entry.accession)
      - entry_description (from signature.entry.description)
      - go_terms: list of dicts, each {'id','name','category'} (可为空列表)
      - start, end, score (from location)
      - model_ac (from match['model-ac'])
      - alignment (from location)
    """
    out = []

    results = result.get('results') or []
    if isinstance(results, dict):
        results = [results]

    for res_item in results:
        matches = res_item.get('matches', []) or []
        for match in matches:
            signature = match.get('signature', {}) or {}
            entry = signature.get('entry', {}) or {}

            # locations 通常是列表
            locations = match.get('locations', []) or []

            # 收集所有可能的 goXRefs 来源（列表）
            candidate_lists = []
            if isinstance(entry.get('goXRefs'), list):
                candidate_lists.append(entry.get('goXRefs'))
            if isinstance(signature.get('goXRefs'), list):
                candidate_lists.append(signature.get('goXRefs'))
            if isinstance(match.get('goXRefs'), list):
                candidate_lists.append(match.get('goXRefs'))

            # 将所有 candidate 合并为统一的 go_terms list（去重）
            go_terms: List[Dict[str, Any]] = []
            seen_go_ids = set()
            for gl in candidate_lists:
                for g in gl:
                    if not isinstance(g, dict):
                        continue
                    # 典型字段
                    gid = g.get('id') or g.get('GO') or g.get('goId')
                    gname = g.get('name') or g.get('term') or g.get('label')
                    gcat = g.get('category') or g.get('namespace') or g.get('type')

                    # 如果 id 非标准，尝试在字典值中正则提取 GO:ddddddd
                    if not (isinstance(gid, str) and gid.startswith('GO:')):
                        for v in g.values():
                            if isinstance(v, str):
                                found = re.findall(r'GO:\d{7}', v)
                                if found:
                                    gid = found[0]
                                    break

                    if isinstance(gid, str) and gid.startswith('GO:'):
                        if gid not in seen_go_ids:
                            seen_go_ids.add(gid)
                            go_terms.append({
                                'id': gid,
                                'name': gname,
                                'category': gcat
                            })

            # 如果没有在 candidate_lists 找到 goXRefs，但 entry 或 signature 里可能包含 GO 字符串（备用查找）
            if not go_terms:
                # 在 entry 的任意字符串字段中找 GO:xxxxxx
                for container in (entry, signature, match):
                    for v in (container or {}).values():
                        if isinstance(v, str) and 'GO:' in v:
                            founds = re.findall(r'GO:\d{7}', v)
                            for f in founds:
                                if f not in seen_go_ids:
                                    seen_go_ids.add(f)
                                    go_terms.append({'id': f, 'name': None, 'category': None})
                # 仍为空则保留空列表

            # 对每个 representative location 生成一条记录
            for loc in locations:
                if not isinstance(loc, dict):
                    continue
                if not loc.get('representative', False):
                    continue

                record = {
                    'accession': signature.get('accession'),
                    'name': signature.get('name'),
                    'entry_accession': entry.get('accession'),
                    'entry_description': entry.get('description'),
                    'go_terms': go_terms.copy(),   # 每条记录独立拷贝
                    'start': loc.get('start'),
                    'end': loc.get('end'),
                    'score': loc.get('score'),
                }
                out.append(record)

    return out

# =====================================================
# 函数4️⃣：富集分析
# =====================================================
def go_analysis(records: List[Dict[str, Any]], obo_file: str = "go-basic.obo"):
    """
    综合 GO 分析：
    - 单蛋白/少量 GO 时：统计并打印 GO Term 分布
    - 多蛋白/足够 GO 时：做富集分析并打印显著 GO Term

    参数:
        records (List[Dict]): extract_main_data 输出的记录列表
        obo_file (str): GO DAG 文件路径
    """
    # 收集 GO 注释
    protein_to_go = {}  # protein_index -> list of GO IDs
    for idx, rec in enumerate(records):
        go_ids = [go['id'] for go in rec.get('go_terms', []) if 'id' in go]
        if go_ids:
            protein_to_go[f"protein_{idx + 1}"] = go_ids

    all_go_ids = [go for gos in protein_to_go.values() for go in gos]
    if not all_go_ids:
        print("[⚠️] 没有找到 GO 注释")
        return

    # 加载 GO DAG
    godag = GODag(obo_file)

    # 判断是否做富集分析
    if len(protein_to_go) < 2 or len(all_go_ids) < 5:
        # 少量蛋白或 GO，直接统计输出
        go_counter = Counter(all_go_ids)
        print(f"[🧬] GO Term 分布统计 (少量蛋白/GO，不做富集分析)：共 {len(go_counter)} 个 GO Term")
        for go_id, count in go_counter.most_common():
            term = godag.get(go_id)
            name = term.name if term else "Unknown"
            namespace = term.namespace if term else "Unknown"
            print(f"{go_id} | {name} | namespace: {namespace} | count: {count}")
        return

    # 多蛋白/GO，做富集分析
    # 构建背景集：所有蛋白的 GO
    population = set(all_go_ids)
    assoc = {p: gos for p, gos in protein_to_go.items()}
    study = list(population)  # 简化处理，用所有 GO 做研究集

    # 初始化富集分析
    goeaobj = GOEnrichmentStudyNS(
        list(population),
        assoc,
        godag,
        propagate_counts=True,
        alpha=0.05,
        methods=["fdr_bh"]
    )

    # 执行分析
    goea_results_all = goeaobj.run_study(study)

    # 显著 GO Term
    sig_results = [r for r in goea_results_all if r.p_fdr_bh < 0.05]

    if not sig_results:
        print("[ℹ️] 没有显著富集 GO Term")
        return

    # 输出结果
    print(f"[🧪] 显著富集 GO Term (FDR < 0.05)：共 {len(sig_results)} 个")
    for r in sig_results:
        print(f"{r.GO} | {r.name} | namespace: {r.NS} | p_fdr_bh: {r.p_fdr_bh:.4g} | study_count: {r.study_count}")


# =====================================================
# 示例用法
# =====================================================
if __name__ == "__main__":
    pdb_path = "3cdd1.pdb"  # 你的PDB文件路径
    # 从PDB中提取序列
    seq = extract_sequence_from_pdb(pdb_path)
    print(f"[🧬] 提取到的序列长度: {len(seq)}")
    # 使用InterPro预测功能
    result = interpro_predict(seq, email="2805238699@qq.com")
    # 打印结果摘要
    print(extract_main_data(result))
    # 进行富集分析
    go_analysis(extract_main_data(result))