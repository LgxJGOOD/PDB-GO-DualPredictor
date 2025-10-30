"""
一体化 PDB 双通道 GO 预测 & 精确+语义相似性对比
主函数：analyze_pdb(pdb_path) -> dict
无需 Wang 模块，仅依赖 goatools 最短路径
"""
from typing import Dict, Set, List, Tuple
import os
import time
import re
import requests
from Bio.PDB import PDBParser, PPBuilder
from goatools.obo_parser import GODag

# ---------------- 常量 ----------------
INTERPRO_WAIT = 20
DEEPFRI_WORKSPACE = "HE93D7"
GO_OBO = "go-basic.obo"     # 需存在于当前目录或指定路径
godag = GODag(GO_OBO)
# -------------------------------------


# ============ 工具函数 ============
def _extract_seq(pdb_path: str) -> str:
    """从 PDB 提取第一条多肽链序列"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("X", pdb_path)
    ppb = PPBuilder()
    seq = ""
    for model in structure:
        for chain in model:
            for pp in ppb.build_peptides(chain):
                seq += str(pp.get_sequence())
        break
    return seq


def _interpro(seq: str, email: str) -> Set[str]:
    """提交 InterProScan 并返回 GO 集合"""
    run_url = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run/"
    resp = requests.post(run_url, data={"email": email, "title": "X", "sequence": seq})
    resp.raise_for_status()
    job_id = resp.text.strip()

    status_url = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status/{job_id}"
    while True:
        status = requests.get(status_url).text.strip()
        if status == "FINISHED":
            break
        if status in {"ERROR", "FAILURE"}:
            raise RuntimeError("InterProScan failed")
        time.sleep(INTERPRO_WAIT)

    result = requests.get(
        f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/{job_id}/json"
    ).json()

    go_set = set()
    results = result.get("results", [])
    if isinstance(results, dict):
        results = [results]
    for res in results:
        for m in res.get("matches", []):
            entry = m.get("signature", {}).get("entry") or {}
            for g in entry.get("goXRefs", []):
                go_id = g.get("id") or ""
                if re.match(r"GO:\d{7}", go_id):
                    go_set.add(go_id)
    return go_set


def _deepfri(pdb_path: str) -> Dict[str, float]:
    """上传 PDB → DeepFRI → 返回 GO:score 字典"""
    session = requests.Session()
    upload_url = f"https://beta.api.deepfri.flatironinstitute.org/workspace/{DEEPFRI_WORKSPACE}/predictions"
    headers = {
        "accept": "application/json",
        "origin": "https://beta.deepfri.flatironinstitute.org",
        "referer": "https://beta.deepfri.flatironinstitute.org/",
        "user-agent": "Mozilla/5.0",
    }
    with open(pdb_path, "rb") as f:
        files = {"file": (os.path.basename(pdb_path), f, "application/octet-stream")}
        data = {"inputType": "structureFile", "tags": ""}
        up_resp = session.post(upload_url, files=files, data=data, headers=headers)
    up_resp.raise_for_status()
    task_name = up_resp.json()["predictions"][0]["name"]

    task_url = f"https://beta.api.deepfri.flatironinstitute.org/workspace/{DEEPFRI_WORKSPACE}/predictions/{task_name}"
    while True:
        info = session.get(task_url, headers=headers).json()
        state = (info.get("prediction") or info.get("predictions", [{}])[0]).get("state", "")
        if state == "finished":
            final_data = (info.get("prediction") or info.get("predictions", [{}])[0]).get("data", {})
            break
        if state == "failed":
            raise RuntimeError("DeepFRI failed")
        time.sleep(2)

    go_score = {}
    for model in ["cnn_cc", "cnn_mf", "cnn_bp"]:
        for p in final_data.get("A", {}).get(model, {}).get("predictions", []):
            go_score[p["go_term"]] = float(p["go_term_score"])
    return go_score


def _go_distance(go1: str, go2: str) -> int:
    """
    计算两个 GO term 在 DAG 中的最短路径距离（边数）。
    若两个 term 相同则距离为 0；
    若无法连通则返回 None。
    """
    if go1 == go2:
        return 0
    # 把字符串转成 GODag 节点
    a = godag.get(go1)
    b = godag.get(go2)
    if a is None or b is None:
        return None

    # 收集 a 的所有祖先（含自己）
    a_anc = {go1} | a.get_all_parents()
    # 收集 b 的所有祖先（含自己）
    b_anc = {go2} | b.get_all_parents()
    # 公共祖先集合
    common = a_anc & b_anc
    if not common:
        return None

    # 对每个公共祖先计算 (a 到祖先步数 + b 到祖先步数)
    def _min_steps(start_go, target_set):
        """BFS 求最短步数"""
        from collections import deque
        queue = deque([(start_go, 0)])
        seen = {start_go}
        while queue:
            cur, step = queue.popleft()
            if cur in target_set:
                return step
            node = godag.get(cur)
            if node is None:
                continue
            for par in node.parents:
                if par.id not in seen:
                    seen.add(par.id)
                    queue.append((par.id, step + 1))
        return None

    best = None
    for ca in common:
        da = _min_steps(go1, {ca})
        db = _min_steps(go2, {ca})
        if da is not None and db is not None:
            best = min(best, da + db) if best is not None else da + db
    return best

def _semantic_pairs(setA: Set[str], setB: Set[str],
                    sim_threshold: float = 0.7) -> Tuple[float, List[Tuple[str, str, float]]]:
    """
    基于自写最短路径的归一化相似度
    返回 (avg_score, [(goA, goB, score), ...] 高于阈值)
    """
    scores = []
    for a in setA:
        best = 0.0
        pair = (a, "", 0.0)
        for b in setB:
            d = _go_distance(a, b)
            sim = 1 / (1 + d) if d is not None else 0.0
            if sim > best:
                best = sim
                pair = (a, b, sim)
        if best > 0:
            scores.append(pair)
    if not scores:
        return 0.0, []
    avg_score = sum(s for _, _, s in scores) / len(scores)
    high_pairs = [(a, b, s) for a, b, s in scores if s >= sim_threshold]
    return avg_score, high_pairs


# ============ 主函数 ============
def analyze_pdb(pdb_path: str, email: str = "2805238699@qq.com",
                sem_sim_threshold: float = 0.7) -> Dict:
    """
    输入 PDB 文件路径，返回含精确+语义一致性的完整对比字典
    """
    seq = _extract_seq(pdb_path)
    ipro_go = _interpro(seq, email)
    df_go = _deepfri(pdb_path)

    # 精确匹配
    common = ipro_go & df_go.keys()
    ipro_only = ipro_go - df_go.keys()
    df_only = df_go.keys() - ipro_go
    jaccard = len(common) / len(ipro_go | df_go.keys()) if (ipro_go | df_go.keys()) else 0.0
    high_conf_common = {go for go in common if df_go[go] >= 0.8}

    # 语义相似性
    avg_sem_sim, sem_pairs = _semantic_pairs(ipro_go, df_go.keys(), sim_threshold=sem_sim_threshold)

    return {
        "interpro_count": len(ipro_go),
        "deepfri_count": len(df_go),
        "common": common,
        "common_count": len(common),
        "high_confidence_common": high_conf_common,
        "high_confidence_count": len(high_conf_common),
        "interpro_only": ipro_only,
        "interpro_only_count": len(ipro_only),
        "deepfri_only": df_only,
        "deepfri_only_count": len(df_only),
        "jaccard": jaccard,
        "deepfri_scores": df_go,
        # 语义指标
        "semantic_avg_sim": avg_sem_sim,
        "semantic_pairs": sem_pairs,
        "semantic_pairs_count": len(sem_pairs),
    }


# ============ 示例调用 ============
if __name__ == "__main__":
    pdb = "1AKI.pdb"  # 换成你的 PDB
    print(analyze_pdb(pdb))