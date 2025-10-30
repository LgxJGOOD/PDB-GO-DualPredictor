#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DeepFri 上传 → JSON → GO 富集分析 一键脚本
"""

import os
import time
import json
import requests
import pandas as pd
from typing import Dict, Any

# GO 分析部分
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy

# ========== 1. 上传 & 轮询 ==========
def upload_and_get_result(file_path: str,
                          workspace_id: str = "HE93D7",
                          max_wait: int = 120,
                          poll_interval: int = 2,
                          save_json: bool = True,
                          output_dir: str = ".") -> Dict[str, Any]:
    session = requests.Session()
    upload_url = f"https://beta.api.deepfri.flatironinstitute.org/workspace/{workspace_id}/predictions"
    headers = {
        "accept": "application/json, text/plain, */*",
        "origin": "https://beta.deepfri.flatironinstitute.org",
        "referer": "https://beta.deepfri.flatironinstitute.org/",
        "user-agent": "Mozilla/5.0"
    }

    with open(file_path, "rb") as f:
        files = {"file": (os.path.basename(file_path), f, "application/octet-stream")}
        data = {"inputType": "structureFile", "tags": ""}
        upload_resp = session.post(upload_url, headers=headers, files=files, data=data)
    upload_resp.raise_for_status()
    predictions = upload_resp.json().get("predictions", [])
    if not predictions:
        raise RuntimeError("上传失败，未返回任务信息")
    task_name = predictions[0]["name"]
    print(f"[上传成功] 任务名: {task_name}")

    task_url = f"https://beta.api.deepfri.flatironinstitute.org/workspace/{workspace_id}/predictions/{task_name}"
    start = time.time()
    while True:
        resp = session.get(task_url, headers=headers)
        resp.raise_for_status()
        info = resp.json().get("prediction") or resp.json().get("predictions", [{}])[0]
        state = info.get("state", "")
        print(f"[轮询] 状态: {state}")
        if state == "finished" and info.get("data"):
            final_data = info["data"]
            break
        if state == "failed":
            raise RuntimeError("任务失败")
        if time.time() - start > max_wait:
            raise TimeoutError("等待超时")
        time.sleep(poll_interval)

    if save_json:
        out_json = os.path.join(output_dir, f"{task_name}_result.json")
        with open(out_json, "w", encoding="utf-8") as f:
            json.dump(final_data, f, ensure_ascii=False, indent=2)
        print(f"[保存] JSON 已写入 {out_json}")
    return final_data


# ========== 2. 提取 GO 列表 ==========
def extract_go_predictions(json_data: Dict[str, Any],
                           score_threshold: float = 0.0) -> pd.DataFrame:
    ns_map = {"cnn_cc": "CC", "cnn_mf": "MF", "cnn_bp": "BP", "cnn_ec": "EC"}
    records = []
    for model, ns in ns_map.items():
        for p in json_data.get("A", {}).get(model, {}).get("predictions", []):
            score = p.get("go_term_score", 0.0)
            if score >= score_threshold:
                records.append({
                    "go_term": p["go_term"],
                    "name": p["go_term_name"],
                    "score": score,
                    "namespace": ns
                })
    df = pd.DataFrame(records).drop_duplicates(subset=["go_term"])
    return df


# ========== 3. GO 富集分析 ==========
def go_enrichment(df: pd.DataFrame,
                  obo_path: str = "go-basic.obo",
                  alpha: float = 0.05) -> pd.DataFrame:
    """
    对 DeepFri 预测的 GO 做富集分析（单蛋白演示）
    返回显著条目 DataFrame
    """
    # 1. 构建 GO  DAG
    godag = GODag(obo_path)

    # 2. 基因列表 & GO 映射（这里仅演示单蛋白）
    study_genes = {"PROTEIN_A"}
    gene2go = {"PROTEIN_A": set(df["go_term"])}
    population = {"PROTEIN_A"}  # 背景同样本，仅演示

    goea = GOEnrichmentStudy(
        population, gene2go, godag,
        propagate_counts=False, alpha=alpha, methods=["fdr_bh"]
    )
    results = goea.run_study(study_genes)

    # 3. 转 DataFrame
    sig = [r for r in results if r.p_fdr_bh < alpha]
    if not sig:
        print("[富集] 无显著条目")
        return pd.DataFrame()
    df_out = pd.DataFrame([
        {
            "GO": r.GO,
            "name": r.name,
            "namespace": r.namespace,
            "p_uncorrected": r.p_uncorrected,
            "p_fdr_bh": r.p_fdr_bh,
            "depth": r.depth
        } for r in sig
    ])
    return df_out


# ========== 4. 主流程 ==========
def main():
    pdb_file = input("请输入 PDB 文件路径: ").strip().strip('"')
    print("1. 上传并等待 DeepFri 结果...")
    json_data = upload_and_get_result(pdb_file)

    print("2. 提取 GO 预测...")
    df_go = extract_go_predictions(json_data, score_threshold=0.5)
    print(df_go.head())

    print("3. 下载 GO 结构文件...")
    obo = "go-basic.obo"
    if not os.path.exists(obo):
        print("正在下载 go-basic.obo...")
        os.system(f"wget -q http://geneontology.org/ontology/{obo}")

    print("4. 执行 GO 富集分析...")
    df_enrich = go_enrichment(df_go)
    if not df_enrich.empty:
        csv_out = "deepfri_go_enrichment.csv"
        df_enrich.to_csv(csv_out, index=False)
        print(f"[富集] 显著条目已保存为 {csv_out}")
        print(df_enrich)
    else:
        print("[富集] 无显著条目")


if __name__ == "__main__":
    main()