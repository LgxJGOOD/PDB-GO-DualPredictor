一体化 PDB 双通道 GO 预测 & 精确+语义相似性对比  
============================================================

**一句话描述**  
上传一个 PDB 文件，即可同时获得 **InterProScan** 与 **DeepFRI** 两套 GO 注释，并在“精确匹配 + 语义相似性”两个维度自动给出可比对的完整报告，无需额外安装 Wang 等外部模块，仅依赖 `goatools` 的 DAG 最短路径。

---

目录  
----  
1. 功能亮点  
2. 安装依赖  
3. 一键运行  
4. 结果解读  
5. 参数速查  
6. 典型场景  
7. 注意事项 & 限速提示  
8. 常见问题  
9. 更新日志  

---

1. 功能亮点  
-------------  
| 特性 | 说明 |
|---|---|
| **双通道预测** | 序列通道：InterProScan；结构通道：DeepFRI |
| **零人工比对** | 自动计算 Jaccard、高置信交集、语义相似度 |
| **纯 Python** | 仅依赖 `goatools` 的 OBO 解析，无需额外 C++ 模块 |
| **最短路径语义相似度** | 基于 GO DAG 自写 BFS，归一化得分 0–1 |
| **离线友好** | 只要提前下载 `go-basic.obo` 即可本地算语义距离 |
| **一行代码** | `analyze_pdb("xxx.pdb")` 返回全部指标 |

---

2. 安装依赖  
-------------  
```bash
# 1. 创建虚拟环境（推荐）
python -m venv venv
source venv/bin/activate   # Windows 用 venv\Scripts\activate

# 2. 安装 Python 包
pip install biopython requests goatools

# 3. 下载 GO 本体（仅需一次）
wget http://purl.obolibrary.org/obo/go/go-basic.obo
# 放在项目根目录，或与 main.py 同路径即可
```

---

3. 一键运行  
-------------  
```python
from main import analyze_pdb

report = analyze_pdb("1aki.pdb", email="you@example.com")
print(report)
```

控制台会实时提示  
```
InterProScan running … job_id=iprscan-R20251031-xxx  
DeepFRI uploading … task_name=1aki_xxx  
…
```

运行完毕返回一个 `dict`，可直接 `json.dump` 落盘。

---

4. 结果解读（返回字段速览）  
----------------------------  
| 字段 | 类型 | 含义 |
|---|---|---|
| `interpro_count` | int | InterProScan 注释到的 GO 数 |
| `deepfri_count` | int | DeepFRI 预测到的 GO 数 |
| `common` | set[str] | 两套注释完全一致的 GO ID |
| `common_count` | int | 精确匹配数量 |
| `high_confidence_common` | set[str] | 在 common 中 DeepFRI 得分 ≥0.8 的子集 |
| `jaccard` | float | 精确匹配 Jaccard 系数 |
| `interpro_only` / `deepfri_only` | set[str] | 各自独有 GO |
| `semantic_avg_sim` | float | 基于最短路径的语义平均相似度 |
| `semantic_pairs` | list[tuple] | 相似度≥阈值的最佳匹配对 (GO1, GO2, score) |
| `semantic_pairs_count` | int | 上述匹配对数量 |

---

5. 参数速查  
-----------  
```python
analyze_pdb(
    pdb_path: str,               # 必填
    email: str = "xxx@qq.com",   # 用于 InterProScan 提交流队列
    sem_sim_threshold: float = 0.7  # 语义匹配阈值，范围 0–1
)
```

---

6. 典型场景  
-----------  
1. **新解析的 PDB 没有功能注释**  
   → 直接跑本脚本，秒得两套 GO 并自动交叉验证。  
2. **文章审稿人要求“结构预测与实验注释一致性”**  
   → 给出 Jaccard + 语义相似度双重指标，一步到位。  
3. **批量评估 AlphaFold2 模型**  
   → 循环调用 `analyze_pdb`，把 `semantic_avg_sim` 画分布图即可。

---

7. 注意事项 & 限速提示  
----------------------  
| 项目 | 说明 |
|---|---|
| **InterProScan** | 免费队列 20–60 min/job，单 IP 并发≤3；请合理间隔。 |
| **DeepFRI** | 官方 beta 接口未公开速率上限，建议 `time.sleep(2)` 以上。 |
| **网络** | 国内用户若无法访问 `ebi.ac.uk` 或 `flatironinstitute.org`，需配置代理。 |
| **OBO 版本** | 不同版本 GO 可能导致最短路径距离略有差异；建议固定版本并注明文章方法。 |

---

8. 常见问题（FAQ）  
----------------  
**Q1: 报错 `RuntimeError: InterProScan failed`**  
A: 检查序列是否<10 aa 或含非法字符；也可先手动在 [https://www.ebi.ac.uk/interpro/](https://www.ebi.ac.uk/interpro/) 提交确认。  

**Q2: `go-basic.obo` 放在哪里？**  
A: 默认读取当前工作目录下的 `go-basic.obo`；也可修改脚本常量 `GO_OBO`。  

**Q3: 想换语义相似度算法？**  
A: 只需替换 `_semantic_pairs` 函数即可，其余模块不动。  

**Q4: 需要批量跑几百个 PDB？**  
A: 建议把 `INTERPRO_WAIT` 调大（≥30 s）并使用缓存，避免被 EBI 封 IP。  

---

9. 更新日志  
-----------  
| 日期 | 版本 | 说明 |
|---|---|---|
| 2025-10-31 | v1.0 | 初版，支持双通道预测 + 最短路径语义相似度 |

---
