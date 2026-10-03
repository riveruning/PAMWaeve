"""Conservative presentation of existing scores; no learned or combined score."""
from __future__ import annotations

import csv
import html
import json
from pathlib import Path


LIMITATIONS = [
    "分数是候选 PAM 的模型兼容度或侧翼证据，不是切割活性、准确率或成功概率。",
    "只在用户提供的候选集合内排序，不能据此声称找到了全部或最优 PAM。",
    "NRG 等 IUPAC 候选代表序列集合；高分不表示集合内每个具体 PAM 的活性相同。",
    "两路证据可能相关；本报告不计算融合总分，不声称融合优于现有方法。",
    "未执行训练集泄漏审计，训练暴露状态为未知；本次运行不是独立泛化评估。",
    "需由用户确认输入为完整 Cas9、Cas-array 配对及 spacer 方向；软件不自动鉴定。",
    "spacer 匹配仅支持替换错配；target contig 数不等于独立病毒数量。",
    "证据状态阈值（40、60）仅用于未校准的人工复核分流，需实验验证。",
]


def evidence_status(protein, spacer, *, orientation_warning=False):
    if spacer is None:
        return "仅蛋白预测；未提供 spacer 证据"
    if orientation_warning:
        return "方向警告；暂停解释 spacer 支持"
    if float(spacer.get("effective_spacer_hits", 0)) < 2:
        return "spacer 证据不足（少于 2 条有效 spacer）"
    p = protein.get("specificity_adjusted_score")
    s = spacer.get("spacer_evidence_score")
    if p is None or s is None:
        return "证据不足"
    p, s = float(p), float(s)
    if (p > 60 and s < 40) or (s > 60 and p < 40):
        return "证据冲突；需人工复核"
    if p > 60 and s > 60:
        return "两路均支持（非独立验证）"
    return "中性或未获两路共同支持"


def make_report(protein, spacer=None, *, mode, provenance):
    spacer_rows = {r["candidate_pam"]: r for r in (spacer or {}).get("candidate_scores", [])}
    rows = []
    for p in protein["scores"]:
        s = spacer_rows.get(p["candidate_pam"])
        rows.append({
            "candidate_pam": p["candidate_pam"],
            "protein_score": p["specificity_adjusted_score"],
            "rank_within_length": p.get("rank_within_length"),
            "spacer_score": s.get("spacer_evidence_score") if s else None,
            "effective_spacers": s.get("effective_spacer_hits") if s else None,
            "support_spacers": s.get("support_spacer_count") if s else None,
            "support_target_contigs": s.get("support_target_count") if s else None,
            "status": evidence_status(p, s, orientation_warning=(spacer or {}).get("orientation_warning", False)),
        })
    # "no clear PAM" is a statement about the model, not about the row count: a
    # partial consensus produces no candidate rows while still having determined
    # some positions, so the two must not be conflated.
    no_clear = bool(protein.get("no_clear_pam"))
    partial = [p.get("predicted_pam") for p in protein["proteins"]
               if p.get("partial_consensus")]
    return {"format_version": 1, "mode": mode,
            "model": protein["model"], "proteins": protein["proteins"],
            "training_exposure": "unknown", "provenance": provenance,
            "limitations": LIMITATIONS, "candidates": rows,
            "no_clear_pam": no_clear,
            "no_candidate_rows": not rows,
            "partial_consensus": partial,
            "candidate_scope_note": protein.get("candidate_scope_note"),
            "spacer_summary": spacer}


def write_report(report, outdir):
    outdir = Path(outdir)
    (outdir / "report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    rows = report["candidates"]
    header = ["candidate_pam", "protein_score", "rank_within_length", "spacer_score",
              "effective_spacers", "support_spacers", "support_target_contigs", "status"]
    with (outdir / "candidates.tsv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    esc = lambda x: html.escape("未提供" if x is None else str(x))
    headings = ["候选 PAM", "蛋白兼容度", "同长度排名", "spacer 分数", "有效 spacer", "支持 spacer", "支持 contig", "状态"]
    table = '<tr>' + ''.join('<th>'+h+'</th>' for h in headings) + '</tr>'
    table += ''.join('<tr>'+''.join('<td>'+esc(row.get(k))+'</td>' for k in header)+'</tr>' for row in rows)
    if not rows:
        table += ('<tr><td colspan="8">没有候选行：模型没有给出可评分的明确 PAM。'
                  '这不是生物学 PAM 结论，也不代表候选筛选成功。</td></tr>')
    demo = "合成演示：没有运行蛋白模型；不是生物学测试结果。" if report["mode"] == "synthetic_demo" else "真实输入推理：未经实验验证，不是独立 benchmark。"
    consensus = '; '.join(str(p.get('predicted_pam') or '无明确共识') for p in report['proteins'])
    no_clear = report.get("no_clear_pam")
    conclusion = ('<aside class="warn">未得到明确 PAM。请勿把本次结果当作该 Cas9 的生物学 PAM，'
                  '也不要把占位符当作候选筛选成功。可在明确候选集合下重跑。'
                  '</aside>' if no_clear else '')
    if not no_clear and report.get("partial_consensus"):
        detail = "; ".join(
            f"{p.get('predicted_pam')}：已判定位置 {p.get('determined_positions')}"
            for p in report["proteins"] if p.get("partial_consensus")
        )
        conclusion = ('<aside class="warn">共识中含有“未判定”位置（N），它描述的是一组 PAM，'
                      '而不是单一序列。只有已判定位置是模型的实际调用，N 不是碱基调用，'
                      '不得按某个碱基解读。' + detail + '。可在指定候选集合后重跑。'
                      '</aside>')
    scope = report.get("candidate_scope_note")
    page = ('<!doctype html><html lang="zh-CN"><meta charset="utf-8">'
            '<title>PAMPRIDICT 候选 PAM 报告</title><style>'
            'body{font-family:system-ui,sans-serif;margin:32px;color:#182938}table{border-collapse:collapse}'
            'td,th{border:1px solid #ccd;padding:10px;text-align:left}aside{padding:16px;background:#fff2cf}'
            '.warn{background:#ffe3e3}pre{white-space:pre-wrap;overflow-wrap:anywhere}</style>'
            '<h1>候选 PAM 与证据报告</h1>'
            '<aside>'+demo+'</aside>'+conclusion
            +'<p>模型：'+esc(report['model'])+'；共识：'+esc(consensus)+'</p>'
            +(('<p>候选范围：'+esc(scope)+'</p>') if scope else '')
            +'<table>'+table+'</table><h2>适用边界</h2><ul>'
            +''.join('<li>'+esc(v)+'</li>' for v in report['limitations'])+'</ul>'
            '<h2>输入与运行记录</h2><pre>'+esc(json.dumps(report['provenance'],ensure_ascii=False,indent=2))+'</pre></html>')
    (outdir / "report.html").write_text(page, encoding="utf-8")
