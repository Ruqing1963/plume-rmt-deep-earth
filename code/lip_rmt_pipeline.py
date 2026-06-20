#!/usr/bin/env python3
"""
═══════════════════════════════════════════════════════════════════════════════
  地幔柱 / 大火成岩省 (LIPs) 喷发节律 RMT 分析流水线
  Mantle Plume / LIP Eruption Rhythm — Random Matrix Theory Pipeline
  ─────────────────────────────────────────────────────────────────
  核心战术: 对抗谱叠加定理 (吸取 Paper 2/3 教训)
    实验A 全球混合态  → 预期退化为 Poisson
    实验B 单源聚集态  → 预期回归 GOE 互斥
  
  Author: Ruqing Chen, GUT Geoservice Inc., Montreal
═══════════════════════════════════════════════════════════════════════════════
"""
import numpy as np
import pandas as pd
from scipy import stats
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid
import warnings
warnings.filterwarnings('ignore')

# ─── 理论分布 ────────────────────────────────────────────────────────────
def wigner_goe(s): return (np.pi/2)*s*np.exp(-np.pi*s**2/4)
def wigner_gue(s): return (32/np.pi**2)*s**2*np.exp(-4*s**2/np.pi)
def poisson_pdf(s): return np.exp(-s)
def make_cdf(f, mx=8, n=10000):
    s=np.linspace(0,mx,n); c=cumulative_trapezoid(f(s),s,initial=0); c/=c[-1]
    return interp1d(s,c,bounds_error=False,fill_value=(0,1))
POI = lambda x: 1-np.exp(-x)
GOE_CDF = make_cdf(wigner_goe)
GUE_CDF = make_cdf(wigner_gue)

# ─── 统计量 ──────────────────────────────────────────────────────────────
def spacing_ratio(sp):
    """间距比 ⟨r⟩: Poisson=0.386, GOE=0.536, GUE=0.603"""
    if len(sp) < 3: return np.nan, np.nan
    r = np.minimum(sp[:-1],sp[1:]) / np.maximum(sp[:-1],sp[1:])
    return np.mean(r), np.std(r)/np.sqrt(len(r))

def beta_from_r(r):
    if np.isnan(r): return np.nan
    if r<=0.386: return 0.0
    elif r<=0.536: return (r-0.386)/0.15
    elif r<=0.603: return 1+(r-0.536)/0.067
    else: return min(2+(r-0.603)/0.1, 3)

def cv(iv): return np.std(iv)/np.mean(iv)

def unfold_intervals(ages):
    """
    展开(Unfolding): 把绝对年龄序列转为归一化间距。
    LIP 喷发率可能随时间变化(地球冷却), 用局部平均率展开。
    简化: 按时间排序取间隔, 再除以全局均值(因样本少, 避免过度展开)。
    """
    a = np.sort(np.asarray(ages))
    iv = np.diff(a)
    iv = iv[iv > 0]
    if len(iv) < 3: return None
    return iv / np.mean(iv)

def rmt_test(ages, label):
    s = unfold_intervals(ages)
    if s is None or len(s) < 4:
        return None
    iv_raw = np.diff(np.sort(ages)); iv_raw=iv_raw[iv_raw>0]
    r, re = spacing_ratio(s)
    ksp, pp = stats.kstest(s, POI)
    kso, po = stats.kstest(s, GOE_CDF)
    ksu, pu = stats.kstest(s, GUE_CDF)
    best = min([('Poisson',ksp),('GOE',kso),('GUE',ksu)], key=lambda x:x[1])[0]
    return {
        'label': label, 'n_events': len(ages), 'n_intervals': len(s),
        'r': r, 'r_err': re, 'beta': beta_from_r(r), 'cv': cv(iv_raw),
        'mean_interval_Ma': np.mean(iv_raw),
        'p_poisson': pp, 'p_goe': po, 'p_gue': pu, 'best': best,
        's': s.tolist()
    }

# ═══════════════════════════════════════════════════════════════════════════
# 主流程
# ═══════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print("="*72)
    print("  地幔柱/LIPs 喷发节律 RMT 分析 — 直插地核的探针")
    print("="*72)

    df = pd.read_csv('/home/claude/lip_database.csv', comment='#')
    df = df.dropna(subset=['age_ma']).sort_values('age_ma').reset_index(drop=True)
    print(f"\n  LIP 数据库: {len(df)} 个事件, {df['age_ma'].min():.0f}–{df['age_ma'].max():.0f} Ma")
    print(f"  来源: Ernst (2014) LIP record + Torsvik LLSVP 重建")

    results = {}

    # ─── 实验 A: 全球混合态 ───
    print("\n" + "="*72)
    print("  实验 A: 全球混合态 (所有 LIPs 不加区分)")
    print("  预期: 谱叠加定理 → 退化为 Poisson")
    print("="*72)
    resA = rmt_test(df['age_ma'].values, 'Global-Mixed')
    results['A_global'] = resA
    print(f"\n  n={resA['n_intervals']} 间隔  平均间隔={resA['mean_interval_Ma']:.1f} Ma")
    print(f"  ⟨r⟩ = {resA['r']:.3f} ± {resA['r_err']:.3f}   β̂ = {resA['beta']:.2f}   CV = {resA['cv']:.2f}")
    print(f"  KS p: Poisson={resA['p_poisson']:.3f}  GOE={resA['p_goe']:.3f}  GUE={resA['p_gue']:.3f}")
    print(f"  → 最优拟合: {resA['best']}")

    # ─── 实验 B: 单源聚集态 (LLSVP 分离) ───
    print("\n" + "="*72)
    print("  实验 B: 单源聚集态 (按 LLSVP 超级地幔柱区分离)")
    print("  预期: 孤立热源 → 回归 GOE 互斥")
    print("="*72)
    for domain in ['African', 'Pacific']:
        sub = df[df['llsvp_domain']==domain]
        res = rmt_test(sub['age_ma'].values, f'{domain}-LLSVP')
        results[f'B_{domain}'] = res
        if res is None:
            print(f"\n  {domain} 域: 样本不足"); continue
        print(f"\n  ─── {domain} 超级地幔柱域 (Tuzo/Jason) ───")
        print(f"  n={res['n_intervals']} 间隔  平均间隔={res['mean_interval_Ma']:.1f} Ma")
        print(f"  ⟨r⟩ = {res['r']:.3f} ± {res['r_err']:.3f}   β̂ = {res['beta']:.2f}   CV = {res['cv']:.2f}")
        print(f"  KS p: Poisson={res['p_poisson']:.3f}  GOE={res['p_goe']:.3f}  GUE={res['p_gue']:.3f}")
        print(f"  → 最优拟合: {res['best']}")

    # ─── 裁决 ───
    print("\n" + "="*72)
    print("  谱叠加对抗实验 — 裁决")
    print("="*72)
    print(f"\n  {'实验':22s}{'n':>4s}{'⟨r⟩':>8s}{'β̂':>7s}{'best':>9s}")
    print("  " + "-"*52)
    for key in ['A_global','B_African','B_Pacific']:
        r = results.get(key)
        if r:
            print(f"  {r['label']:22s}{r['n_intervals']:4d}{r['r']:8.3f}{r['beta']:7.2f}{r['best']:>9s}")

    rA = results['A_global']['r']
    rAf = results['B_African']['r'] if results.get('B_African') else np.nan
    rPa = results['B_Pacific']['r'] if results.get('B_Pacific') else np.nan
    print(f"\n  全球混合 ⟨r⟩ = {rA:.3f}")
    print(f"  单源平均 ⟨r⟩ = {np.nanmean([rAf,rPa]):.3f}")
    if np.nanmean([rAf,rPa]) > rA + 0.03:
        print(f"\n  ✓ 假说支持: 单源互斥强于全球混合 → 谱叠加确实抹平了热阴影信号!")
    else:
        print(f"\n  ✗ 假说未获支持: 单源未显著回归互斥")

    import json
    with open('/home/claude/lip_rmt_results.json','w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n  ✓ 已保存 lip_rmt_results.json")
