# lme4

普通 LMM/GLMM 拟合。补齐引擎缺口：registry 只有 GEE 和 statsmodels
路径，缺 lme4 语义的随机截距/斜率模型（主计划 §16.1 供者随机效应、
§15.5-4 切片嵌套供者的嵌套结构）。

## 节点

| kind | 输入 | 输出 |
|---|---|---|
| `lme4_fit` | 长 TSV | `lme4_fixed_effects.tsv`（term/estimate/std.error/df/statistic/p.value/conf.low/conf.high）、`lme4_random_effects.tsv`（VarCorr）、`lme4_model.txt` |

参数：`formula`（必填，lme4 语法）、`family`（gaussian|binomial|poisson）、
`reml`（默认 true，glmer 忽略）、`ddf`（satterthwaite|kenward）、
`nagq`、`n_boot`（>0 时 bootMer 置信区间）、`seed`。

## 布局

```text
lme4/
├── manifest.toml
├── scripts/lme4_fit.R
├── Dockerfile
└── README.md
```

镜像：`localhost/autonomics/lme4-r@sha256:e224e7fa972342afafac184edad050bd7aec48d8cc206ddf5be1cd30e9d7e18d`
（rocker/r-ver:4.4.1 + lme4/lmerTest/pbkrtest/broom.mixed/emmeans/data.table，P3M 二进制装包）。

## 口径

- `lmerTest::lmer`（含 ddf）或 `glmer`（bobyqa）；`broom.mixed::tidy`
  统一输出。
- `isSingular` 检出即写解释性注释进 log：随机方差→0 是**模型不可辨识**，
  不是"变异为零"的证据（奇异拟合不得按显著结果解读）。
- Kenward-Roger 自由度代价高，默认 satterthwaite；KR 适合小样本
  终稿复核。

## 冒烟记录

真实斜率 0.30 的合成数据：估计 0.2989；受试者 SD 4.80（真值 5）；
bootMer CI 覆盖真值。
