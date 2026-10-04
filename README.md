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

## 口径

- `lmerTest::lmer`（含 ddf）或 `glmer`（bobyqa）；`broom.mixed::tidy`
  统一输出。**ddf 参数真正传入 `summary(fit, ddf=)` 与
  `broom.mixed::tidy(..., ddf=)`**（2026-10-05 审计 F10 修复：旧实现
  ddf=kenward 只写 `fit_info$ddf` 标签、从未传参，KR 自由度实际从未生效）；
  `ddf` 取值白名单 satterthwaite|kenward，非法值退出 2。
- `lme4_model.txt` 末尾有机器可读 status 块（审计 F10 收敛门控）：
  `status: ok|singular|not-converged|not-converged+singular`，附
  `converged` / `singular_fit` / `optimizer_code` / `lme4_messages` 四行。
  判据 = optimizer 退出码 + `optinfo$conv$lme4$messages` + 拟合期
  warnings（lmerTest 的 "eigenvalue close to zero" 类提示只走 warnings()
  不进 optinfo，已一并捕获）+ `isSingular`；"boundary (singular) fit"
  归 singular 通道不算收敛失败。非收敛/奇异时三件输出仍写（退出 0、
  stderr 提示），但 status 显式标注，不得当干净 RUN_PASS 解读。
- `isSingular` 检出即写解释性注释进 model 文本：随机方差→0 是**模型不可辨识**，
  不是"变异为零"的证据（奇异拟合不得按显著结果解读）。
- Kenward-Roger 自由度代价高，默认 satterthwaite；KR 适合小样本
  终稿复核。

## 布局

```text
lme4/
├── manifest.toml
├── scripts/lme4_fit.R
├── Dockerfile
├── README.md
└── tests/test_lme4_fit.R   # R 测试（宿主跑，不进镜像）
```

镜像：`localhost/autonomics/lme4-r@sha256:e224e7fa972342afafac184edad050bd7aec48d8cc206ddf5be1cd30e9d7e18d`
（rocker/r-ver:4.4.1 + lme4/lmerTest/pbkrtest/broom.mixed/emmeans/data.table，P3M 二进制装包；
KR 依赖的 pbkrtest 镜像里已有，本次修复无需重建镜像）。

## 测试

```bash
Rscript tests/test_lme4_fit.R   # 23 checks all passed（2026-10-05，宿主 R 4.6.1）
```

覆盖：ddf=kenward 输出与同一 R 会话内官方 `lmerTest::lmer +
summary(ddf="Kenward-Roger")` 的 df/p 完全一致（KR=17 vs ST≈16.9997，
两路径确实不同）；奇异 fixture（受试者层 x + `(x|g)`）status=singular
且 converged=true；不收敛 fixture（Response×1e9 触发 lmerTest
eigenvalue 警告）status=not-converged、输出仍写；非法 ddf 退出 2。
口径：脚本与官方实现在同一 R 会话内自洽对拍（官方实现自身比对，
非独立金标准）。宿主依赖：`install.packages(c("lme4","lmerTest",
"pbkrtest","broom.mixed"), lib="~/R/library")`（2026-10-05 已装）。

## 冒烟记录

真实斜率 0.30 的合成数据：估计 0.2989；受试者 SD 4.80（真值 5）；
bootMer CI 覆盖真值。
