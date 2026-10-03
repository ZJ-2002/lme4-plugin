# 中文注释：普通 LMM/GLMM 拟合节点（lme4 插件家族）。
# 输入构念：长表 TSV，一行一个观测；供者/动物/切片等分组因子是数据列。
# 规则来源：方案 §16.1——重复测量用连续时间混合模型，保留信号×时间项；
# 单独主效应只反映水平差异。统计口径：gaussian → lmerTest::lmer
# （Satterthwaite/Kenward ddf），其余 family → glmer Wald z。
# 失败处理：公式或列缺失 → 退出 2；不收敛/奇异拟合不中止，全文入
# lme4_model.txt，由研究者判读（解释边界：随机方差→0=不可识别，非无变异）。
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(data.table)
})

formula_str <- Sys.getenv("LME4_FORMULA")
family_str  <- Sys.getenv("LME4_FAMILY", "gaussian")
reml_flag   <- Sys.getenv("LME4_REML", "true") == "true"
ddf_method  <- Sys.getenv("LME4_DDF", "satterthwaite")
nagq        <- as.integer(Sys.getenv("LME4_NAGQ", "1"))
n_boot      <- as.integer(Sys.getenv("LME4_NBOOT", "0"))
seed        <- as.integer(Sys.getenv("LME4_SEED", "42"))

input_path  <- Sys.getenv("AUTONOMICS_INPUT0")
out_fixed   <- Sys.getenv("AUTONOMICS_OUTPUT0")
out_random  <- Sys.getenv("AUTONOMICS_OUTPUT1")
out_model   <- Sys.getenv("AUTONOMICS_OUTPUT2")

fail <- function(message) {
  message(message)
  quit(status = 2)
}

if (nzchar(formula_str) == FALSE) fail("LME4_FORMULA is empty")
if (!family_str %in% c("gaussian", "binomial", "poisson")) {
  fail(paste("unsupported family:", family_str))
}

frame <- fread(input_path, sep = "\t", na.strings = c("NA", ""))
parsed <- tryCatch(as.formula(formula_str), error = function(e) NULL)
if (is.null(parsed)) fail(paste("cannot parse formula:", formula_str))
needed <- all.vars(parsed)
missing_cols <- setdiff(needed, names(frame))
if (length(missing_cols) > 0) {
  fail(paste("columns missing from input:", paste(missing_cols, collapse = ", ")))
}

set.seed(seed)
fit_info <- list()
if (family_str == "gaussian") {
  fit <- tryCatch(
    lmerTest::lmer(parsed, data = frame, REML = reml_flag),
    error = function(e) fail(paste("lmer failed:", conditionMessage(e)))
  )
  # Kenward 需 pbkrtest；Satterthwaite 是 lmerTest 默认。
  if (ddf_method == "kenward") {
    fit <- lmerTest::as_lmerModLmerTest(fit)
    fit_info$ddf <- "kenward"
  } else {
    fit_info$ddf <- "satterthwaite"
  }
} else {
  fit <- tryCatch(
    glmer(parsed, data = frame, family = family_str, nAGQ = nagq,
          control = glmerControl(optimizer = "bobyqa")),
    error = function(e) fail(paste("glmer failed:", conditionMessage(e)))
  )
  fit_info$ddf <- "wald-z"
}

fixed <- broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE)
if (n_boot > 0 && family_str == "gaussian") {
  boot_fit <- bootMer(fit, FUN = fixef, nsim = n_boot, seed = seed, .progress = "none")
  boot_ci <- t(apply(boot_fit$t, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
  fixed$boot_low <- boot_ci[, 1]
  fixed$boot_high <- boot_ci[, 2]
  fit_info$bootstrap <- n_boot
}
fwrite(fixed, out_fixed, sep = "\t")

random <- broom.mixed::tidy(fit, effects = "ran_pars")
fwrite(random, out_random, sep = "\t")

convergence_notes <- character(0)
if (family_str == "gaussian" && length(fit@optinfo$conv$lme4$messages) > 0) {
  convergence_notes <- fit@optinfo$conv$lme4$messages
}
is_singular <- isSingular(fit)

model_text <- c(
  capture.output(summary(fit)),
  "",
  sprintf("formula: %s", formula_str),
  sprintf("family: %s | ddf: %s | REML: %s | nAGQ: %d | seed: %d",
          family_str, fit_info$ddf, reml_flag, nagq, seed),
  sprintf("n_obs: %d | singular_fit: %s", nobs(fit), is_singular),
  if (length(convergence_notes) > 0) paste("convergence:", convergence_notes) else "convergence: no lme4 messages",
  if (is_singular) "NOTE singular fit: a random-effect variance collapsed to ~0; that effect is not identifiable here - do not read as zero variation (plan §16.1)." else character(0)
)
writeLines(model_text, out_model)

if (is_singular || length(convergence_notes) > 0) {
  message("fit has singular/convergence notes; see lme4_model.txt")
}
