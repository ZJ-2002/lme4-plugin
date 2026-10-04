# 中文注释：普通 LMM/GLMM 拟合节点（lme4 插件家族）。
# 输入构念：长表 TSV，一行一个观测；供者/动物/切片等分组因子是数据列。
# 规则来源：方案 §16.1——重复测量用连续时间混合模型，保留信号×时间项；
# 单独主效应只反映水平差异。统计口径：gaussian → lmerTest::lmer
# （Satterthwaite/Kenward ddf，ddf 参数真正传入 summary 与 broom.mixed::tidy
# ——2026-10-05 审计 F10：旧实现 kenward 只写标签未传参），其余 family →
# glmer Wald z。ddf 取值白名单：satterthwaite|kenward（小样本终稿复核用 KR）。
# 失败处理：公式或列缺失 → 退出 2；不收敛/奇异拟合不中止，输出仍写，
# 但 lme4_model.txt 末尾的机器可读 status 块显式标记
# （status: ok|singular|not-converged|not-converged+singular，含 converged/
# singular_fit/optimizer_code/lme4_messages），不得当干净 RUN_PASS 解读
# （解释边界：随机方差→0=不可识别，非无变异）。
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(data.table)
})

formula_str <- Sys.getenv("LME4_FORMULA")
family_str  <- Sys.getenv("LME4_FAMILY", "gaussian")
reml_flag   <- Sys.getenv("LME4_REML", "true") == "true"
ddf_method  <- tolower(Sys.getenv("LME4_DDF", "satterthwaite"))
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
if (!ddf_method %in% c("satterthwaite", "kenward")) {
  fail(paste("unsupported ddf:", ddf_method))
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
# 拟合期警告捕获：lmerTest/lme4 的部分不可信拟合提示（如 "eigenvalue close
# to zero"）走 warnings() 而不进 optinfo，收敛门控必须一并收口（审计 F10）。
fit_warn_msgs <- character(0)
fit_with_warnings <- function(fit_call, label) {
  withCallingHandlers(
    tryCatch(fit_call, error = function(e) fail(paste(label, "failed:", conditionMessage(e)))),
    warning = function(w) {
      fit_warn_msgs <<- c(fit_warn_msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
}
if (family_str == "gaussian") {
  # lmerTest::lmer 本身返回 lmerModLmerTest；summary/tidy 的 ddf= 参数
  # 是 Kenward-Roger/Satterthwaite 真正生效的唯一通道（审计 F10）。
  ddf_arg <- if (ddf_method == "kenward") "Kenward-Roger" else "Satterthwaite"
  fit <- fit_with_warnings(lmerTest::lmer(parsed, data = frame, REML = reml_flag), "lmer")
  fit_info$ddf <- ddf_arg
} else {
  ddf_arg <- NA_character_
  fit <- fit_with_warnings(
    glmer(parsed, data = frame, family = family_str, nAGQ = nagq,
          control = glmerControl(optimizer = "bobyqa")),
    "glmer")
  fit_info$ddf <- "wald-z"
}

fixed <- if (family_str == "gaussian") {
  broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE, ddf = ddf_arg)
} else {
  broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE)
}
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

# 收敛门控结构化（审计 F10）：optimizer 退出码 + lme4 消息 + 拟合期警告 +
# isSingular。"boundary (singular) fit" 是奇异性提示不是收敛失败，归 singular 通道。
conv_code <- tryCatch(as.integer(fit@optinfo$conv$opt), error = function(e) NA_integer_)
conv_msgs <- fit@optinfo$conv$lme4$messages
if (is.null(conv_msgs)) conv_msgs <- character(0)
conv_warn_msgs <- c(conv_msgs[!grepl("boundary \\(singular\\) fit", conv_msgs)],
                    fit_warn_msgs[!grepl("boundary \\(singular\\) fit", fit_warn_msgs)])
singular <- isSingular(fit)
converged <- identical(conv_code, 0L) && length(conv_warn_msgs) == 0L

status <- if (!converged && singular) {
  "not-converged+singular"
} else if (!converged) {
  "not-converged"
} else if (singular) {
  "singular"
} else {
  "ok"
}

fit_summary <- if (family_str == "gaussian") {
  summary(fit, ddf = ddf_arg)
} else {
  summary(fit)
}

model_text <- c(
  capture.output(fit_summary),
  "",
  sprintf("formula: %s", formula_str),
  sprintf("family: %s | ddf: %s (applied: %s) | REML: %s | nAGQ: %d | seed: %d",
          family_str, ddf_method, fit_info$ddf, reml_flag, nagq, seed),
  sprintf("n_obs: %d", nobs(fit)),
  "## status (machine-readable; audit F10)",
  sprintf("status: %s", status),
  sprintf("converged: %s", tolower(as.character(converged))),
  sprintf("singular_fit: %s", tolower(as.character(singular))),
  sprintf("optimizer_code: %s", ifelse(is.na(conv_code), "NA", as.character(conv_code))),
  if (length(conv_msgs) > 0) paste("lme4_messages:", paste(conv_msgs, collapse = " | ")) else "lme4_messages: none",
  if (singular) "NOTE singular fit: a random-effect variance collapsed to ~0; that effect is not identifiable here - do not read as zero variation (plan §16.1)." else character(0),
  if (!converged) "NOTE not converged: optimizer did not reach a trustworthy optimum; treat estimates as exploratory, re-fit with adjusted control/data before inference (plan §16.1)." else character(0)
)
writeLines(model_text, out_model)

if (!converged || singular) {
  message(sprintf("fit status: %s; outputs written but not a clean pass; see lme4_model.txt", status))
}
