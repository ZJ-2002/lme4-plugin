# 中文注释：审计 F10 的回归锁定套件（2026-10-05）。
# 运行：Rscript tests/test_lme4_fit.R   （需 lme4/lmerTest/broom.mixed/data.table）
# 覆盖：
#   1. ddf=kenward 真正生效——固定效应自由度与 p 值在同一 R 会话内与
#      lmerTest 官方 summary(ddf="Kenward-Roger") 完全一致，且与默认
#      Satterthwaite 路径不同（KR=17 vs ST≈16.9997，sleepstudy 经典差异）；
#   2. lme4_model.txt 的机器可读 status 块：ok / singular / not-converged
#      三种病灶各自锁定，非收敛时输出仍写但 status 显式标注；
#   3. ddf 取值白名单：非法值退出 2。
# 口径：脚本 vs 官方实现同会话自洽对拍（官方实现自身比对，非独立金标准）。
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(data.table)
})

script <- file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = 1))))),
                    "scripts", "lme4_fit.R")
work <- tempfile("lme4_test_")
dir.create(work)
failures <- character(0)

check <- function(label, condition) {
  if (isTRUE(condition)) {
    cat("PASS", label, "\n")
  } else {
    failures <<- c(failures, label)
    cat("FAIL", label, "\n")
  }
}

run_node <- function(input_csv, formula_str, ddf = "satterthwaite", family = "gaussian") {
  out_fixed <- file.path(work, paste0("fixed_", ddf, "_", family, ".tsv"))
  out_random <- file.path(work, paste0("random_", ddf, "_", family, ".tsv"))
  out_model <- file.path(work, paste0("model_", ddf, "_", family, ".txt"))
  env <- paste(
    sprintf("LME4_FORMULA=%s", shQuote(formula_str)),
    sprintf("LME4_FAMILY=%s", shQuote(family)),
    "LME4_REML=true",
    sprintf("LME4_DDF=%s", shQuote(ddf)),
    "LME4_NAGQ=1", "LME4_NBOOT=0", "LME4_SEED=42",
    sprintf("AUTONOMICS_INPUT0=%s", shQuote(input_csv)),
    sprintf("AUTONOMICS_OUTPUT0=%s", shQuote(out_fixed)),
    sprintf("AUTONOMICS_OUTPUT1=%s", shQuote(out_random)),
    sprintf("AUTONOMICS_OUTPUT2=%s", shQuote(out_model)))
  stderr_file <- file.path(work, "stderr.txt")
  exit_code <- system2("Rscript", c(shQuote(script)), stdout = FALSE,
                       stderr = stderr_file, env = env)
  list(exit_code = exit_code, fixed = out_fixed, random = out_random,
       model = out_model, stderr = paste(readLines(stderr_file, warn = FALSE), collapse = "\n"))
}

write_input <- function(data_frame, name) {
  path <- file.path(work, name)
  fwrite(data_frame, path, sep = "\t")
  path
}

# ---- fixture 1：sleepstudy 类数据（官方数据集，KR/ST 自由度可区分）----
ok_csv <- write_input(as.data.frame(lme4::sleepstudy), "sleepstudy.tsv")
reference <- lmerTest::lmer(Reaction ~ Days + (Days | Subject), data = lme4::sleepstudy)
ref_kr <- summary(reference, ddf = "Kenward-Roger")$coefficients
ref_st <- summary(reference, ddf = "Satterthwaite")$coefficients
check("reference KR and Satterthwaite df differ",
      !isTRUE(all.equal(ref_kr[, "df"], ref_st[, "df"])))

run_kr <- run_node(ok_csv, "Reaction ~ Days + (Days | Subject)", ddf = "kenward")
run_st <- run_node(ok_csv, "Reaction ~ Days + (Days | Subject)", ddf = "satterthwaite")
check("kenward run exits 0", run_kr$exit_code == 0)
check("satterthwaite run exits 0", run_st$exit_code == 0)

kr <- fread(run_kr$fixed, sep = "\t")
st <- fread(run_st$fixed, sep = "\t")
check("fixed table has df/p/conf columns",
      all(c("df", "p.value", "conf.low", "conf.high") %in% names(kr)))
check("kenward df equals official KR summary",
      isTRUE(all.equal(kr$df, as.numeric(ref_kr[, "df"]), tolerance = 1e-10)))
check("kenward p equals official KR summary",
      isTRUE(all.equal(kr$p.value, as.numeric(ref_kr[, "Pr(>|t|)"]), tolerance = 1e-10)))
check("satterthwaite df equals official ST summary",
      isTRUE(all.equal(st$df, as.numeric(ref_st[, "df"]), tolerance = 1e-10)))
check("kenward df differs from satterthwaite run",
      !isTRUE(all.equal(kr$df, st$df)))
model_kr <- paste(readLines(run_kr$model, warn = FALSE), collapse = "\n")
check("model txt records applied Kenward-Roger",
      grepl("ddf: kenward \\(applied: Kenward-Roger\\)", model_kr))
check("ok run status is ok", grepl("status: ok", model_kr))
check("ok run converged true singular false",
      grepl("converged: true", model_kr) && grepl("singular_fit: false", model_kr))

# ---- fixture 2：奇异拟合（受试者层 x + (x|g)，斜率方差塌缩到 0）----
set.seed(7)
g <- rep(1:20, each = 5); x <- rep(rnorm(20), each = 5)
singular_csv <- write_input(data.frame(y = 0.5 * x + rnorm(100), x = x, g = g), "singular.tsv")
run_sg <- run_node(singular_csv, "y ~ x + (x | g)")
model_sg <- paste(readLines(run_sg$model, warn = FALSE), collapse = "\n")
check("singular run exits 0", run_sg$exit_code == 0)
check("singular status is singular",
      grepl("status: singular", model_sg) && !grepl("status: not-converged", model_sg))
check("singular flags converged true singular_fit true",
      grepl("converged: true", model_sg) && grepl("singular_fit: true", model_sg))
check("singular note explains identifiability",
      grepl("NOTE singular fit", model_sg))
check("singular run mentions status on stderr",
      grepl("fit status: singular", run_sg$stderr))

# ---- fixture 3：不收敛（Response 放大 1e9 → lmerTest eigenvalue 警告）----
nc_csv <- write_input(transform(lme4::sleepstudy, Reaction = Reaction * 1e9), "not_converged.tsv")
run_nc <- run_node(nc_csv, "Reaction ~ Days + (Days | Subject)")
model_nc <- paste(readLines(run_nc$model, warn = FALSE), collapse = "\n")
check("not-converged run exits 0 but is marked",
      run_nc$exit_code == 0 && grepl("status: not-converged", model_nc))
check("not-converged flags converged false",
      grepl("converged: false", model_nc))
check("not-converged keeps lme4_messages line",
      grepl("lme4_messages:", model_nc))
check("not-converged note present",
      grepl("NOTE not converged", model_nc))
check("not-converged run mentions status on stderr",
      grepl("fit status: not-converged", run_nc$stderr))
check("not-converged outputs still written",
      file.exists(run_nc$fixed) && nrow(fread(run_nc$fixed, sep = "\t")) > 0)

# ---- fixture 4：ddf 白名单 ----
run_bad <- run_node(ok_csv, "Reaction ~ Days + (Days | Subject)", ddf = "bogus")
check("invalid ddf exits 2", run_bad$exit_code == 2)

if (length(failures)) {
  cat("\nFAILED:", paste(failures, collapse = ", "), "\n")
  quit(status = 1)
}
cat("\nall checks passed\n")
