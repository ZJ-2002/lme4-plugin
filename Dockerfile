# R 混合模型镜像：lme4 插件家族专用。
# rocker/r-ver:4.4.1 自带编译工具链且默认仓库指向 Posit P3M 的 jammy 二进制源，
# lme4/lmerTest/pbkrtest/broom.mixed/emmeans/data.table 均以二进制安装，无需源码编译。
FROM docker.io/rocker/r-ver:4.4.1

RUN Rscript -e 'install.packages(c("lme4","lmerTest","pbkrtest","broom.mixed","emmeans","data.table"))'

# 中文注释：运行期只读根文件系统；随机种子由节点参数传入并在脚本内 set.seed。
