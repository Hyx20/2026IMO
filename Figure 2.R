library(tidyverse)
library(data.table)
library(ggplot2)
library(ggalluvial)
library(ggsci)
# library(ggvenn)

####F2.1 adonis=====
adonis_df <- fread("data/adonis_results.csv")

var_map <- data.frame(
  Variable = c("educagroup","age","gender","income1","marrygroup",
               "smoke","drinkgroup3","BMI","occupation_g","citytype",
               "crop","fruit","vegetable","meat",
               "metawaist1","metaDM2","TG2","HDL2","metaHP2",
               "PM25","NO2","TEMP",
               "vil_name","Town_Name","District"),
  Label = c("Education","Age","Gender","Income","Marital status",
            "Smoking","Drinking","BMI","Occupation","Urban/Rural",
            "Crop intake","Fruit intake","Vegetable intake","Meat intake",
            "Abdominal obesity","High FBG","High TG","Low HDL","High BP",
            "PM2.5","NO2","Temperature",
            "Village","Town","District"),
  Group = c(
    rep("Demographic", 5),
    rep("Lifestyle", 5),
    rep("Diet", 4),
    rep("Health status", 5),
    rep("Environmental exposure", 3),
    rep("Geographical factor", 3)
  ),
  stringsAsFactors = FALSE
)

adonis_df <- adonis_df %>%
  left_join(var_map, by = "Variable")%>%
  arrange(Group, desc(R2))

F2.1 <- ggplot(adonis_df, aes(x = reorder(Label, -R2), 
                              y = R2, 
                              fill = Group)) +
  geom_col(width = 0.7) +
  theme_bw(base_size = 14) +
  scale_fill_npg(alpha = 0.9) +   # 每个分组一个颜色
  theme(
    axis.text.y = element_text(face = "bold", size = 12, colour = "black"),
    axis.text.x = element_text(face = "bold", size = 12, colour = "black", angle = 45, hjust = 1),
    axis.title  = element_text(face = "bold", size = 14, colour = "black"),
    legend.title = element_text(face = "bold", size = 13),
    legend.text  = element_text(size = 11),
    legend.position = "right"
  ) +
  labs(x = "", 
       y = "Explained variance (R²)", 
       fill = "Group", 
       title = "")

F2.1

# ggsave("Figure 2/Figure2a.pdf", F2.1, width = 12, height = 6)

####F2.2 关联热图=====
#load data
m_genues <- fread("data/lefse_genus_all.csv")
fc <- fread("data/fc_all.csv")
otu <- fread("data/otu.csv")

meta <- fread("data/umap.csv")

meta$bmi_g <- cut(meta$BMI,
                  breaks = c(-Inf, 18.5, 24, 28, Inf),
                  labels = c("Underweight", "Normal", "Overweight", "Obese"))

table(meta$bmi_g)

# 处理OTU表，提取分类信息
tax <- otu %>%
  select(Kingdom, Phylum, Class, Order, Family, Genus) %>%
  distinct()

tax$otu_id <- paste0("OTU",1:nrow(tax))
tax <- tax %>% column_to_rownames("otu_id")

tax <- tax %>%
  mutate(
    Genus = ifelse(duplicated(Genus) | duplicated(Genus, fromLast = TRUE),
                   paste0(Family, "_", Genus),
                   Genus)
  )

#计算相对丰度
otu_rel <- otu %>%
  select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus) %>%
  as.data.frame()

otu_rel <- sweep(otu_rel, 2, colSums(otu_rel), FUN = "/")

otu_rel <- t(otu_rel) %>% as.data.frame()
colnames(otu_rel) <- row.names(tax)

###并集
m_genues <- unique(m_genues$Genus)

fc$Genus <- paste0("g__", fc$Genus)
fc_genus    <- unique(fc$Genus)



# 并集
genus_union <- union(m_genues, fc_genus)

tax_filtered <- tax[tax$Genus %in% genus_union, ]

otu_rel_filtered <- otu_rel[, colnames(otu_rel) %in% rownames(tax_filtered)]
colnames(otu_rel_filtered) <- tax_filtered$Genus[match(colnames(otu_rel_filtered), rownames(tax_filtered))]

otu_rel_filtered <- otu_rel_filtered %>% as.data.frame()%>% rownames_to_column("SampleID")

meta <- meta %>% left_join(otu_rel_filtered, by = "SampleID")  

env_factors <- c("PM25","NO2","TEMP")

meta <- as.data.frame(meta)
iqr_vals <- sapply(meta[env_factors], IQR, na.rm = TRUE)  # 计算 IQR

for (v in env_factors) {
  meta[[v]] <- meta[[v]] / iqr_vals[[v]]
}

####model==
covars <- c("citytype","agegroup","gender","educagroup",
            "marrygroup","occupation_g","income1","fuelg",
            "smoke","drinkgroup3","total_ng","crop","vegetable","fruit","meat","season")
res <- data.frame()
for (i in env_factors) {
  for (g in tax_filtered$Genus) {
    
    fml <- as.formula(
      paste0("log(`", g, "` + 1e-6) ~ ", paste(c(i, covars), collapse = " + "))
    )
    
    model <- lm(fml, data = meta)
    # 提取系数和p值
    
    coefficients <- summary(model)$coefficients[i, 1]
    p_value <- summary(model)$coefficients[i, 4]
    
    # 将结果存储到数据框中
    res <- rbind(res, data.frame(
      Genus = g,
      EnvFactor = i,
      Coefficient = coefficients,
      PValue = p_value
    ))
  }
}

#FDR矫正
res$FDR <- p.adjust(res$PValue, method = "fdr")

res$log10_FDR <- -log10(res$FDR)*sign(res$Coefficient)

#标记显著生成 +/-
plot_df <- res %>%
  mutate(
    sig       = FDR < 0.05,
    sym       = case_when(
      sig  & Coefficient > 0 ~ "+",
      sig  & Coefficient < 0 ~ "-",
      TRUE ~ ""
    ),
    EnvFactor = factor(EnvFactor, levels = c("NO2","PM25","TEMP","RH"))
  )


plot_sig <- plot_df %>%
  group_by(Genus) %>%
  filter(any(sig, na.rm = TRUE)) %>%
  ungroup()


wide_sig <- plot_sig %>%
  mutate(Coeff_sig = ifelse(sig, Coefficient, 0)) %>%
  select(Genus, EnvFactor, Coeff_sig) %>%
  tidyr::pivot_wider(names_from = EnvFactor, values_from = Coeff_sig, values_fill = 0)

mat <- as.matrix(wide_sig[,-1])
rownames(mat) <- wide_sig$Genus

#判定每个属的“总体方向”：显著系数求和 >0 记为正，<0 记为负，=0 视为“混合”
n_pos <- rowSums(mat >  0, na.rm = TRUE)  # 显著后保留的矩阵 mat，非显著已置0
n_neg <- rowSums(mat <  0, na.rm = TRUE)

grp <- ifelse(n_pos > 0 & n_neg == 0, "Positive",
              ifelse(n_neg > 0 & n_pos == 0, "Negative", "Mixed"))
names(grp) <- rownames(mat)


F2.2 <- ggplot(
  plot_sig,
  aes(x = EnvFactor, y = Genus, fill = Coefficient)
) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = sym), size = 4.2, fontface = "bold") +
  scale_fill_gradient2(
    name = "Coefficient",
    low = "#3B4CC0", mid = "white", high = "#B40426",
    limits = c(-1.5, 1), midpoint = 0,
    breaks = c(-1.5,-1,-0.5,0,0.5,1)
  ) +
  scale_x_discrete(labels = c("NO2"="NO2","PM25"="PM2.5","TEMP"="Temperature","RH"="Relative Humidity"),
                   guide = guide_axis(angle = 45)) +
  labs(x = "Environmental factor", y = NULL) +
  coord_fixed(ratio = 0.25, clip = "off") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid   = element_blank(),
    axis.ticks   = element_blank(),
    axis.text.x  = element_text(color = "black", size = 10,hjust = 1, vjust = 1),
    axis.text.y  = element_text(color = "black", size = 9),
    legend.title = element_text(face = "bold", size = 13),
    legend.text  = element_text(size = 12),
    plot.margin  = margin(5, 8, 5, 8)
  )

F2.2

# ggsave("Figure 2/Figure 2.2.pdf", F2.2, width = 6, height = 12)

####F2.3 环境因素与代谢综合征的关联====
covars <- c("citytype","agegroup","gender","educagroup","nationg",
            "marrygroup","occupation_g","income1","fuelg",
            "smoke","drinkgroup3","total_ng","crop","vegetable","fruit","meat","season")

env_factors <- c("PM25","NO2","TEMP")

meta <- as.data.frame(meta)


##IQR
iqr_vals <- sapply(meta[env_factors], IQR, na.rm = TRUE)  # 计算 IQR

for (v in env_factors) {
  meta[[v]] <- meta[[v]] / iqr_vals[[v]]
}

res <- data.frame()
for (i in env_factors) {
  
  fml <- as.formula(paste0("mets ~ ", paste(c(i, covars), collapse = " + ")))
  model <- glm(fml, data = meta, family = binomial())
  
  # 提取回归结果
  beta <- coef(model)[[i]]
  se   <- summary(model)$coefficients[i, 2]
  pval <- summary(model)$coefficients[i, 4]
  
  OR     <- exp(beta)
  CI_low <- exp(beta - 1.96 * se)
  CI_high<- exp(beta + 1.96 * se)
  
  res <- rbind(res, data.frame(
    EnvFactor = i,
    OR = OR, CI_low = CI_low, CI_high = CI_high,
    logOR = beta, SE = se, PValue = pval
  ))
}

#绘图
F2.3 <- ggplot(res, aes(x = EnvFactor, y = OR)) +
  geom_point(size = 3, color = "#B40426") +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                width = 0.2, linewidth = 0.8, color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  ylim(0.7, 1.4)+
  labs(x = "Environmental Factor (per IQR increase)",
       y = "Odds Ratio (95% CI)",
       title = "") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(face = "bold", color = "black",angle = 45,hjust = 1, vjust = 1),
    axis.text.y = element_text(face = "bold", color = "black"),
    plot.title  = element_text(face = "bold"),
    plot.margin = margin(t = 5.5, r = 20, b = 5.5, l = 20, unit = "pt")
  )

F2.3

# ggsave("Figure 2/Figure2.3_OR.pdf", F2.3, width = 6, height = 6)


####F2.4中介分析=====
res_A <- fread("Figure 4/res_a.csv")
res_B <- fread("Figure 4/res_b.csv")

res_A <- res_A %>%  rename(otuid = Genus)
tax <- tax %>% rownames_to_column("otuid")
res_A <- res_A %>% left_join(tax,by = "otuid")

res_A <- res_A[res_A$Genus %in% genus_union, ]


###按照环境因子分组矫正FDR
res_A <- res_A %>%
  group_by(EnvFactor) %>%
  mutate(FDR_ACME = p.adjust(ACME_p, method = "BH"),
         FDR_ADE = p.adjust(ADE_p,   method = "BH"),
         FDR_Total = p.adjust(Total_p, method = "BH")) %>%
  ungroup() %>%
  filter(FDR_ACME < 0.1,
        FDR_Total < 0.1,
         ACME * TotalEffect > 0)

res_A <- res_A %>%
  filter(EnvFactor %in% c("NO2","PM25","TEMP"))

###plot
df_pairs <- res_A%>%
  mutate(
    EnvFactor = factor(EnvFactor,
                       levels = c("NO2","PM25","TEMP"),
                       labels = c("NO2","PM2.5","Temperature")),
    Genus = factor(Genus)
  ) %>%
  group_by(EnvFactor, Genus) %>%
  mutate(pair_id = cur_group_id()) %>%
  ungroup()


df_pairs <- df_pairs %>%
  mutate(fill_group = EnvFactor)               # 这列会在转长后保留

lodes_df <- to_lodes_form(
  data = df_pairs,
  key  = "axis",                               # 轴名列
  axes = c("EnvFactor", "Genus"),              # 左→右顺序
  id   = "pair_id"                             # alluvium ID
)                                              

## 3) 画图
env_cols <- c(
  "NO2"              = "#E64B35FF",  # NPG red
  "PM2.5"            = "#4DBBD5FF",  # NPG blue
  "Temperature"      = "#00A087FF"  # NPG green
)

genus_levels <- levels(df_pairs$Genus)

# 至少取 8 个颜色，再根据菌属数量截取前 n 个
genus_palette <- ggsci::pal_npg("nrc")(max(length(genus_levels), 8))

genus_cols <- setNames(
  alpha(genus_palette, 0.6)[seq_along(genus_levels)],  # 降低饱和度/透明度
  genus_levels
)

## 3. 合并颜色
fill_values <- c(env_cols, genus_cols)

# 合并所有颜色，并让图例只显示环境因子
fill_values <- c(env_cols, genus_cols)

F2.4 <- ggplot(
  lodes_df,
  aes(x = axis, stratum = stratum, alluvium = pair_id, y = PropMediated)
) +
  # 条带按环境因子着色
  geom_alluvium(aes(fill = fill_group), knot.pos = 0.4, alpha = 0.9) +
  
  # 左侧：环境因子柱——用环境因子颜色
  geom_stratum(
    data = subset(lodes_df, axis == "EnvFactor"),
    aes(fill = fill_group),
    width = 1/5, color = "grey40"
  ) +
  
  # 右侧：属柱——按属各自颜色（不进图例）
  geom_stratum(
    data = subset(lodes_df, axis == "Genus"),
    aes(fill = after_stat(stratum)),
    width = 1/4, color = "grey40",
    show.legend = FALSE
  ) +
  
  # 左右两侧标签
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 3
  ) +
  
  # 统一的填充尺度：包含“环境因子色 + 属色”，但图例只展示环境因子
  scale_fill_manual(
    values = fill_values,
    breaks = names(env_cols),   # 只显示环境因子图例项
    name = NULL
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid   = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom"
  )

F2.4

ggsave("Figure 2/Figure 2.4.pdf", F2.4, width = 12, height = 6)


