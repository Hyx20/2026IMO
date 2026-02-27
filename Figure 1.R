library(tidyverse)
library(ggplot2)
library(data.table)
library(ggsci)
library(tidytext)
library(showtext)
library(scales)
library(ggnewscale)
library(vegan)
library(ape)
library(umap)
library(slingshot)
library(gghalves)
library(ggpubr)
library(ggforce)  # geom_arc_bar
library(ggh4x)
showtext_auto() 

####Figure 2.1 Alpha多样性====
alpha <- fread("data/F2.1_alpha.csv")
alpha$mets_group <- factor(alpha$mets,levels = c("Non-MetS", "MetS"), labels = c("Non-MetS", "MetS"))

group_colors <- c("Non-MetS" = "#1f77b4", "MetS" = "#d62728")


alpha$index <- factor(alpha$index,
                      levels = unique(alpha$index))

idx_lvls   <- levels(alpha$index)  # 
strip_cols <- c("#4DBBD5","#E64B35","#EFC000","#00A087")[seq_along(idx_lvls)]


F2.a <- ggplot(alpha, aes(x = mets_group, y = value)) +
  ## 半小提琴：左=Non-MetS，右=MetS
  geom_half_violin(
    data = subset(alpha, mets_group == "Non-MetS"),
    aes(fill = mets_group),
    side = "l", width = 0.8, trim = FALSE,
    alpha = 0.7, color = NA
  ) +
  geom_half_violin(
    data = subset(alpha, mets_group == "MetS"),
    aes(fill = mets_group),
    side = "r", width = 0.8, trim = FALSE,
    alpha = 0.7, color = NA
  ) +
  geom_point(
    data = subset(alpha, mets_group == "Non-MetS"),
    aes(x = as.numeric(mets_group) + 0.15, color = mets_group),
    position = position_jitter(width = 0.07),
    alpha = 0.7
  ) +
  geom_point(
    data = subset(alpha, mets_group == "MetS"),
    aes(x = as.numeric(mets_group) - 0.15, color = mets_group),
    position = position_jitter(width = 0.07),
    alpha = 0.7
  ) +
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    comparisons = list(c("Non-MetS", "MetS")),
    tip.length = 0.01,
    label.y.npc = 0.98,
    size = 2.8
    # 如需强制显著性文字为黑色，可加：color = "black"
  ) +
  facet_wrap2(
    ~ index, scales = "free_y", nrow = 1,
    strip = strip_themed(
      background_x = elem_list_rect(fill = strip_cols, color = NA),
      text_x = elem_list_text(
        color = "black",   # ← 改这里
        face = "plain",
        size = 8,
        family = "Arial"
      )
    )
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(x = "", y = "") +
  coord_cartesian(clip = "off") +
  theme_bw(base_family = "Arial", base_size = 8) +
  theme(
    text = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    strip.text = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    axis.text  = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    axis.title = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    legend.position = "none",
    strip.background = element_blank()
  )

F2.a
# ggsave("Figure 2/Figure 2a.pdf", F2.a,width = 12, height = 6, dpi = 300)

####Figure 2.2 Beta多样性====
otu <- fread("data/otu.csv")

# 处理OTU表，提取分类信息
tax <- otu %>%
  select(Kingdom, Phylum, Class, Order, Family, Genus) %>%
  distinct()

tax$otu_id <- paste0("OTU",1:nrow(tax))

#计算相对丰度
otu_rel <- otu %>%
  select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus) %>%
  as.data.frame()

otu_rel <- sweep(otu_rel, 2, colSums(otu_rel), FUN = "/")

otu_rel <- t(otu_rel) %>% as.data.frame()
colnames(otu_rel) <- tax$otu_id

# 计算 Bray-Curtis 距离矩阵
bc_dist <- vegdist(otu_rel, method = "bray")

###PCoA
pcoa_result <-pcoa(bc_dist)

pcoa_mat <- as.data.frame(pcoa_result$vectors[, 1:50])
rownames(pcoa_mat) <- rownames(pcoa_result$vectors)

set.seed(123)
umap_res <- umap(pcoa_mat,
                 n_neighbors = 15,
                 min_dist = 0.1,
                 metric = "euclidean")

umap_df <- as.data.frame(umap_res$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")

wss <- sapply(1:10, function(k) kmeans(umap_df, centers = k, nstart = 10)$tot.withinss)
plot(1:10, wss, type = "b", pch = 16, xlab = "Cluster Number", ylab = "Within-cluster SS")


# 使用 kmeans 聚类
km_res <- kmeans(umap_df, centers = 4, nstart = 20)
umap_df$Cluster <- factor(km_res$cluster)

umap_df$SampleID <- rownames(umap_df)
umap_df <- merge(umap_df, meta, by = "SampleID")

table(umap_df$Cluster)
####改变聚类的名字顺序
mapping <- c("1" = "2",
             "2" = "3",
             "3" = "4",
             "4" = "1")

umap_df$Cluster <- mapping[umap_df$Cluster]

table(umap_df$Cluster)

umap_df$Cluster <- as.factor(umap_df$Cluster)

center_coords <- umap_df %>%
  group_by(Cluster) %>%
  summarise(x = mean(UMAP1), y = mean(UMAP2))

pie_data <- umap_df %>%
  group_by(Cluster, mets) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  arrange(Cluster, mets) %>%
  group_by(Cluster) %>%
  mutate(
    start = cumsum(lag(freq, default = 0)) * 2 * pi,
    end = cumsum(freq) * 2 * pi
  ) %>%
  left_join(center_coords, by = "Cluster")

F2.b <- ggplot() +
  # 背景点
  geom_point(
    data = umap_df,
    aes(x = UMAP1, y = UMAP2, color = Cluster),
    size = 1
  ) +
  # 添加每个类中心的饼图
  geom_arc_bar(
    data = pie_data,
    aes(
      x0 = x, y0 = y, r0 = 0, r = 0.6,
      start = start, end = end, fill = factor(mets)
    ),
    color = "white", size = 0.2, alpha = 0.9
  ) +
  scale_fill_manual(
    values = c("0" = "#4575B4", "1" = "#D73027"),
    name = "MetS"
  ) +
  scale_color_npg() +
  coord_fixed() +
  theme_minimal(base_family = "Arial", base_size = 8) +
  theme(
    text = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    axis.text = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    axis.title = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    axis.line  = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    legend.title = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    plot.title = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    plot.subtitle = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    plot.caption = element_text(family = "Arial", size = 8, face = "plain", color = "black")
  )

F2.b

#拟时序分析
umap_df$Cluster <- as.factor(umap_df$Cluster)
sds <- slingshot::slingshot(umap_df[, c("UMAP1", "UMAP2")],
                            clusterLabels = umap_df$Cluster,
                            start.clus = "1") 

lineage_df1 <- as.data.frame(slingCurves(sds)$Lineage1$s)
colnames(lineage_df1) <- c("UMAP1", "UMAP2")

arrow_df1 <- tail(lineage_df1, 2)

# 合并箭头数据
arrow_df <- rbind(
  data.frame(x = arrow_df1$UMAP1[1], y = arrow_df1$UMAP2[1],
             xend = arrow_df1$UMAP1[2], yend = arrow_df1$UMAP2[2])
)


F2.b <- F2.b+geom_path(data = lineage_df1, aes(x = UMAP1, y = UMAP2),
                       color = "black", size = 1.2)+
  geom_segment(data = arrow_df,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.25, "inches")),
               color = "black", size = 1.2)
F2.b

# write(umap,"data/umap.csv")
# ggsave("Figure 2/Figure 2b.pdf", F2.b,width = 8, height = 6, dpi = 300)

####Figure 2.3 cluster FC====
meta <- fread("data/umap.csv")
clusters <- 1:4

stopifnot(all(rownames(otu_rel) %in% meta$SampleID))
otu_rel <- otu_rel[match(meta$SampleID, rownames(otu_rel)), , drop = FALSE]

res_list <- lapply(clusters, function(cl){
  in_idx  <- meta$Cluster == cl
  out_idx <- !in_idx
  
  mean_in  <- colMeans(otu_rel[in_idx, , drop = FALSE],  na.rm = TRUE)
  mean_out <- colMeans(otu_rel[out_idx, , drop = FALSE], na.rm = TRUE)
  
  # 对每个 OTU 做 Wilcoxon 检验
  pvals <- vapply(colnames(otu_rel), function(g){
    x <- otu_rel[in_idx,  g, drop = TRUE]
    y <- otu_rel[out_idx, g, drop = TRUE]
    
    wilcox.test(x, y, exact = FALSE)$p.value}, numeric(1))
  
  tibble(
    Cluster  = as.character(cl),
    otu_id   = names(mean_in),
    mean_in  = as.numeric(mean_in),
    mean_out = as.numeric(mean_out),
    log2FC   = log2((mean_in + 1e-6) / (mean_out + 1e-6)),
    p_value  = pvals
  ) %>%
    arrange(desc(log2FC))
})

fc_table <- bind_rows(res_list)

fc_table <- fc_table %>%
  mutate(q_value = p.adjust(p_value, method = "fdr"))

leaders_per_cluster <- fc_table %>%
  group_by(Cluster) %>%
  filter(log2FC > 1.5) %>% 
  arrange(desc(log2FC)) %>%
  ungroup()

leaders_per_cluster <- leaders_per_cluster %>% 
  left_join(tax, by = "otu_id") %>%
  select(Cluster, otu_id, Genus, log2FC, p_value, q_value)

###去掉未注释
leaders_per_cluster <- leaders_per_cluster %>%
  filter(Genus != "g__")

leaders_per_cluster <- leaders_per_cluster %>%
  group_by(Cluster) %>%
  mutate(Genus = fct_reorder(Genus, log2FC, .desc = TRUE)) %>%
  ungroup()

F2.c <- ggplot(leaders_per_cluster, aes(x = Genus, y = Cluster)) +
  geom_point(aes(size = log2FC), color = "steelblue", alpha = 0.8) +
  scale_size_continuous(range = c(3, 10)) +
  labs(
    x = "Genus",
    y = "Cluster",
    size = "log2FC",
    title = "All genera FDR < 0.001"
  ) +
  theme_bw(base_family = "Arial", base_size = 8) +
  theme(
    text = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    axis.text.x = element_text(
      family = "Arial", size = 8, face = "plain", color = "black",
      angle = 45, hjust = 1
    ),
    axis.text.y = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    axis.title  = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    plot.title  = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    legend.text = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    legend.title = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    axis.line  = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

F2.c
# ggsave("Figure 2/Figure 2c.pdf", F2.c,width = 16, height = 8, dpi = 300)


####Figure 2.4 不同cluster临床表型====
library(ggh4x)
library(rstatix)

dt <- fread("data/umap.csv")
dt$Cluster <- as.factor(dt$Cluster)
clinical_vars <- c("glu","waist","tg","hdl","sbp","dbp")

dt_long <- dt %>%
  select(Cluster, all_of(clinical_vars)) %>%
  pivot_longer(all_of(clinical_vars), names_to = "feature", values_to = "value") %>%
  filter(!is.na(value))

desired_order <- c("Waist","GLU","TG","HDL","SBP","DBP")
name_map <- c(GLU="GLU", Waist="Waist", TG="TG", HDL="HDL", SBP="SBP", DBP="DBP")

dt_long <- dt_long %>%
  mutate(feature = recode(feature,
                          bmi="BMI", glu="GLU", waist="Waist", tg="TG",
                          hdl="HDL", sbp="SBP", dbp="DBP",
                          .default = feature))
dt_long$feature <- factor(dt_long$feature, levels = desired_order)

feat_lvls  <- levels(dt_long$feature)
strip_cols <- c("#4DBBD5", "#E64B35", "#EFC000", "#00A087",
                "#3C5488", "#7E6148", "#91D1C2")[seq_along(feat_lvls)]

#每个分面单独的 y 轴上限
caps <- dt_long %>%
  group_by(feature) %>%
  summarise(
    q3  = quantile(value, 0.75, na.rm = TRUE),
    iqr = IQR(value, na.rm = TRUE),
    p98 = quantile(value, 0.98, na.rm = TRUE),
    cap = pmin(q3 + 2*iqr, p98) * 1.25,  # 上限= min(Q3+2IQR, P98) 再留5%余量
    .groups = "drop"
  )

#计算两两比较 + y 位置
stat_sig <- dt_long %>%
  group_by(feature) %>%
  pairwise_t_test(value ~ Cluster, p.adjust.method = "BH") %>%
  filter(p.adj < 0.05) %>%
  ungroup() %>%
  mutate(
    group1 = as.character(group1),
    group2 = as.character(group2)
  ) %>%
  group_by(feature) %>%
  arrange(group1, group2, .by_group = TRUE) %>%
  mutate(
    rank_in_facet = row_number(),
    n_layer       = n()
  ) %>%
  left_join(caps %>% select(feature, q3, iqr, cap), by = "feature") %>%
  mutate(
    top_margin = 0.03 * (cap - q3),
    y_top      = cap - top_margin,
    
    y_min      = q3 + 1.2 * iqr, 
    
    span       = pmax(0.05 * iqr, y_top - y_min),
    base_step  = span / (n_layer + 1),
    y_step     = pmax(0.01 * iqr, 0.6 * base_step),
    y_raw      = y_min + rank_in_facet * y_step,
    shift      = 0.10 * (cap - q3),
    y.position = pmin(y_top, y_raw + shift)
  ) %>%
  ungroup() %>%
  select(-q3, -iqr, -cap, -top_margin, -y_top,
         -y_min, -span, -base_step, -y_step, -y_raw, -shift)

y_scales <- lapply(seq_len(nrow(caps)), function(i){
  f  <- as.character(caps$feature[i])
  up <- caps$cap[i]
  eval(parse(text = sprintf('feature == "%s" ~ scale_y_continuous(limits = c(NA, %f))', f, up)))
})

#plot
F2.d <- ggplot(dt_long, aes(x = Cluster, y = value, fill = Cluster)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, size = 0.35) +
  stat_pvalue_manual(
    stat_sig,
    label = "p.adj.signif",
    bracket.size = 0.2,
    size = 2.8,
    color = "black",   # 显著性文字/括号颜色更稳
    tip.length = 0
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  coord_cartesian(clip = "off") +
  facet_wrap2(
    ~ feature, ncol = 3, scales = "free_y",
    labeller = as_labeller(name_map),
    strip = strip_themed(
      background_x = elem_list_rect(fill = strip_cols, color = NA),
      text_x = elem_list_text(
        color = "black",   # ← 改成黑色
        face = "plain",
        size = 8,
        family = "Arial"
      )
    )
  ) +
  facetted_pos_scales(y = y_scales) +
  scale_fill_npg() +
  labs(x = NULL, y = NULL, fill = "Cluster") +
  theme_bw(base_family = "Arial", base_size = 8) +
  theme(
    text = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    axis.text = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    axis.title = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    strip.text = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    legend.text = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    legend.title = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    plot.title = element_text(family = "Arial", size = 8, face = "plain", color = "black"),
    
    axis.line  = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.spacing.x = unit(2, "mm"),
    panel.spacing.y = unit(2, "mm"),
    legend.position = "top"
  )

F2.d
# ggsave("Figure 2/Figure 2d.pdf", F2.d, width = 8, height = 8, dpi = 300)

####Figure 2.5 网络图====
library(phyloseq)
library(microbiomeMarker)
library(Hmisc)      # rcorr
library(igraph)
library(ggraph)
library(scales)

meta <- fread("data/umap.csv")

meta$bmi_g <- cut(meta$BMI,
                  breaks = c(-Inf, 18.5, 24, 28, Inf),
                  labels = c("Underweight", "Normal", "Overweight", "Obese"))

table(meta$bmi_g)

otu <- fread("data/otu.csv")

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

###lefse 差异分析==
#总lesfe
otu_rel <- as.matrix(otu_rel)
tax <- as.matrix(tax)
meta <- meta %>% column_to_rownames("SampleID")

otu_tab <- otu_table(otu_rel, taxa_are_rows = F)
tax_tab <- tax_table(tax)
sample_tab <- sample_data(meta)

ps <- phyloseq(otu_tab, tax_tab, sample_tab)

###lefse分析前进行cpm标准化
ps_st <- transform_sample_counts(ps, function(x) x * 1e6)

marker_res <- run_lefse(
  ps_st,
  group = "mets",        # 分组变量
  lda_cutoff = 2,        # LDA阈值
  norm = "none",          # 归一化方法
)

m <- marker_table(marker_res) %>% as.matrix() %>% as.data.frame()

m <- m %>%
  separate(feature, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = "\\|", fill = "right")

m_genues <- m %>%
  filter(grepl("g__", Genus))

m_genues <- m_genues %>%
  filter(
    grepl("^g__", Genus),                         # 属级注释
    !grepl("g__$|g__.*_g[_]*$|g__.*f__.*", Genus)  # 去掉 g__、g__xxx_f__g_ 这类伪属
  )

#预处理：按 LDA 值排序 + 转换 p 值为 log10
m_genues <- m_genues %>%
  mutate(
    ef_lda = as.numeric(ef_lda),
    log10_FDR = -log10(as.numeric(pvalue))
  ) %>%
  group_by(enrich_group) %>%
  arrange(ef_lda) %>%
  mutate(Genus_ordered = factor(Genus, levels = unique(Genus))) %>%
  ungroup()

m_genues$enrich_group <- factor(m_genues$enrich_group, levels = c(0, 1), labels = c("Non-MetS", "MetS"))
# fwrite(m_genues, "Figure 2/lefse_genus_all.csv")

m_genues <- fread("Figure 2/lefse_genus_all.csv")
m_genues<- m_genues %>%
  left_join(tax %>% as.data.frame() %>%  rownames_to_column("otu_id") 
            %>% select(otu_id, Genus), by = "Genus")

###标准化
otu_rel <- log(otu_rel + 1e-6)

###Spearman correlation==
rc <- rcorr(as.matrix(otu_rel), type = "spearman")
R  <- rc$r
P  <- rc$P

edges <- as.data.frame(as.table(R)) %>%
  dplyr::rename(from = Var1, to = Var2, rho = Freq) %>%
  filter(as.character(from) < as.character(to))

pvals <- as.data.frame(as.table(P)) %>%
  dplyr::rename(from = Var1, to = Var2, p = Freq) %>%
  filter(as.character(from) < as.character(to))

edges <- left_join(edges, pvals, by = c("from","to")) %>%
  mutate(q = p.adjust(p, method = "BH")) %>%
  filter(!is.na(q), abs(rho) >= 0.3, q < 0.05) #过滤阈值 cor > 0.3, q < 0.05

##构建绘图数据
nodes <- data.frame(name = rownames(tax), stringsAsFactors = FALSE)

lefse_tags <- m_genues %>%
  transmute(
    name = otu_id,
    lefse_cat = case_when(
      enrich_group == "MetS"     ~ "MetS",
      enrich_group == "Non-MetS" ~ "Non-MetS",
      TRUE                       ~ NA_character_
    ),
    LDA = as.numeric(ef_lda)
  )

nodes <- nodes %>%
  left_join(lefse_tags, by = "name") %>%
  mutate(
    lefse_cat = ifelse(is.na(lefse_cat), "Non-significant", lefse_cat),
    lefse_cat = factor(lefse_cat, levels = c("MetS","Non-MetS","Non-significant")),
    LDA = replace_na(LDA, 0)
  )

nodes <- nodes %>%
  left_join(tax %>% as.data.frame() %>% rownames_to_column("name"), by = "name")

nodes$Genus <- gsub("^g__", "", nodes$Genus)

###igraph network ==
# 根据边表 (edges) 和点表 (nodes) 构建 igraph 对象 g，网络是无向的
g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)

# 计算每个节点的度数 (degree)，即连接边的数量，并存到节点属性 degree
V(g)$degree <- degree(g)

# 如果图里至少有一条边
if (ecount(g) > 0) {
  # 用 Louvain 算法进行社区发现 (聚类)，得到每个节点所属的社区
  com <- cluster_louvain(g)
  # 把社区编号作为节点属性 community
  V(g)$community <- factor(membership(com))
} else {
  # 如果没有边，就把 community 设置为 NA
  V(g)$community <- factor(NA)
}

# 把边的相关系数绝对值做 0.3–2.0 的线性缩放，用于绘图时控制边宽
E(g)$w_abs <- rescale(abs(E(g)$rho), to = c(0.3, 2.0))

##plot==
set.seed(123)
F2.e <- ggraph(g, layout = "fr") +
  # ---- edges: 只保留颜色图例，隐藏“宽度”图例 
  geom_edge_link(
    aes(color = rho, width = w_abs),
    alpha = 0.6,
    show.legend = c(edge_colour = TRUE, edge_width = FALSE)
  ) +
  scale_edge_color_gradient2(
    low = "#08519c",      # 深蓝，负相关
    mid = "white",        # 中性
    high = "#a50f15",     # 深红，正相关
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman"
  ) +
  scale_edge_width(range = c(0.5, 2.5), guide = "none") +
  
  # nodes: 大小=LDA, 填充=LEfSe 分类；边框固定黑色
  geom_node_point(
    aes(size = LDA, fill = lefse_cat),
    shape = 21, stroke = 0.6, color = "black",
    show.legend = TRUE
  ) +
  scale_size_continuous(
    name   = "LDA",
    limits = c(0, 4.5),
    breaks = 0:4,
    range  = c(2, 8),
    guide  = guide_legend(
      order = 3,
      override.aes = list(shape = 21, fill = "black", color = "black", linetype = 0)
    )
  ) +
  scale_fill_manual(
    values = c(
      "MetS" = "#a50f15",
      "Non-MetS" = "#08519c",
      "Non-significant" = "#bdbdbd"
    ),
    name = "Enriched group",
    guide = guide_legend(
      override.aes = list(shape = 21, size = 5, color = "black")
    )
  ) +
  # 把“size=LDA”的图例强制成圆点，去掉线段
  guides(
    colour = guide_colourbar(order = 1),   # Spearman (边颜色)
    fill   = guide_legend(order = 2, override.aes = list(
      shape = 21, size = 5, color = "black"
    )),
    size   = guide_legend(order = 3, override.aes = list(
      shape = 21, fill = "black", color = "black", linetype = 0
    ))
  ) +
  
  # 仅显示显著的标签
  geom_node_text(
    data = function(d) d %>% dplyr::filter(lefse_cat != "Non-significant"),
    aes(label = Genus),
    repel = TRUE,
    family = "Arial",
    size = 2.8,  # 约等于 8 pt 的视觉大小
    point.padding = grid::unit(0.3, "lines"),
    box.padding   = grid::unit(0.4, "lines"),
    force = 1,
    segment.size = 0.2,
    segment.alpha = 0.5
  ) +
  theme_void(base_family = "Arial", base_size = 8) +
  theme(
    text = element_text(family = "Arial", size = 8, face = "plain"),
    legend.position = "right",
    legend.text = element_text(family = "Arial", size = 8, face = "plain"),
    legend.title = element_text(family = "Arial", size = 8, face = "plain"),
    plot.title = element_text(
      family = "Arial", size = 8, face = "plain", hjust = 0.5
    )
  )

F2.e

# ggsave("Figure 2/Figure 2e.network.pdf", width = 10, height = 10, dpi = 300)

#####拼图=====
library(patchwork)

F2.up <- (F2.a | F2.b) / (F2.c) +
  plot_layout(widths = c(2, 1), heights = c(1.5, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "plain", size = 10, family = "Arial"))

F2.down <- (F2.d | F2.e) +
  plot_layout(widths = c(1, 1.3)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "plain", size = 10, family = "Arial"))

F2 <- F2.up / F2.down +
  plot_layout(heights = c(1.5, 1,3)) +   
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "plain", size = 10, family = "Arial"))

F2

ggsave("Figure 2/Figure 1.pdf", plot = F2, width = 12, height = 14, dpi = 300)
