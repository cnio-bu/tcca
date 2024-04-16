library(vegan)
library(tidyverse)

ith_by_sample <- data.table::fread(
    input  = "results/modules/annotated/metacom_proportions_primary_wide.tsv"
)

source("src/figures/TCCA_palette.R")

## start with naive primaries
ith_primaries_naive <- ith_by_sample %>%
    filter(
        study != "cell_lines_gabriella_kinker",
        sample_type == "p",
        treated == FALSE
    ) %>%
    filter(
        !(study == "brca_bhupinder_pal" & tumor_subtype == "predicted_tumour")
    ) %>%
    select(
        sample_study,
        tumor_type,
        sample_type,
        treated,
        n.prop_metacom_untreated_1,
        n.prop_metacom_untreated_2,
        n.prop_metacom_untreated_3,
        n.prop_metacom_untreated_4,
        n.prop_metacom_untreated_5,
        n.prop_metacom_untreated_6,
        best_metacom
    ) 

## comm data
ith_primaries_naive_mat <- ith_primaries_naive %>%
    select(sample_study, n.prop_metacom_untreated_1:n.prop_metacom_untreated_6) %>%
    as.data.frame()

rownames(ith_primaries_naive_mat) <- ith_primaries_naive_mat$sample_study
ith_primaries_naive_mat$sample_study <- NULL

ith_primaries_naive_mat <- as.matrix(ith_primaries_naive_mat)
ith_primaries_naive_mat[is.na(ith_primaries_naive_mat)] <- 0

shan <- vegan::diversity(
    x = ith_primaries_naive_mat,
    index = "shannon",
    MARGIN = 1,
    #   groups = ith_primaries_naive$tumor_type,
    #  equalize.groups = TRUE
)
ith_primaries_naive$shan <- shan
ith_primaries_naive$tumor_type <- fct_reorder(ith_primaries_naive$tumor_type, ith_primaries_naive$shan)

## genomic ith
genomic_ith <- read_tsv("results/cna/genomic_ith.tsv")

## add n.cells by sample
tcs <- read_tsv("results/annotation/beyondcell_with_therapeutic_clusters.tsv") %>%
    mutate(
        sample_study = paste0(study, "__", sample)
    ) %>%
    group_by(sample_study) %>%
    mutate(
        n.cells = n()
    ) %>%
    filter(
        !is.na(cell)
    )

genomic_ith <- genomic_ith %>%
    left_join(
        y = tcs[, c("sample_study", "n.cells")],
        by = c("study__sample" = "sample_study")
    ) %>%
    distinct() 

genomic_ith_rates <- genomic_ith %>%
    mutate(
        n.prop = nclones * 1000 / n.cells,
        sample_study = paste0(sample, study, sep = "_")
    )

ith_primaries_naive$shan <- shan

genomic_therapeutic_ith <- ith_primaries_naive %>%
    left_join(
        y = genomic_ith_rates,
        by = "sample_study"
    ) %>%
    filter(
        !is.na(n.prop),
        !is.na(shan)
    ) %>%
    group_by(tumor_type) %>%
    mutate(
        mean_shan = median(shan),
        avg.clones = mean(n.prop),
        n.samples = n()
    ) %>%
    filter(
        n.samples > 5
    )



mean_shan_levels <- ggplot(genomic_therapeutic_ith,
                           aes(x = tumor_type,
                               y = shan,
                               color = tumor_type
                           )
) +
    geom_segment(aes(x =  tumor_type,
                     xend = tumor_type,
                     y = median(genomic_therapeutic_ith$shan),
                     yend = mean_shan
    )
    ) +
    geom_jitter(size = 3, alpha = 0.25, width = 0.2) +
    geom_hline(aes(yintercept = median(genomic_therapeutic_ith$shan)), 
               color = "gray70",
               size = 0.6
    ) +
    stat_summary(fun = median, geom = "point", size = 5) +
    stat_summary(aes(label = sample),
                 geom = "text",
                 fun.y = function(y){ o <- jitter.stats(y)$out; if (o >= 3000) o else NA}
    ) +
    annotate("text",
             x = 5.5,
             y = 2000,
             size = 3.8,
             color = "gray20",
             lineheight = .9,
             label = "Pan-cancer therapeutic heterogeneity"
    ) + 
    coord_flip() +
    scale_y_continuous(limits = c(1, 2), expand = c(0,0)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "none",
          text = element_text(family = "Arial"),
          plot.caption = element_text(size = 9, color = "gray50"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()
    ) +
    labs(y = "Therapeutic heterogeneity index",
         x = ""
    )

ggsave(
    plot = mean_shan_levels,
    filename = "results/figures/mean_shan_levels.png",
    dpi = 100,
    height = 7,
    width = 7
)


mean_clonal_levels <- ggplot(genomic_therapeutic_ith,
                             aes(x = tumor_type,
                                 y = n.prop,
                                 color = tumor_type
                             )
) +
    geom_segment(aes(x =  tumor_type,
                     xend = tumor_type,
                     y = mean(genomic_therapeutic_ith$n.prop),
                     yend = avg.clones
    )
    ) +
    geom_jitter(size = 3, alpha = 0.25, width = 0.2) +
    geom_hline(aes(yintercept = mean(genomic_therapeutic_ith$n.prop)), 
               color = "gray70",
               size = 0.6
    ) +
    stat_summary(fun = mean, geom = "point", size = 5) +
    stat_summary(aes(label = sample),
                 geom = "text",
                 fun.y = function(y){ o <- jitter.stats(y)$out; if (o >= 3000) o else NA}
    ) +
    annotate("text",
             x = 5.5,
             y = 2000,
             size = 3.8,
             color = "gray20",
             lineheight = .9,
             label = "Pan-cancer clonal heterogeneity"
    ) + 
    coord_flip() +
    scale_y_continuous(limits = c(0, 60), expand = c(0,0)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "none",
          text = element_text(family = "Arial"),
          plot.caption = element_text(size = 9, color = "gray50"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()
    ) +
    labs(y = "Clone proportion by 1000 cells",
         x = ""
    )

ggsave(
    plot = mean_clonal_levels,
    filename = "results/figures/mean_clonal_levels.png",
    dpi = 100,
    height = 7,
    width = 7
)



correlation_wrap <- ggplot(data = genomic_therapeutic_ith, aes(x = shan, y = n.prop)) +
    geom_point() +
    geom_smooth(method = "lm") +
    ggpubr::stat_cor(method = "spearman") +
    facet_wrap(~tumor_type, ncol = 5, nrow = 5) +
    theme_bw() + 
    theme(panel.grid = element_blank(),
          legend.position = "none",
          text = element_text(family = "Arial"),
          plot.caption = element_text(size = 9, color = "gray50"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()
    )


ggsave(
    plot = correlation_wrap,
    filename = "results/figures/clonal_shan_correlation_cancer_type.png",
    dpi = 300,
    height = 21,
    width = 21
)

avg_shan_by_metacom <- ith_primaries_naive %>%
    group_by(best_metacom) %>%
    summarise(
        avg_shan = median(shan)
    )
## Average shannon and clonal heterogeneity by therapeutic module
metacom_by_shannon <- ggplot(
    data = ith_primaries_naive,
    aes(x = best_metacom, y = shan)
) +
    geom_boxplot() +
    scale_x_discrete(labels= c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6")) +
    xlab("") +
    ylab("Shannon diversity") +
    ggpubr::stat_pwc(method = "wilcox.test", p.adjust.method = "BH", hide.ns = TRUE) +
    theme_minimal() +
    theme(
        panel.grid = element_blank(),
        legend.position = "none",
        text = element_text(family = "Arial"),
        plot.caption = element_text(size = 9, color = "gray50"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
    )


ggsave(
    filename = "results/figures/metacom_by_shannon.png",
    plot = metacom_by_shannon,
    dpi = 100,
    height = 7,
    width = 7
    )

## Average shannon and clonal heterogeneity by therapeutic module
avg_clonal_by_metacom <- genomic_therapeutic_ith %>%
    group_by(best_metacom) %>%
    summarise(
        avg_metacom_clonal = median(avg.clones)
    )

metacom_by_clonal <- ggplot(
    data = genomic_therapeutic_ith,
    aes(x = best_metacom, y = avg.clones)
) +
    geom_point() +
    geom_boxplot() +
    scale_x_discrete(labels= c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6")) +
    xlab("") +
    ylab("Clonal diversity") +
    ggpubr::stat_pwc(method = "wilcox.test", p.adjust.method = "BH", hide.ns = TRUE) +
    theme_minimal() +
    theme(
        panel.grid = element_blank(),
        legend.position = "none",
        text = element_text(family = "Arial"),
        plot.caption = element_text(size = 9, color = "gray50"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
    )


ggsave(
    filename = "results/figures/metacom_by_clones.png",
    plot = metacom_by_clonal,
    dpi = 100,
    height = 7,
    width = 7
)

## test single plot
genomic_clona_ith_long <- genomic_therapeutic_ith %>%
    pivot_longer(
        cols = c("mean_shan", "avg.clones"),
        names_to = "diversity",
        values_to = "val"
        )

metacom_clone_ith_single <- ggplot(
    data = genomic_clona_ith_long,
    aes( x= best_metacom,
         y = val,
         fill = diversity
         )
    ) + geom_boxplot()



## Final figures

## Exclude LAML since scEVAN performed notably bad
genomic_therapeutic_ith_valid <- genomic_therapeutic_ith %>%
    filter(tumor_type != "LAML")


module_clonal_diversity <- ggplot(
    data = genomic_therapeutic_ith_valid,
    aes(x = best_metacom, y = n.prop, fill = best_metacom)
) +
    geom_boxplot() +
    geom_jitter() +
    scale_x_discrete(
        labels= c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6")
        ) +
    scale_y_continuous(n.breaks = 10) + 
    scale_fill_manual(values = unname(module_colors)) +
    xlab("") +
    ylab("Clonal diversity") +
    ggpubr::stat_pwc(
        method = "wilcox.test",
        p.adjust.method = "BH",
        hide.ns = TRUE
        ) +
    theme_minimal() +
    theme(
        panel.grid = element_blank(),
        legend.position = "none",
        text = element_text(family = "Arial"),
        plot.caption = element_text(size = 9, color = "gray50"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
    )


ggsave(
    filename = "results/figures/module_clonal_diversity.png",
    plot = module_clonal_diversity,
    dpi = 100,
    height = 7,
    width = 7
       )

## Average shannon and clonal heterogeneity by therapeutic module
module_ith_diversity <- ggplot(
    data = genomic_therapeutic_ith_valid,
    aes(x = best_metacom, y = shan, fill = best_metacom)
) +
    geom_boxplot() +
    geom_jitter() +
    scale_x_discrete(labels= c("TM 1", "TM 2", "TM 3", "TM 4", "TM 5", "TM 6")) +
    scale_y_continuous(n.breaks = 10) + 
    scale_fill_manual(values = unname(module_colors)) +
    xlab("") +
    ylab("Shannon diversity") +
    ggpubr::stat_pwc(method = "wilcox.test", p.adjust.method = "BH", hide.ns = TRUE) +
    theme_minimal() +
    theme(
        panel.grid = element_blank(),
        legend.position = "none",
        text = element_text(family = "Arial"),
        plot.caption = element_text(size = 9, color = "gray50"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
    )

ggsave(
    filename = "results/figures/module_ith_diversity.png",
    plot = module_ith_diversity,
    dpi = 100,
    height = 7,
    width = 7
)
