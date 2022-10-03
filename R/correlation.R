# correlation between normalised counts of IDO
library("ggrepel")
library("viridis")
# gene_id == "ENSG00000131203"

# get normalised counts
countsNormalized <- DESeq2::counts(dds, normalized = TRUE) %>% data.frame()
countsNormalized$gene_id <- row.names(countsNormalized)

temp <- countsNormalized %>%
  filter(gene_id == "ENSG00000131203") %>%
  select(-gene_id) %>%
  t() %>% as.data.frame()

temp$sample <- row.names(temp)
temp <- temp %>%
  extract(sample, remove = FALSE,
          c("samp", "donor", "treatment", "idoresp"),
          "(S[0-9]{2})_(P[0-9]{5})_([a-z]{6})_([a-z]{2})") %>%
  mutate(names = str_replace(sample,"TMP_", ""))
coldata$names <- as.character((coldata$names))

temp <- coldata %>% select(names, fcido) %>%
  merge(temp, by = "names")

ido_cor <- temp %>%
  filter(treatment == "stimul") %>%
  ggplot(aes(x = fcido, y = ENSG00000131203, colour = idoresp)) +
  geom_point(size = 3) +
  xlab("IDO (PCR)") +
  ylab("(DESeq2) Normalised counts of IDO") +
  scale_color_manual(name = "IDO response\n class",
                       values = viridis(end = 0.6, n = 2),
                     labels = c("High", "Low")) +
  geom_text_repel(aes(label = donor)) +
  ggtitle("Correlation between PCR'd IDO and IDI transcripts") +
  theme_classic()

ggsave(filename = "reports/figures/ido-correlation.png",
       plot = ido_cor,
       units = "in",
       height = 5,
       width = 5,
       dpi = 300,
       device = "png")

temp <- temp %>%
  filter(treatment == "stimul")
cor.test(temp$fcido, temp$ENSG00000131203)



# correlation between normalised counts of IDO
