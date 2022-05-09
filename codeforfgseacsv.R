ibrary(tidyverse)
res <- read.csv("postmort_expres.csv")
head(res)
row.names(res) <-res$X
res <- select(res, -X)
head(res)
colnames(res)
rownames(res)


#just want to care about stat and symbol 

res2 <- res %>% dplyr::select(symbol, stat) %>% na.omit() %>% distinct() %>% group_by(symbol) %>% summarize(stat = mean(stat))
res2
ranks <- deframe(res2)
pathways.hallmark <- gmtPathways("h.all.v7.5.1.symbols.gmt.txt")

pathways.hallmark %>% head() %>% lapply(head)
#run the fgsea with 1000 permutations 
fgseaRes <- fgsea(pathways=pathways.hallmark, stats = ranks)
head(fgseaRes)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


help <- fgseaResTidy %>% select(-leadingEdge, -ES) %>% arrange(padj) 

write.csv(help, file = "fgsea.csv")
