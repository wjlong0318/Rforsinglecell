#install.packages("immunarch")
rm(list = ls()) 
library(immunarch)
file_path = "H:\\shanghai\\singlecell\\P-MATRIX\\P-164-B"
immdata_10x <- repLoad(file_path)
immdata_10x$data
immdata=immdata_10x
div_chao <- repDiversity(immdata$data$filtered_contig_annotations, "chao1")

# Hill numbers
div_hill <- repDiversity(immdata$data, "hill")

# D50
div_d50 <- repDiversity(immdata$data, "d50")

# Ecological diversity measure
div_div <- repDiversity(immdata$data, "div")



p1 <- vis(div_chao)
p2 <- vis(div_chao, .by = c("Status", "Sex"), .meta = immdata$meta)
p3 <- vis(div_hill, .by = c("Status", "Sex"), .meta = immdata$meta)

p4 <- vis(div_d50)
p5 <- vis(div_d50, .by = "Status", .meta = immdata$meta)
p6 <- vis(div_div)


repExplore(immdata$data, .method = "volume")

exp_len <- repExplore(immdata$data, .method = "len", .col = "aa")
exp_cnt <- repExplore(immdata$data, .method = "count")
exp_vol <- repExplore(immdata$data, .method = "volume")

p1 <- vis(exp_len)
p2 <- vis(exp_cnt)
p3 <- vis(exp_vol)

target <- c("CARDKSMVRGVEAFDIW", "CASSDSSGGANEQFF", "CASSDSSGSTDTQYF", "CASSLAGGYNEQFF", "CASSDSAGGTDTQYF", "CASSLDSYEQYF", "CASSSAGGYNEQFF")
tc <- trackClonotypes(immdata$data, target, .col = "aa")
vis(tc)
################################################################
#SD  P-177  P-179   P-184  P-185  P-186  P-187
#PR  P-164  P-165   P-176   P-181   P-189   P-191   P-198  P-199
#PD  P-174   P-175   P-178   P-180   P-192   P-193   P-197   P-202
# sample= c("P-177","P-179","P-184","P-185","P-186","P-187",#SD  
#  "P-164","P-165","P-176","P-181","P-189","P-191","P-198","P-199",#PR 
# "P-174","P-175","P-178","P-180","P-192","P-193","P-197","P-202")#PD
rm(list = ls()) 
library(immunarch)
setwd("H:\\shanghai\\singlecell\\P-MATRIX\\")
meta0=read.csv("../meta.csv",header=T,sep=",")
meta=meta0[meta0$sample!="P-150"&meta0$Response!="Na"&meta0$Treatment!=5,]
meta$sample
meta$RS=paste(meta$Response,meta$sample,sep="_")
imm=NULL
for(s in c(1:dim(meta)[1])){
  #s=1
  file_path = paste("H:\\shanghai\\singlecell\\P-MATRIX\\",meta$sample[s],"-B",sep="")
  immdata_10x <- repLoad(file_path)
  imm$data[[meta$RS[s]]]=immdata_10x$data$filtered_contig_annotations
}
saveRDS(imm, file="Tall_immunarch.rds")
#pbmc<- readRDS("sample.combined_T10_others_singleR_ref3_main_fine_celltype.rds")

div_chao <- repDiversity(imm$data,"chao1")
vis(div_chao)#.by = c("Response","Gender"),.meta=meta
div_hill <- repDiversity(imm$data, "hill")
vis(div_hill)
div_d50 <- repDiversity(imm$data, "d50")
vis(div_d50,.by = meta$Response)
div_div <- repDiversity(imm$data, "div")
vis(div_div)#,.by = meta$Response)


exp_vol <- repExplore(imm$data, .method = "volume")
vis(exp_vol)#,.by = meta$Response)
exp_len <- repExplore(imm$data, .method = "len", .col = "aa")

vis(exp_len,.by = "Response",.meta=meta)
fixVis(p1)
names(imm) 
exp_cnt <- repExplore(imm$data, .method = "count")
exp_vol <- repExplore(imm$data, .method = "volume")
vis(exp_cnt,.by = "Response",.meta=meta)
vis(exp_vol)#,.by = meta$Response)
imm_pr <- repClonality(imm$data, .method = "clonal.prop")
imm_top <- repClonality(imm$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
imm_top

imm_rare <- repClonality(imm$data, .method = "rare")
imm_rare
vis(imm_rare)
imm_hom <- repClonality(imm$data,
                        .method = "homeo",
                        .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1)
)
imm_hom
vis(imm_hom)

vis(imm_top)+facet_wrap(~.Response)
vis(imm_top, .by ="Response", .meta=new_meta)

new_meta$sample=new_meta$RS
new_meta=as_tibble(meta)

imm_top1 <- repClonality(immdata$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
vis(imm_top1, .by = "Status", .meta = immdata$meta)


imm_ov1 <- repOverlap(imm$data, .method = "public", .verbose = F)
imm_ov2 <- repOverlap(imm$data, .method = "morisita", .verbose = F)

p1 <- vis(imm_ov1)
p2 <- vis(imm_ov2, .text.size = 2)

vis(imm_ov1, "heatmap2")
p1 <- vis(imm_ov2, .text.size = 2.5, .signif.digits = 1)
p2 <- vis(imm_ov2, .text.size = 2, .signif.digits = 2)
repOverlapAnalysis(imm_ov1, "mds")
repOverlapAnalysis(imm_ov1, "tsne")
repOverlapAnalysis(imm_ov1, "tsne") %>% vis()
repOverlapAnalysis(imm_ov1, "tsne+kmeans") %>% vis()

imm_gu <- geneUsage(imm$data[c(1, 2)], "hs.trbv", .norm = T)

vis(imm_gu)

imm_gu <- geneUsage(imm$data, "hs.trbv", .norm = T)

imm_gu_js <- geneUsageAnalysis(imm_gu, .method = "js", .verbose = F)
imm_gu_cor <- geneUsageAnalysis(imm_gu, .method = "cor", .verbose = F)

p1 <- vis(imm_gu_js, .title = "Gene usage JS-divergence", .leg.title = "JS", .text.size = 1.5)
p2 <- vis(imm_gu_cor, .title = "Gene usage correlation", .leg.title = "Cor", .text.size = 1.5)

imm_gu_js[is.na(imm_gu_js)] <- 0

vis(geneUsageAnalysis(imm_gu, "cosine+hclust", .verbose = F))

imm_cl_pca <- geneUsageAnalysis(imm_gu, "js+pca+kmeans", .verbose = F)
imm_cl_mds <- geneUsageAnalysis(imm_gu, "js+mds+kmeans", .verbose = F)
imm_cl_tsne <- geneUsageAnalysis(imm_gu, "js+tsne+kmeans", .perp = .01, .verbose = F)
## Perplexity should be lower than K!
p1 <- vis(imm_cl_pca, .plot = "clust")
p2 <- vis(imm_cl_mds, .plot = "clust")
p3 <- vis(imm_cl_tsne, .plot = "clust")

p1 <- vis(spectratype(imm$data[[1]], .quant = "id", .col = "nt"))
p2 <- vis(spectratype(imm$data[[1]], .quant = "count", .col = "aa+v"))

imm_raref <- repDiversity(imm$data, "raref", .verbose = F)

p1 <- vis(imm_raref)
p2 <- vis(imm_raref, .by = "Response", .meta = new_meta)

vdjdb = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/vdjdb.slim.txt.gz", "vdjdb", .species = "HomoSapiens", .chain = "TRB", .pathology = "CMV")
vdjdb
dbAnnotate(imm$data, vdjdb, "CDR3.aa", "cdr3")
vdjdb_st = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/SearchTable-2019-10-17%2012_36_11.989.tsv.gz", "vdjdb-search", .species = "HomoSapiens", .chain = "TRB", .pathology = "CMV")
vdjdb_st
mcpas = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/McPAS-TCR.csv.gz", "mcpas", .species = "Human", .chain = "TRB", .pathology = "Cytomegalovirus (CMV)")
mcpas$CDR3.beta.aa
a=dbAnnotate(imm$data, mcpas, c("CDR3.aa", "V.name"), c("CDR3.beta.aa", "TRBV"))
names(imm$data$`PR_P-164`$CDR3.aa)
summary(a)
tbadb = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/TBAdb.xlsx", "tbadb", .species = "Homo Sapiens", .chain = c("TRB", "TRA-TRB"), .pathology = "CMV")
tbadb

kmers <- getKmers(imm$data[[1]], 5)
vis(kmers)

p1 <- vis(kmers, .head = 5)
p2 <- vis(kmers, .head = 10)
p3 <- vis(kmers, .head = 30)
p3 <- vis(kmers, .head = 5, .position = "dodge")

dbAnnotate(imm$data, vdjdb, "CDR3.aa", "cdr3")

kp <- kmer_profile(kmers, "self")
p1 <- vis(kp)
p2 <- vis(kp, .plot = "seq")
