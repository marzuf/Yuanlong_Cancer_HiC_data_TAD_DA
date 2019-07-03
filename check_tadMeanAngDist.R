dt_0 <- eval(parse(text = load("COSINE_BETWEEN_WITHIN_ALL/tad_angDist_fc_coexpr_DT.Rdata")))
dt_1 <- eval(parse(text = load("COSINE_BETWEEN_WITHIN_ALL_DIAG1//tad_angDist_fc_coexpr_DT.Rdata")))

colnames(dt_0)[colnames(dt_0) == "tad_meanAngDist"] <- "tad_meanAngDist_0"
colnames(dt_1)[colnames(dt_1) == "tad_meanAngDist"] <- "tad_meanAngDist_1"

colnames(dt_0)[colnames(dt_0) == "nbrGenes"] <- "nbrGenes_0"
colnames(dt_1)[colnames(dt_1) == "nbrGenes"] <- "nbrGenes_1"

cmp_dt <- merge(dt_0[,c("dataset", "region", "tad_meanAngDist_0", "nbrGenes_0")], 
                dt_1[,c("dataset", "region", "tad_meanAngDist_1", "nbrGenes_1")],
                by = c("dataset", "region"))

stopifnot(nrow(cmp_dt) == nrow(dt_0))
stopifnot(nrow(dt_1) == nrow(dt_0))

stopifnot(cmp_dt$nbrGenes_1 == cmp_dt$nbrGenes_0)
cmp_dt$nbrGenes <- cmp_dt$nbrGenes_0
cmp_dt$nbrGenes_0 <- NULL
cmp_dt$nbrGenes_1 <- NULL

plot(x = cmp_dt$tad_meanAngDist_1,
     y = cmp_dt$tad_meanAngDist_0
     )

oneThird_dt <- cmp_dt[cmp_dt$tad_meanAngDist_0 == (1/3),]
plot(x = oneThird_dt$tad_meanAngDist_1,
     y = oneThird_dt$tad_meanAngDist_0
)
plot(x = oneThird_dt$nbrGenes,
     y = oneThird_dt$tad_meanAngDist_0
)
