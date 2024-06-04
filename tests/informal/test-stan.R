# Load rosace functions
devtools::load_all(".")

# Load data
load("data/oct1_rosace.rda")
rosace <- oct1_rosace
key <- "1SM73"

# test SLR
rosace <- runSLR(rosace, name = key, type = "AssaySet")
rosace <- runSLR(rosace, name = paste(key, "1", sep = "_"), type = "Assay")
head(OutputScore(rosace, name = "1SM73_SLR"))
head(OutputScore(rosace, name = "1SM73_1_SLR"))
head(OutputScore.Score(rosace@scores[[1]]))
head(OutputScore.Score(rosace@scores[[2]]))

# small sample test
rosace@assays$`1SM73_1`@norm.counts <- rosace@assays$`1SM73_1`@norm.counts[1:100, ]
rosace@assays$`1SM73_1`@norm.var.names <- rosace@assays$`1SM73_1`@norm.var.names[1:100]
rosace@assay.sets$`1SM73`@raw.counts <- rosace@assay.sets$`1SM73`@raw.counts[1:100, ]
rosace@assay.sets$`1SM73`@combined.counts <- rosace@assay.sets$`1SM73`@combined.counts[1:100, ]
rosace@assay.sets$`1SM73`@var.names <- rosace@assay.sets$`1SM73`@var.names[1:100]

# test RunRosace Growth
# Assay (1 replicate)
# position model
rosace <- RunRosace(object = rosace,
                    name = "1SM73_1",
                    type = "Assay",
                    savedir = "tests/results/stan/assay",
                    pos.col = "position",
                    ctrl.col = "type",
                    ctrl.name = "synonymous",
                    stop.col = "type",
                    stop.name = "deletion",
                    install = FALSE)
head(OutputScore.Score(rosace@scores[[3]]))
head(OutputScore(rosace, name = "1SM73_1_ROSACE"))

OutputScore(rosace, name = "1SM73_1_ROSACE", pos.info = TRUE)

# non-position model
rosace <- RunRosace(object = rosace,
                    name = "1SM73_1",
                    type = "Assay",
                    savedir = "tests/results/stan/assay_nopos",
                    install = FALSE)
head(OutputScore.Score(rosace@scores[[3]]))
head(OutputScore(rosace, name = "1SM73_1_ROSACE"))

# AssaySet (multiple replicates)
# position model
rosace <- RunRosace(object = rosace,
                    name = "1SM73",
                    type = "AssaySet",
                    savedir = "tests/results/stan/assayset",
                    pos.col = "position",
                    ctrl.col = "type",
                    ctrl.name = "synonymous",
                    install = FALSE)
head(OutputScore.Score(rosace@scores[[4]]))
head(OutputScore(rosace, name = "1SM73_ROSACE"))

oct1_rosace_scored <- rosace
save(oct1_rosace_scored, file = "oct1_rosace_scored.rda", compress = "xz")

# non-position model
rosace <- RunRosace(object = rosace,
                    name = "1SM73",
                    type = "AssaySet",
                    savedir = "tests/results/stan/assayset_nopos",
                    install = FALSE)

# debug_list <- RunRosace(object = rosace,
#                         name = "1SM73",
#                         type = "AssaySet",
#                         savedir = "tests/results/stan/assayset",
#                         pos.col = "position",
#                         ctrl.col = "type",
#                         ctrl.name = "synonymous",
#                         debug = TRUE,
#                         install = FALSE)

# save(rosace, file = "tests/testdata/rosace_1SM73_eval.RData")


