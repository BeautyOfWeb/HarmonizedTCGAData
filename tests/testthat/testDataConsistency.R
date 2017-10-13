library(HarmonizedTCGAData)
context("Data consistency")

library(ExperimentHub)
eh <- ExperimentHub()
myfiles <- query(eh, "HarmonizedTCGAData")
Wall <- myfiles[[1]]
project_ids <- myfiles[[2]]
surv.plot <- myfiles[[3]]


test_that("names of Wall", {
    cancer_types <- c("adrenal_gland", "lung",
                      "uterus", "kidney", "colorectal")
    feature_types <- c("raw.all", "raw.sel", "log.all",
                       "log.sel", "vst.sel", "normalized")
    data_types <- c("fpkm", "mirnas", "methy450")

    expect_equal(names(Wall), cancer_types)
    expect_equal(names(Wall[[1]]), feature_types)
    expect_equal(names(Wall[[1]][[1]]), data_types)
})


test_that("rownames of matrices in Wall", {
    rowname.correct <- all(sapply(Wall, function(x) sapply(x, function(y) {
        identical(rownames(y[[1]]), rownames(y[[2]])) &
            identical(rownames(y[[1]]), rownames(y[[3]])) &
            unique(nchar(rownames(y[[1]]))) == 12
    })))
    eval(bquote(expect_equal(rowname.correct, TRUE)))
})


test_that("colnames of matrices in Wall", {
    colname.correct <- all(sapply(Wall, function(x) sapply(x, function(y) {
        case.ids <- rownames(y[["fpkm"]])
        fpkm <- colnames(y[["fpkm"]])
        mirnas <- colnames(y[["mirnas"]])
        methy450 <- colnames(y[["methy450"]])

        nullnames.correct <- is.null(names(fpkm)) & !is.null(names(mirnas)) &
            !is.null(names(methy450))

        names(mirnas) <- NULL
        names(methy450) <- NULL
        match.col.correct <- identical(substr(fpkm, 1, 15), substr(mirnas, 1, 15)) &
            identical(substr(fpkm, 1, 15), substr(methy450, 1, 15)) &
            identical(substr(fpkm, 1, 12), case.ids)

        sample.type.correct <- unique(substr(fpkm, 14, 15)) == "01"

        length.correct <- unique(nchar(fpkm)) == 28 & unique(nchar(mirnas)) == 28 &
            unique(nchar(methy450)) == 28

        aliquot.type.correct <- unique(substr(fpkm, 20, 20)) == "R" &
            all(unique(substr(mirnas, 20, 20)) %in% c("R", "T", "H")) &
            unique(substr(methy450, 20, 20)) == "D"

        return(nullnames.correct & match.col.correct & length.correct &
                   sample.type.correct & aliquot.type.correct)
    })))
    eval(bquote(expect_equal(colname.correct, TRUE)))
})

test_that("Content of project_ids", {
    expect_equal(length(project_ids), 14551)
    expect_equal(sum(duplicated(names(project_ids))), 0)
    expect_equal(sum(grepl("^TCGA-", names(project_ids))), 11315)
    expect_equal(unique(nchar(grep("^TCGA-", names(project_ids),
                                   value = TRUE))), 12)
})

test_that("Content of surv.plot", {
    expect_equal(nrow(surv.plot), 12899)
    expect_equal(colnames(surv.plot), c("survivalEstimate", "id",
                                        "censored", "time"))
    expect_equal(sum(surv.plot$censored), 8712)
    expect_equal(sum(grepl("^TCGA-", rownames(surv.plot))), 11001)
})
