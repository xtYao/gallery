alt = dir("~/ALT/JaBbA/", "jabba.simple.rds", full.names=T, all.files=T, recursive=T)
alt = data.table(fn = alt,
                 id = dir("~/ALT/JaBbA/"))
setkey(alt, "id")

jabba = alt[, lapply(fn, readRDS)]

saveRDS(alt, "alt.jabba.simple.rds")
