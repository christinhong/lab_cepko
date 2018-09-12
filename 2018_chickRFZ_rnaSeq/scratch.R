#### Add lanes together ####
# Since the lanes look good, may as well merge those rows together to simplify the analysis and visualization
## from https://stackoverflow.com/questions/32462640/sum-rows-in-data-frame-based-on-two-columns

# Remove Lane and Sample columns (note that the rownames for this data table are from the Sample column)
fc2 <- dplyr::select(fc, -Lane, -Sample)

fc2[1:12,1:8]


# Sum columns after grouping on Tissue, Retina, and Batch
fc3 <- data.table::setDT(fc2)[, lapply(.SD, sum), by = c("Tissue", "Retina", "Batch")]
# Went with data.table instead of dplyr because data.table maintains significant figures while dplyr rounds. data.table also runs notably faster.

fc3[1:3,1:8]
tail(fc3[,1:8])
nrow(fc3) # OK, looks good.


