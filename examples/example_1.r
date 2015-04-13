region <- list(NULL, c("REGION_1", "REGION_2"))

raw.data <- load.data.snapp("examples/example_snp.snapp", 
                            "examples/example_coordinates.csv", 
                            sep=',', ploidy=ploidy)                                

pop <- make.pop(raw.data, ploidy)
psi <- get.all.psi(pop)

res <- run.regions(region=region, pop=pop, psi=psi, xlen=10,ylen=20)

summary(res)
plot(res)
