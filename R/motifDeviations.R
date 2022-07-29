

### ChromVar-inspired analysis for how specific peaksets use less or more 
### Designed to be used with scMACs normalization protocol. Use with other normalization schemes may not control for technical bias. 



######## getMotifCoverage: Function that takes in Sample-specific coverage data and finds 
## the average coverage for each motif location within a list of motif locations
## @covFiles: GRangesList of coverage for each sample
## @MotifLocations: GRangesList of Motif locations
## @numCores: number of cores to parallelize over. 

getMotifCoverages <- function(covFiles, MotifLocations, numCores=1){
    
  allMotifs <- stack(MotifLocations)
  counts <- mclapply(covFiles, function(x){ 
                	x %>% plyranges::mutate(NewScore = score) %>%
                    plyranges::join_overlap_intersect(allMotifs) %>%
                    plyranges::mutate(WeightedScore = NewScore*width(.)) %>%
		    plyranges::group_by(name) %>%
                    plyranges::reduce_ranges(score = mean(WeightedScore)) 
    
    	}, mc.cores = numCores)
 return(counts)
}

    






normalized_wilcoxon <- function(data, group, point.mass=0, test='wilcoxon'){
    # Function for calculating two-part statistics
    Index1 <- c(group==1)
    Group1 <- data[Index1]
    Group0 <- data[!Index1]
    n1 <- length(Group1)
    n2 <- length(Group0)
    obs <- c(n1, n2)
    success <- c(sum(Group1!=point.mass), sum(Group0!=point.mass))
    pointmass <- obs-success
    if (sum(success)==0) {
        T2 <- 0
        B2 <- 0
        }
    else if ((success[1]==0)|(success[2]==0)) {
        T2 <- 0
        B2 <- prop.test(pointmass, obs)$statistic
    }
    else if ((success[1]==1)|(success[2]==1)){
        T2 <- 0
        B2 <- prop.test(pointmass, obs)$statistic
    }
    else {
        uniq1 <- length(unique(Group1[Group1!=point.mass]))
        uniq2 <- length(unique(Group0[Group0!=point.mass]))
    if ((uniq1 < 2) & (uniq2 < 2)){
        T2 <- 0
    if (sum(pointmass)==0)
        B2 <- 0
    else
        B2 <- prop.test(pointmass, obs)$statistic
    }
    else if (sum(pointmass)==0){
        B2 <- 0
    if (test=="t.test")
        T2 <- t.test(data~group)$statistic
    if (test=="wilcoxon") {
        W <- wilcox.test(data~group, exact=FALSE)$statistic
        mu <- (n1*n2)/2
        sigma <- sqrt((n1*n2*(n1+n2+1))/12)
        T2 <- ((W-mu)-0.5)/sigma
    }
    }
        else {
            B2 <- prop.test(pointmass, obs)$statistic
            contIndex <- data!=point.mass
            cont <- data[contIndex]
            cGroup <- group[contIndex]
            n1c <- sum(cGroup==1)
            n2c <- sum(cGroup==0)
            if (test=="t.test")
                T2 <- t.test(cont~cGroup)$statistic
            if (test=="wilcoxon") {
                W <- wilcox.test(cont~cGroup, exact=FALSE)$statistic
                mu <- (n1c*n2c)/2
                sigma <- sqrt((n1c*n2c*(n1c+n2c+1))/12)
                T2 <- ((W-mu)-0.5)/sigma
            }
        }
    }
    return(T2)
}
Example
a = sapply(1:100000, function(x) 
    { data = rpois(n=100, lambda=20)
    group = c(rep(1,50),rep(0,50))
    normalized_wilcoxon(data,group,test='wilcoxon')
     })
png('tmp.png')
hist(a, breaks=100)
dev.off()