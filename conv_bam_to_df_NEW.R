#this code should get a bam file and then it converts it into readable dataframe for further analysis

#how to launch, an example:
#R CMD BATCH --vanilla "--args induction_res_FragSeq/mapped.bam induction_res_FragSeq/mapped_bam.csv" conv_bam_to_df_NEW.R

#read arguments from the command line
myargs = commandArgs(trailingOnly=TRUE)
path_to_bam = myargs[1]
path_to_csv = myargs[2]


#load library
library(Rsamtools)
#read in entire BAM file

bam <- scanBam(path_to_bam)

#names of the BAM fields
names(bam[[1]])
#distribution of BAM flags
table(bam[[1]]$flag)
.unlist <- function (x){
    ## do.call(c, ...) coerces factor to integer, which is undesired
    x1 <- x[[1L]]
    if (is.factor(x1)){
       structure(unlist(x), class = "factor", levels = levels(x1))
    } else {
       do.call(c, x)
    }
 }
#store names of BAM fields
bam_field <- names(bam[[1]])

#go through each BAM field and unlist
list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))
               
#store as data frame
bam_df <- do.call("DataFrame", list)
names(bam_df) <- bam_field
dim(bam_df)
#View(bam_df)
write.table(bam_df, path_to_csv )
#head(bam_df)

