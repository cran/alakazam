## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
library(alakazam)
library(dplyr)
library(airr)

db <- read_rearrangement(system.file("extdata", "example_quality.tsv", package="alakazam"))
fastq_file <- system.file("extdata", "example_quality.fastq", package="alakazam")

## -----------------------------------------------------------------------------
original_cols <- colnames(db)
db <- readFastqDb(db, fastq_file, style="both", quality_sequence=TRUE)
new_cols <- setdiff(colnames(db), original_cols)
db[,new_cols] %>% head()

## -----------------------------------------------------------------------------
quality <- getPositionQuality(db, sequence_id="sequence_id", 
                              sequence="sequence_alignment",
                              quality_num="quality_alignment_num")
head(quality)

## ---- fig.cap="Sequence quality per IMGT position for one sequence.", fig.asp=0.25----
min_pos <- min(quality$position)
max_pos <- max(quality$position)

ggplot(quality, aes(x=position,
                    y=quality_alignment_num,
                    color=nt)) +
  geom_point() +
  coord_cartesian(xlim=c(110,120)) +
  xlab("IMGT position") +
  ylab("Sequencing quality") +
  scale_fill_gradient(low = "light blue",  high = "dark red") +
  scale_x_continuous(breaks=c(min_pos:max_pos)) +
  alakazam::baseTheme()

## -----------------------------------------------------------------------------
db <- maskPositionsByQuality(db, min_quality=70,
                             sequence="sequence_alignment",
                             quality="quality_alignment_num")

