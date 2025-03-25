require(tidyverse)

setwd("D:/BMC/Box/Helping_with_analysis/Carina/TAPS_summary/Extract_Segments")
load('allRegions.Rdata')
# Probe specific datasets
# alf : allele frequency, Log2, regs: Approximated CN regions (exhaustive)
load('TAPS_plot_output.Rdata')

# calculating mBAFs of each probe
alf$mBAF <- 0.5+abs(0.5 - alf$Value)

{#===============DNAcopy===============#
# regs <- regs %>%
#   rowwise() %>%
#   mutate(Low_mBAF = {
#     mBAF_values <- alf %>%
#       filter(Start >= start & Start <= end) %>%
#       pull(mBAF)
#     mean(mBAF_values[mBAF_values < 0.9])
#   })
# 
# regs <- regs %>%
#   rowwise() %>%
#   mutate(probe_density = {
#     probe_number <- alf %>%
#       filter(Start >= start & Start <= end) %>%
#       pull(mBAF)
#       (end-start)/length(probe_number)
#   })
}
#===============DNAcopy : Segment Prediction===============#
require(DNAcopy)
alf_noHomo <- alf %>%  filter (Value < 0.95 & Value > 0.05)
segment_by_AF <- segment(CNA(alf_noHomo$mBAF, alf_noHomo$Chromosome, alf_noHomo$Start))
segment_by_LR <- segment(CNA(Log2$Value, Log2$Chromosome, Log2$Start))

chr.analysis='chr1'

#===============functions for Visualization and Segment Prediction===============#
SNP.stat <- function(base.SNP.data, segment_data=segment_data, loc.start = loc.start, loc.end = loc.end){
#'*base.SNP.data is the plotting base data where the probe stats are stored. It can be ALF, Log2 etc*
  chr.dat <<- base.SNP.data %>% 
  filter(Chromosome == chr.analysis)
  nr.bins <<- nrow(chr.dat)
  chr.length.mb <<- round((max(chr.dat$End) - min (chr.dat$Start))/1e6,1)
  bp.start.lim <- min(chr.dat$Start)
  bp.end.lim <- max(chr.dat$End)
  
#'*segment_data is where segment predictions are saved*
#'*loc.start and loc.end are starting and ending base pair positions*
  test_segment <- segment_data %>% 
  filter(chrom == chr.analysis & loc.start > bp.start.lim-1 & loc.end < bp.end.lim+1) %>%
    mutate(orig.loc.start = loc.start, orig.loc.end = loc.end) %>% 
  mutate(loc.start = loc.start / 1e6, loc.end = loc.end / 1e6)  %>% 
  mutate (loc.end = loc.end * (nr.bins / chr.length.mb)) 
  test_segment$loc.start[-1] = test_segment$loc.start[-1] * (nr.bins / chr.length.mb)
  test_segment <<- test_segment
}
Allele.freq.plot <- function(segment_data=segment_data){
  SNP.stat(base.SNP.data = alf, segment_data)
  plot(chr.dat$mBAF,ylim = c(0,1), col = alpha('grey', 0.6), cex = 1, xaxt='n', xlab='Basepair', ylab = 'Allele frequency', xaxs="i")
  axis(1, at=1:nr.bins, labels = round(seq(0, chr.length.mb, by = chr.length.mb/nr.bins),0)[-1])
  segments(test_segment$loc.start,test_segment$seg.mean,test_segment$loc.end,test_segment$seg.mean, col='darkorchid', lwd=3)
  alf_segments<<-test_segment
}
Log.ratio.plot <- function(segment_data=segment_data){
  SNP.stat(base.SNP.data = Log2, segment_data)
  plot(chr.dat$Value,ylim = c(-1.5,1.5), col = alpha('grey', 0.6), cex = 1, xaxt='n', xlab='Basepair', ylab = 'Log ratio', xaxs="i")
  axis(1, at=1:nr.bins, labels = round(seq(0, chr.length.mb, by = chr.length.mb/nr.bins),0)[-1])
  segments(test_segment$loc.start,test_segment$seg.mean,test_segment$loc.end,test_segment$seg.mean, col='coral3', lwd=3)
  Log2_segments<<-test_segment
}
Merge_segment <- function(log2_mean, seg_length_kb) {
  n <- length(log2_mean)
  merged_segment <- numeric(n)
  merged_segment[1] <- 1
  for (i in 2:n) {
    if (abs(log2_mean[i] - log2_mean[i-1]) <= 0.2 | seg_length_kb[i] <= 10) {
      merged_segment[i] <- merged_segment[i-1]
    } else if (seg_length_kb[i-1] <= 10) {
      merged_segment[i] <- merged_segment[i-1]
    } else {
      merged_segment[i] <- merged_segment[i-1] + 1
    }
  }
  return(merged_segment)
}

#===============Visualization===============#
par(mfrow=c(5,1))
Allele.freq.plot(segment_data = segment_by_AF[['output']])
Log.ratio.plot(segment_data = segment_by_LR[['output']])

#===============Optional : overlay TAPS Segment Prediction===============#
{chr.dat <<- Log2 %>% 
    filter(Chromosome == chr.analysis)
  nr.bins <<- nrow(chr.dat)
  chr.length.mb <<- round((max(chr.dat$End) - min (chr.dat$Start))/1e6,1)
  bp.start.lim <- min(chr.dat$Start)
  bp.end.lim <- max(chr.dat$End)
  
  test_segment <- segments %>% 
    filter(Chromosome == chr.analysis & Start > bp.start.lim-1 & End < bp.end.lim+1) %>%
    mutate(loc.start = Start / 1e6, loc.end = End / 1e6)  %>% 
    mutate (loc.end = loc.end * (nr.bins / chr.length.mb)) 
  test_segment$loc.start[-1] = test_segment$loc.start[-1] * (nr.bins / chr.length.mb)
  test_segment <<- test_segment
segments(test_segment$loc.start,test_segment$Value,test_segment$loc.end,test_segment$Value, col='aquamarine4', lwd=3, lty=2)
}
rm(list = c("nr.bins", "chr.length.mb", "chr.dat", "test_segment"))

#===============Segment Curation based on Log2 ratio after Visualization===============#
Log2_segments <- Log2_segments %>% mutate(seg_length_kb = round((orig.loc.end - orig.loc.start)/1e3,1)) %>% 
  rowwise() %>%
  mutate(alf_subset = list(subset(alf_noHomo, Chromosome == chr.analysis & 
                                    Start > orig.loc.start & End < orig.loc.end))) %>%
  mutate(log2_subset = list(subset(Log2, Chromosome == chr.analysis & 
                                     Start > orig.loc.start & End < orig.loc.end))) %>%
  mutate(alf_variance = round(var(alf_subset$Value*10),1)) %>% 
  mutate(log2_mean = round(mean(log2_subset$Value),3)) %>% 
  select(-alf_subset,-log2_subset)%>%
  mutate(alf_variance = ifelse(is.na(alf_variance), 0, alf_variance))%>%
  mutate(log2_mean = ifelse(is.na(log2_mean), 0, log2_mean)) 

  Log2_segments$merged_segment <- Merge_segment(Log2_segments$log2_mean, Log2_segments$seg_length_kb)

#'*We are ignoring but not removing segments that are less than 10kb long*
#'*Most amplicons are on an avg few hundred to about 5kb long*


Segment_summary_LR <- Log2_segments %>% group_by(merged_segment) %>% 
  summarise(chrom = unique(chrom),
            loc.start = min(orig.loc.start),
            loc.end = max(orig.loc.end),
            length = (loc.end - loc.start),
            num.mark = sum(num.mark),
            uniq.segs = length(seg.mean),
            seg.mean = round(median(log2_mean),3))
Log.ratio.plot(segment_data = Segment_summary_LR)

#===============Segment Curation based on Allele freq ratio after Visualization===============#
alf_noHomo <- alf %>%  
  filter (Value < 0.9 & Value > 0.1) %>% 
  filter(Chromosome == chr.analysis) 
alf_noHomo <- alf_noHomo %>%
  mutate(bp_start := Start) %>%
  mutate(bp_end := Start+100000) %>% 
  rowwise() %>%
  mutate(noHomo_AF = {
    mBAF_values <- alf_noHomo %>%
      filter(Start >= bp_start & End <= bp_end) %>%
      pull(mBAF)
    mean(mBAF_values)
    probe_number <- length(mBAF_values)
  })

alf_noHomo <- alf_noHomo %>% 
  mutate(merged_segment = 1)
{
# for (i in 2:nrow(alf_noHomo)) {
#   # If the absolute difference between the current and previous mBAF values is less than 0.3,
#   # keep the same segments value, otherwise increment it by 1
#   if (abs(1 - alf_noHomo$mBAF[i] / alf_noHomo$mBAF[i-1]) < 0.3) {
#     alf_noHomo$merged_segment[i] <- alf_noHomo$merged_segment[i-1]
#   } else {
#     alf_noHomo$merged_segment[i] <- alf_noHomo$merged_segment[i-1] + 1
#   }
# }
  }
for (i in 2:nrow(alf_noHomo)) {
  # If the absolute difference between the current and previous mBAF values is less than 0.3,
  # keep the same segments value, otherwise increment it by 1
  if (abs(1 - alf_noHomo$noHomo_AF[i] / alf_noHomo$noHomo_AF[i-1]) < 1.6) {
    alf_noHomo$merged_segment[i] <- alf_noHomo$merged_segment[i-1]
  } else {
    alf_noHomo$merged_segment[i] <- alf_noHomo$merged_segment[i-1] + 1
  }
}
Segment_summary_ALF <- alf_noHomo %>% group_by(merged_segment) %>% 
  summarise(chrom = unique(Chromosome),
            loc.start = min(Start),
            loc.end = max(End),
            length = (loc.end - loc.start),
            uniq.segs = length(mBAF),
            Value.mean = round(median(Value),3),
            seg.mean = round(median(mBAF),3)) %>% 
  filter(uniq.segs > 5)
Allele.freq.plot(segment_data = Segment_summary_ALF)
