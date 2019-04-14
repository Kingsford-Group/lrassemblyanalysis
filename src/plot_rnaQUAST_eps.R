#
# plot_rnaQUAST_eps.R
#
# Laura Tung
#
# Plotting box plots of rnaQUAST results in EPS.
#
# IMPORTANT:
#	Please replace <DATASETS_LOCATION> in the path_vector (at the end of this file) 
#       with your actual datasets location (full-path directory) before running this code.
#
# Usage:
#        Rscript plot_rnaQUAST_eps.R <Organism>
#
# <Organism>: human or mouse
# Run process_rnaQUAST_results.sh first before running this code.

plot_box_plots <- function(path_vector, organism) {

   Scallop_isoform_0_50 <- c()
   Stringtie_isoform_0_50 <- c()
   Isoseq_isoform_0_50 <- c()

   Scallop_isoform_50_75 <- c()
   Stringtie_isoform_50_75 <- c()
   Isoseq_isoform_50_75 <- c()

   Scallop_isoform_75_95 <- c()
   Stringtie_isoform_75_95 <- c()
   Isoseq_isoform_75_95 <- c()

   Scallop_isoform_95_100 <- c()
   Stringtie_isoform_95_100 <- c()
   Isoseq_isoform_95_100 <- c()

   Scallop_matched_0_50 <- c()
   Stringtie_matched_0_50 <- c()
   Isoseq_matched_0_50 <- c()

   Scallop_matched_50_75 <- c()
   Stringtie_matched_50_75 <- c()
   Isoseq_matched_50_75 <- c()

   Scallop_matched_75_95 <- c()
   Stringtie_matched_75_95 <- c()
   Isoseq_matched_75_95 <- c()

   Scallop_matched_95_100 <- c()
   Stringtie_matched_95_100 <- c()
   Isoseq_matched_95_100 <- c()

   for (i in 1:length(path_vector)) {
        rna_data <- read.table(paste0(path_vector[i], "/rnaQUAST_bins"), header = TRUE)

       Scallop_isoform_0_50[i] <- rna_data$Scallop[1]       
       Stringtie_isoform_0_50[i] <- rna_data$StringTie[1]
       Isoseq_isoform_0_50[i] <- rna_data$Isoseq[1]

       Scallop_isoform_50_75[i] <- rna_data$Scallop[2]
       Stringtie_isoform_50_75[i] <- rna_data$StringTie[2]
       Isoseq_isoform_50_75[i] <- rna_data$Isoseq[2]

       Scallop_isoform_75_95[i] <- rna_data$Scallop[3]
       Stringtie_isoform_75_95[i] <- rna_data$StringTie[3]
       Isoseq_isoform_75_95[i] <- rna_data$Isoseq[3]

       Scallop_isoform_95_100[i] <- rna_data$Scallop[4]
       Stringtie_isoform_95_100[i] <- rna_data$StringTie[4]
       Isoseq_isoform_95_100[i] <- rna_data$Isoseq[4]

       Scallop_matched_0_50[i] <- rna_data$Scallop[5]
       Stringtie_matched_0_50[i] <- rna_data$StringTie[5]
       Isoseq_matched_0_50[i] <- rna_data$Isoseq[5]

       Scallop_matched_50_75[i] <- rna_data$Scallop[6]
       Stringtie_matched_50_75[i] <- rna_data$StringTie[6]
       Isoseq_matched_50_75[i] <- rna_data$Isoseq[6]

       Scallop_matched_75_95[i] <- rna_data$Scallop[7]
       Stringtie_matched_75_95[i] <- rna_data$StringTie[7]
       Isoseq_matched_75_95[i] <- rna_data$Isoseq[7]

       Scallop_matched_95_100[i] <- rna_data$Scallop[8]
       Stringtie_matched_95_100[i] <- rna_data$StringTie[8]
       Isoseq_matched_95_100[i] <- rna_data$Isoseq[8]
   }

   print(length(Scallop_isoform_0_50))

   # plot box plot for assembled isoforms
   title <- paste0("Box plots of assembled isoforms for Scallop-LR, StringTie, Iso-Seq (", organism, " datasets)")
   xlabel <- "Range of assembled fraction (x-y%) of isoform" 
   ylabel <- "Number of x-y% assembled isoforms"
   
   setEPS()
   postscript(paste0("rnaQUAST_box_plots_", organism, "_assembled_isoforms.eps"), paper="special", horizontal=FALSE, width=11, height=9)
   boxplot(Scallop_isoform_0_50, Stringtie_isoform_0_50, Isoseq_isoform_0_50, Scallop_isoform_50_75, Stringtie_isoform_50_75, Isoseq_isoform_50_75, Scallop_isoform_75_95, Stringtie_isoform_75_95, Isoseq_isoform_75_95, Scallop_isoform_95_100, Stringtie_isoform_95_100, Isoseq_isoform_95_100, main=title, xlab=xlabel, ylab=ylabel, at=c(1,2,3,6,7,8,11,12,13,16,17,18), names=c("", "0-50%", "", "", "50-75%", "", "", "75-95%", "", "", "95-100%", ""), col=c("green", "red", "blue"), border="black")
   legend(x="top", fill=c("green", "red", "blue"), legend=c("Scallop-LR", "StringTie", "Iso-Seq"), col=c("green", "red", "blue"))
   dev.off()

   # plot box plot for matched transcripts
   title <- paste0("Box plots of matched transcripts for Scallop-LR, StringTie, Iso-Seq (", organism, " datasets)")
   xlabel <- "Range of matched fraction (x-y%) of transcript"
   ylabel <- "Number of x-y% matched transcripts"

   setEPS()
   postscript(paste0("rnaQUAST_box_plots_", organism, "_matched_transcripts.eps"), paper="special", horizontal=FALSE, width=11, height=9)
   boxplot(Scallop_matched_0_50, Stringtie_matched_0_50, Isoseq_matched_0_50, Scallop_matched_50_75, Stringtie_matched_50_75, Isoseq_matched_50_75, Scallop_matched_75_95, Stringtie_matched_75_95, Isoseq_matched_75_95, Scallop_matched_95_100, Stringtie_matched_95_100, Isoseq_matched_95_100, main=title, xlab=xlabel, ylab=ylabel, at=c(1,2,3,6,7,8,11,12,13,16,17,18), names=c("", "0-50%", "", "", "50-75%", "", "", "75-95%", "", "", "95-100%", ""), col=c("green", "red", "blue"), border="black")
   legend(x="top", fill=c("green", "red", "blue"), legend=c("Scallop-LR", "StringTie", "Iso-Seq"), col=c("green", "red", "blue"))
   dev.off()

}


plot_mean_box_plots <- function(path_vector, organism) {

   Scallop_mean_isoform <- c()
   Stringtie_mean_isoform <- c()
   Isoseq_mean_isoform <- c()

   Scallop_mean_transcript <- c()
   Stringtie_mean_transcript <- c()
   Isoseq_mean_transcript <- c()

   for (i in 1:length(path_vector)) {
        mean_isoform <- read.table(paste0(path_vector[i], "/rnaQUAST_output/mean_isoform"), header = FALSE)
        Scallop_mean_isoform[i] <- mean_isoform[1, 4]
        Stringtie_mean_isoform[i] <- mean_isoform[1, 5]
        Isoseq_mean_isoform[i]  <- mean_isoform[1, 6]

        mean_transcript <- read.table(paste0(path_vector[i], "/rnaQUAST_output/mean_transcript"), header = FALSE)
        Scallop_mean_transcript[i] <- mean_transcript[1, 6]
        Stringtie_mean_transcript[i] <- mean_transcript[1, 7]
        Isoseq_mean_transcript[i] <- mean_transcript[1, 8]
   }

   print(length(Scallop_mean_isoform))

   # plot box plot for mean fractions of isoform and transcript
   title <- paste0("Box plots of mean fractions of isoform assembled & transcript matched for Scallop-LR, StringTie, Iso-Seq (", organism, " data)")
   ylabel <- "Mean fractions"

   setEPS()
   postscript(paste0("rnaQUAST_box_plots_", organism, "_mean_fractions.eps"), paper="special", horizontal=FALSE, width=11.5, height=9)
   boxplot(Scallop_mean_isoform, Stringtie_mean_isoform, Isoseq_mean_isoform, Scallop_mean_transcript, Stringtie_mean_transcript, Isoseq_mean_transcript, main=title, ylab=ylabel, at=c(1,2,3,7,8,9), names=c("", "Mean isoform assembly", "", "                                          Mean fraction of transcript matched", "", ""), col=c("green", "red", "blue"), border="black")
   legend(x="top", fill=c("green", "red", "blue"), legend=c("Scallop-LR", "StringTie", "Iso-Seq"), col=c("green", "red", "blue"))
   dev.off()

}


args = commandArgs(trailingOnly=TRUE)

organism <- args[1]

if (organism == "human") {
    cat("human\n")
    path_vector <- c("<DATASETS_LOCATION>/human/ERP015321/BioSamples/SAMN00001694/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/human/ERP015321/BioSamples/SAMN00001695/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/human/ERP015321/BioSamples/SAMN00001696/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/human/ERP015321/BioSamples/SAMN00006465/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/human/ERP015321/BioSamples/SAMN00006466/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/human/ERP015321/BioSamples/SAMN00006467/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/human/ERP015321/BioSamples/SAMN00006579/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/human/ERP015321/BioSamples/SAMN00006580/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/human/ERP015321/BioSamples/SAMN00006581/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/human/SRP126849/BioSamples/SAMN08182059/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/human/SRP126849/BioSamples/SAMN08182060/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/human/SRP071928/BioSamples/SAMN04563763/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/human/SRP068953/BioSamples/SAMN04169050/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/human/SRP098984/BioSamples/SAMN07611993/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/human/SRP065930/BioSamples/SAMN04251426_1/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/human/SRP065930/BioSamples/SAMN04251426_2/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/human/SRP065930/BioSamples/SAMN04251426_3/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/human/SRP065930/BioSamples/SAMN04251426_4/ccs_flnc_and_nfl/minimap2")
} else {
    cat("mouse\n")
    path_vector <- c("<DATASETS_LOCATION>/mouse/ERP010189/BioSamples/SAMEA3374575/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/mouse/ERP010189/BioSamples/SAMEA3374576/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/mouse/ERP010189/BioSamples/SAMEA3374577/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/mouse/ERP010189/BioSamples/SAMEA3374578/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/mouse/ERP010189/BioSamples/SAMEA3374579/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/mouse/ERP010189/BioSamples/SAMEA3374580/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/mouse/ERP010189/BioSamples/SAMEA3374581/ccs_flnc_and_nfl/minimap2",
                     "<DATASETS_LOCATION>/mouse/ERP010189/BioSamples/SAMEA3374582/ccs_flnc_and_nfl/minimap2")
}

plot_box_plots(path_vector, organism)

plot_mean_box_plots(path_vector, organism)



