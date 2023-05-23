requireNamespace('here')
requireNamespace('ggplot2')


# Read data
data <- read.csv(here::here('results','public_tcrs','public_counts.tsv'), sep='\t')


# Create plot
public_plot<-ggplot2::ggplot(data=data, ggplot2::aes(x=nr_volunteers, y=CDR3_beta, fill=epitope,label=epitope)) +
  ggplot2::geom_bar(stat="identity", position = "dodge") +
  ggplot2::geom_text(ggplot2::aes(label = CDR3_beta),  vjust = -0.4, size=5, position = ggplot2::position_dodge(0.9)) +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size=12),
        panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor.x = ggplot2::element_blank())+
  ggplot2::ylab("Number of public CDR3s shared by different volunteers") +
  ggplot2::xlab("Number of volunteers having the shared CDR3 sequence")+
  ggplot2::labs(title = "Publicness of WT1 CDR3 beta sequences")

# Show plot
public_plot

# Save plot
pdf(file=here::here('results','public_tcrs','public_tcrs_plot.pdf'))
public_plot
dev.off()



