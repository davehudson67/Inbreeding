library(related)


badgers <- readgenotypedata("Gene_dataEDIT.txt")

QGestimator <- coancestry(badgers$gdata, quellergt = 1)

sim <- familysim(badgers$freqs, 100)
output <- coancestry(sim, wang = 1)
simrel <- cleanuprvals(output$relatedness, 100)

relvalues <- simrel[, 6]

label1 <- rep("PO", 100)
label2 <- rep("FULL", 100)
label3 <- rep("Half", 100)
label4 <- rep("Unrelated", 100)

labels <- c(label1, label2, label3, label4)

plot(as.factor(labels), relvalues, ylab = "Relatedness Value", xlab = "Relatedness")
qplot(as.factor(labels), relvalues, geom = "boxplot", ylab = "Relatedness Value", xlab = "Relatedness")

Relationship <- labels
newdata <- as.data.frame(cbind(Relationship, relvalues))
newdata$relvalues <- as.numeric(as.character(newdata$relvalues))
qplot(relvalues, ..density.., data = newdata, geom = "density", colour = as.factor(Relationship), ylab = "Relatedness Value", xlab = "Relatedness")


compareestimators(badgers, 100)


hist(CH$hom)
ggplot(CH) +
  geom_histogram(bins = 14, colour = "black", fill = "white", aes(x = hom)) +
  geom_vline(xintercept = 0.3181818, colour = "red" ) +
  labs(x = "Multi-locus homozygosity", y = "Count") +
  theme(text = element_text(size = 20))

ggplot(CH) + 
  geom_point(aes(x = ))

median(CH$hom)
IQR(CH$hom)
max(CH$hom)
