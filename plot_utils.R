plot.neg10logp <- function(pval.vec, alpha){
  # obtaining p values from spatial betareg Model
  pval.df <- as.data.frame(pval.vec)
  pval.df$Covariates <- rownames(pval.df) 
  rownames(pval.df)<- NULL
  pval.df <- pval.df %>% rename( "P_values"= "pval.vec")%>%arrange(P_values)
  pval.df$P_values <- -log10(pval.df$P_values)
  
  #betareg_P <- ggplot(pval.df, aes(P_values,Covariates))
  #betareg_P + geom_point() + theme(text = element_text(size=10)) + ggtitle("Spatial betareg P Values")
  
  # Descending Order spatial P values
  
  pval.plot <- pval.df %>%
    arrange(P_values) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
    mutate(name=factor(Covariates, levels=Covariates)) %>%   # This trick update the factor levels
    ggplot( aes(x=name, y=P_values)) +
    geom_segment( aes(xend=name, yend=0)) +
    geom_point( size=4, color="red") +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed")+
    # geom_hline(yintercept = -log10(0.1/20),linetype = "dotted")+
    annotate("text", x = c(2,2), y = c(-log10(alpha + 0.02)), 
             label = c("0.05", "0.1"), parse = T)+
    coord_flip() +
    theme_minimal() +
    labs(y = "-log10 p values",
         x = "", 
         title = "Beta regression parameters (P-values)",
         subtitle = "NSAA Data",
         caption = "Vertical lines indicate threshold of significance \n at 5% (dashed) and 10% (dotted)")
  
  return(pval.plot)
}

plot.beta <- function(beta.vec){
  # obtaining beta values from spatial betareg Model
  beta.df <- as.data.frame(beta.vec)
  beta.df$Covariates <- rownames(beta.df) 
  rownames(beta.df)<- NULL
  beta.df <- beta.df %>% rename( "betahat"= "beta.vec")%>%arrange(betahat)
  # beta.df$P_values <- -log10(beta.df$P_values)
  
  #betareg_P <- ggplot(beta.df, aes(P_values,Covariates))
  #betareg_P + geom_point() + theme(text = element_text(size=10)) + ggtitle("Spatial betareg P Values")
  
  # Descending Order spatial P values
  
  beta.plot <- beta.df %>%
    arrange(betahat) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
    mutate(name=factor(Covariates, levels=Covariates)) %>%   # This trick update the factor levels
    ggplot( aes(x=name, y=betahat)) +
    geom_segment( aes(xend=name, yend=0)) +
    geom_point( size=4, color="red") +
    geom_hline(yintercept = 0, linetype = "dashed")+
    # geom_hline(yintercept = -log10(0.1),linetype = "dotted")+
    annotate("text", x = c(2), y = c(0), 
             label = c("0"), parse = T)+
    coord_flip() +
    theme_minimal() +
    labs(y = "Z-scores",
         x = "", 
         title = "Beta regression parameters (Z-scores)",
         subtitle = "NSAA Data",
         caption = "Vertical line at 0")
  
  return(beta.plot)
}


plotTheme <- function() {
  theme(
    plot.title = element_text(size = 14, family = "sans", face = "plain", hjust = 0),
    plot.subtitle=element_text(size = 11, family = "sans", hjust = 0),
    plot.caption=element_text(size = 10, family = "sans", face = "italic", hjust = 0), 
    axis.title.x = element_text(size = 10, family = "sans", face = "plain", hjust = 1, vjust = -0.5),
    axis.title.y = element_text(size = 10, family = "sans", face = "plain", hjust = 1, vjust = 1),
    axis.text = element_text(size = 9, family = "sans", face = "plain"),
    panel.background = element_blank(),
    panel.grid.minor = element_line(colour = "gray"),
    panel.grid.major = element_line(colour = "gray"),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 10, family = "sans"),
    legend.text = element_text(size = 9, family = "sans"),
    axis.line = element_blank()
  )
}

