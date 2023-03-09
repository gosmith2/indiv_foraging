
theme_dark_black <- function(base_size=14, base_family="sans") {
   library(grid)
   library(ggthemes)
   (theme_foundation(base_size=base_size, base_family=base_family)
      + theme(plot.title = element_text(face = "bold", colour = '#ffffb3',
                                        size = rel(1.2), hjust = 0.5, margin = margin(0,0,20,0)),
              text = element_text(),
              panel.background = element_rect(colour = NA, fill = 'black'),
              plot.background = element_rect(colour = NA, fill = 'black'),
              panel.border = element_rect(colour = NA),
              axis.title = element_text(face = "bold",size = rel(1), colour = 'white'),
              axis.title.y = element_text(angle=90,vjust =2),
              axis.title.x = element_text(vjust = -0.2),
              axis.text = element_text(colour = 'white'),
              axis.line.x = element_line(colour="white"),
              axis.line.y = element_line(colour="white"),
              axis.ticks = element_line(colour="white"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "none",
              legend.background = element_rect(fill ='black'),
              legend.text = element_text(color = 'white'),
              legend.key = element_rect(colour = NA, fill = 'black'),
              ## legend.position = "bottom",
              legend.direction = "horizontal",
              legend.box = "vetical",
              legend.key.size= unit(0.5, "cm"),
              #legend.margin = unit(0, "cm"),
              legend.title = element_text(face="italic", colour = 'white'),
              plot.margin=unit(c(10,5,5,5),"mm"),
              strip.background=element_rect(colour="#2D3A4C",fill="black"),
              strip.text = element_text(face="bold", colour =
                                                         'white')
      ))
}






theme_ms <- function(base_size=14, base_family="sans") {
   library(grid)
   library(ggthemes)
   (theme_foundation(base_size=base_size, base_family=base_family)
      + theme(plot.title = element_text(face = "bold", colour = '#ffffb3',
                                        size = rel(1.2), hjust = 0.5,
                                        margin = margin(0,0,20,0)),
              text = element_text(),
              panel.background = element_rect(colour = NA, fill = 'white'),
              plot.background = element_rect(colour = NA, fill = 'white'),
              panel.border = element_rect(colour = NA),
              axis.title = element_text(face = "bold",size = rel(1), colour = 'black'),
              axis.title.y = element_text(angle=90,vjust =2),
              axis.title.x = element_text(vjust = -0.2),
              axis.text = element_text(colour = 'black'),
              axis.line.x = element_line(colour="black"),
              axis.line.y = element_line(colour="black"),
              axis.ticks = element_line(colour="black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              ##legend.position = "none",
              legend.background = element_rect(fill ='white'),
              legend.text = element_text(color = 'black'),
              legend.key = element_rect(colour = NA, fill = 'white'),
              legend.position = "right",
              legend.direction = "vertical",
              legend.box = "vetical",
              legend.key.size= unit(0.5, "cm"),
              #legend.margin = unit(0, "cm"),
              legend.title = element_text(face="italic", colour = 'black'),
              plot.margin=unit(c(10,5,5,5),"mm"),
              strip.background=element_rect(colour="#2D3A4C",fill="white"),
              strip.text = element_text(face="bold", colour =
                                                         'black')
      ))
}
