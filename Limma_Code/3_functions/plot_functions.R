#install.packages("ggthemes", lib = "C:/Program Files/R/R-3.5.2/library")

theme_Publication <- function(base_size=14) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  

}


# see https://community.rstudio.com/t/how-to-automatically-add-text-annotations-or-tags-outside-of-faceted-plots/13700/5

tag_facet2 <-  function(p, open=c("(",""), close = c(")","."),
                        tag_fun_right = function(i) letters[i],
                        x = c(0,0), y = c(0.5, 1),
                        hjust = c(0,0), vjust = c(0.5,1), 
                        fontface = c(2,2), ...){
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  nm <- names(gb$layout$facet$params$rows)
  
  
  tags_right <- paste0(open[2],tag_fun_right(unique(lay$ROW)),close[2])
  
  
  rl <- lapply(tags_right, grid::textGrob, x=x[2], y=y[2],
               hjust=hjust[2], vjust=vjust[2], gp=grid::gpar(fontface=fontface[2],...))
  
  
  g <- ggplot_gtable(gb)
 
  
  g <- gtable::gtable_add_cols(g, grid::unit(2,"line"), pos = -1)
  t <- unique(g$layout[grepl("panel",g$layout$name), "t"])
  g <- gtable::gtable_add_grob(g, grobs = rl, t=t, l=ncol(g))
  
  grid::grid.newpage()
  grid::grid.draw(g)
}


tag_facet1 <-  function(p, open=c("(",""), close = c(")","."),
                        tag_fun_top = function(i) letters[i],
                        tag_fun_right = utils::as.roman,
                        x = c(0,0), y = c(0.5, 1),
                        hjust = c(0,0), vjust = c(0.5,1), 
                        fontface = c(2,2), ...){
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  nm <- names(gb$layout$facet$params$rows)
  
  tags_top <- paste0(open[1],tag_fun_top(unique(lay$COL)),close[1])
  tags_right <- paste0(open[2],tag_fun_right(unique(lay$ROW)),close[2])
  
  tl <- lapply(tags_top, grid::textGrob, x=x[1], y=y[1],
               hjust=hjust[1], vjust=vjust[1], gp=grid::gpar(fontface=fontface[1], ...))
  
  
  g <- ggplot_gtable(gb)
  
  g <- gtable::gtable_add_rows(g, grid::unit(1,"line"), pos = 0)
  l <- unique(g$layout[grepl("panel",g$layout$name), "l"])
  g <- gtable::gtable_add_grob(g, grobs = tl, t=1, l=l)
  
  
  
  grid::grid.newpage()
  grid::grid.draw(g)
}

