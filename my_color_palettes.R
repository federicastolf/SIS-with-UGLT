library(fields)

my_color_palettes = function(theme = "Instagram", reverse = FALSE){
  
  if(tolower(theme) == "unipd"){
    
    set_sparse = c("#9B0014","#9B0014","brown3","firebrick2","indianred1", 
                   "rosybrown1","white","azure",
                   "lightsteelblue2","lightskyblue2","steelblue2","#0093D5", "#0093D5")
    
    set_heatmap = c("#484F59","snow4","mistyrose4",
                    "mistyrose3","rosybrown2",
                    "indianred3","brown3","#9B0014")
    
    set_scatter = c("#0093D5", "lightskyblue2", "turquoise2",
                    "lightsteelblue4","#484F59","mistyrose4",
                    "orangered1","rosybrown2","#9B0014")
    
  }else{
    if(tolower(theme) != "instagram"){
      warning("You can choose between 'Instagram' theme or 'Unipd' theme")
    }
    
    set_sparse = c("orangered2","orangered2","orange2","goldenrod2","gold",
                   "lightgoldenrod","lightyellow","white","azure","lightsteelblue2",
                   "slateblue1","slateblue3","royalblue3","royalblue4", "royalblue4")
      
    set_heatmap = c("yellow", "gold","goldenrod2","orange2",
                    "orangered2","red2","violetred2","magenta3",
                    "mediumorchid3","purple3","slateblue3","royalblue4")
    
    set_scatter = c("gold","orange2","red2",
                    "hotpink1","magenta3", "purple3",
                    "aquamarine2","dodgerblue2", "royalblue4")
  }
  
  if(reverse){
    return(list( sparse = colorRampPalette(rev(set_sparse)),
                 heatmap = colorRampPalette(rev(set_heatmap)),
                 scatter = colorRampPalette(rev(set_scatter))))
  }else{
    return(list( sparse = colorRampPalette(set_sparse),
                 heatmap = colorRampPalette(set_heatmap),
                 scatter = colorRampPalette(set_scatter)))
  }

}


display_my_palettes = function(ncolors = 10){
  
  par(mfrow=c(2,3))
  
  dat = as.matrix(iris[,1:4])
  sc_dat = scale(dat)
  sc_dat[which(abs(sc_dat)<0.5)] = 0
  km = kmeans(iris[,1:2], centers = ncolors)
  
  pal1 = my_color_palettes("Instagram")
  pal2 = my_color_palettes("Unipd")
  
  ncolors_star = ifelse(ncolors%%2==0, ncolors+1, ncolors)
  image.plot(t(sc_dat), axes = F, col = pal1$sparse(ncolors_star),
             breaks = c(seq(min(sc_dat),-0.5,length = floor(ncolors/2)),
                        -0.25,0.25, seq(0.5,max(sc_dat),length = floor(ncolors/2))),
             main = "Sparse"
             )
  mtext("Instagram palette", side=2, adj=0, cex.lab=2,las=3)
  
  image.plot(t(dat), axes = F, col = pal1$heatmap(ncolors), 
             main = "Heatmap")
  
  plot(iris[,1:2], xlab ="", col = pal1$scatter(ncolors)[km$cluster], pch = 20, cex = 2,
       main = "Scatter")
  
  
  image.plot(t(sc_dat), axes = F, col = pal2$sparse(ncolors_star),
             breaks = c(seq(min(sc_dat),-0.5,length = floor(ncolors/2)),
                        -0.25,0.25, seq(0.5,max(sc_dat),length = floor(ncolors/2)))
             #, main = "Sparse"
  )
  mtext("Unipd palette", side=2, adj=0, cex.lab=2,las=3)
  
  image.plot(t(dat), axes = F, col = pal2$heatmap(ncolors))
  
  plot(iris[,1:2], xlab ="",
       col = pal2$scatter(ncolors)[km$cluster], pch = 20, cex = 2)
  
  par(mfrow=c(1,1))
}


my_colors = function(ncolors=10,theme="Instagram", type="Heatmap",
                     reverse = FALSE){
  
  f_col = my_color_palettes(theme = theme, reverse = reverse)
  if(tolower(type) =="sparse"){
    return(f_col$sparse(ncolors))
  }else if(tolower(type) =="heatmap"){
    return(f_col$heatmap(ncolors))
  }else if(tolower(type) =="scatter"){
    return(f_col$scatter(ncolors))
  }else{
    warning("Palettes type are 'Heatmap','Sparse', or 'Scatter'.")
    return(f_col$heatmap(ncolors))
  }
  
}



display_my_palettes(10)
