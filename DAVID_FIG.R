# Tips: only English in file path
###################
library(ggplot2)
library(scatterpie)
library(stringr)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Read parameters
parameter <- read.csv('set parameters.csv', row.names =1)

IS_log_change <- parameter[1,1]
fc <- parameter[2,1]
p_threshold <- parameter[3,1]
p_threshold_path <- parameter[4,1]
pie_radius <- parameter[5,1]
fontsize <- parameter[6,1]

# Read DAVID info
FILE <- list.files(pattern = '.txt')
data_all <- read.delim(FILE, header = FALSE, sep = "\t")

catogory_list = unique(data_all$V1)
catogory_list = catogory_list[!str_detect(catogory_list, 'Annotation|Category')]

for (aa in 1: length(catogory_list)) {
  choosedata <- data_all[,1]== catogory_list[aa]
  
  if (sum(choosedata) > 1) {
    data <- data_all[data_all[,1]==catogory_list[aa],]
    data <- unique(data)
    if (sum(data$V13 <= 0.05) >0) {
      # Read RNA-seq info
      data2 <- read.csv('data.csv')
      colnames(data2)[1] <- 'Gene_name'
      colnames(data2)[2] <- 'fold_change'
      colnames(data2)[3] <- 'pvalue'
      data2$Gene_name <- as.character(data2$Gene_name)
      data2$fold_change <- as.numeric(as.character(data2$fold_change))
      data2$pvalue <- as.numeric(as.character(data2$pvalue))
      data2[is.na(data2)] <- 1
      
      
      if (IS_log_change == 0) {
        data2$fold_change <- log2(data2$fold_change)
      }
      
      
      PATH_NUMBER <- dim(data)[1]
      
      for (a in 1:PATH_NUMBER) {
        
        x <- str_split(data[a,2], ':', n = Inf, simplify = FALSE)
        Pathway <- x[[1]][2]
        
        pvalue <- data[a,13]
        pvalue <- as.numeric(as.character(pvalue))
        
        GENE_LIST <-  str_split(data[a,6], ', ', n = Inf, simplify = FALSE)[[1]]
        Up = 0
        Down = 0
        for (b in 1:length(GENE_LIST)) {
          GENE_NAME <- GENE_LIST[b]
          row = which(data2[,1] == GENE_NAME)
          
          if (length(row) != 0) {
            if (data2[row,2] >= fc & data2[row,3]<= p_threshold) {
              Up = Up + 1 } else if (data2[row,2] <= -fc & data2[row,3]<= p_threshold) {
                Down = Down + 1}
          }
          

        }
        Number = Up + Down
        output1 <- cbind(Pathway, Up, Down, Number, pvalue)
        if (a ==1) {
          output <- output1
        } else {output <- rbind(output,output1)}
        
      }
      
      write.csv(output,paste0('output_',catogory_list[aa],'.csv'))

      data3= as.data.frame(output)
      data3 <-  data3[data3$pvalue <= p_threshold_path,]
      
      
      
      #colnames(data3) <- c('KEGG pathway','Number of genes',"-Log2(pvalue)")
      
      data3$pvalue <- -log(as.numeric(as.character(data3$pvalue)),10)
      data3 <- data3[order(data3$pvalue),]
      
      pie_x_position <- as.integer(max(data3$pvalue))+2
      data3$pie[1:dim(data3)[1]] <- pie_x_position
      data3$pie <- as.numeric(data3$pie)
      data3$radius[1:dim(data3)[1]] <- pie_radius
      data3$ord[1:dim(data3)[1]] <- (1:dim(data3)[1])
      
      data3_2 = sapply(2:ncol(data3), function(j) (data3[,j] = as.numeric(as.character(data3[,j]))))
      data3[,2:ncol(data3)] = data3_2
      ggplot(data3, aes(x=pvalue, y= reorder(Pathway,pvalue), size = Number))+
        geom_point()+
        labs(x = "-Log10(pvalue)", y = "KEGG pathway", size = "Number of genes")
      
      
      p <- ggplot(data3, aes(x=pvalue, y= reorder(Pathway,pvalue), size = Number))+
        geom_point()+ theme(text = element_text(size = fontsize))+
        labs(x = "-Log10(pvalue)", y = "KEGG pathway", size = "Number of genes")
      p
      
      q <- geom_scatterpie(aes(x=pie, y= ord, r= radius ), data=data3,cols=c('Up','Down'), legend_name = "Gene_expression")
      q
      
      pdf(paste0('enrichment_',catogory_list[aa],'.pdf'))
      print(p+q)
      dev.off()
    }
    
    
  }
}




