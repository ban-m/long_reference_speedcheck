library("tidyverse")
setwd("~/work/long_reference_speedcheck/")

rawdata <- read_csv("result/speedcheck.csv")

sub <- rawdata %>% filter(bandwidth==0)
generalplot <- function(g,name){
    pdf(paste0("./pdf/",name,".pdf"),width=21,height=21)
    plot(g)
    dev.off()
    png(paste0("./png/",name,".png"),width = 480*2,height = 480*2)
    plot(g)
    dev.off()
}    


g <- rawdata %>%
    filter(bandwidth != 0) %>%
    mutate(bandwidth = as.factor(bandwidth)) %>% 
    ggplot(mapping=aes(x = refsize,y = time,colour = bandwidth)) + geom_line()+
    geom_line(data = sub,mapping = aes(x =refsize,y = time),colour = "red",size = 2) +
    geom_hline(yintercept = 0.3,colour = "blue",size = 2)

g <- sub%>% ggplot(mapping = aes(x = refsize,y = time)) + geom_line() +
    geom_hline(yintercept = 0.3,colour = "blue",size = 1)+
    labs(x = "reference length(bp)",y="elapsed time(sec)") +
    theme_bw(base_size = 28)

generalplot(g,"elapsed_time")
