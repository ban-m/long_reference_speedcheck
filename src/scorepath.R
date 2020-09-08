library("tidyverse")
setwd("~/work/long_reference_speedcheck/")
rawdata <- read_csv("./result/scorepath.csv")

generalplot <- function(g,name){
    pdf(paste0("./pdf/",name,".pdf"))
    plot(g)
    dev.off()
    png(paste0("./png/",name,".png"))
    plot(g)
    dev.off()
}    

plot_each_bandwidth <- function(id,opt_score,opt_position,bandwidth,df){
    name <- paste0(id,"band-",formatC(bandwidth,width = 2,flag="0"))
    wide <- df %>%filter(bandwidth !=0)
    g <- wide %>%
        ggplot(mapping = aes(x = position,y = score)) + 
        geom_line(show.legend = FALSE) +
        geom_line(mapping=aes(x = position,y = kl_divergence),colour = "green",alpha=0.3) + 
        geom_vline(xintercept = opt_position,colour = "red",size = 1,alpha = 0.2) +
        geom_hline(yintercept = opt_score,colour = "red",size = 1,alpha = 0.2) 
    generalplot(g,name)
    narrow <- df %>% filter(bandwidth !=0) %>%
        filter(opt_position -100 < position,position <opt_position+100)
    g2 <- narrow %>%
        ggplot(mapping = aes(x = position,y = score)) + 
        geom_line(mapping=aes(x = position,y = kl_divergence),colour = "green",alpha=0.3) + 
        geom_line(show.legend = FALSE) + 
        geom_vline(xintercept = opt_position,colour = "red",size = 1,alpha = 0.2) +
        geom_hline(yintercept = opt_score,colour = "red",size = 1,alpha = 0.2) 
    generalplot(g2,paste0(name,"neighbor"))
}

plot_path <- function(id,df){
    print(id)
    name <- paste0(id,"all")
    optimal <- df %>% filter(bandwidth == 0)
    opt_score <- optimal %>% pull(score)
    opt_position <- optimal %>% pull(position)
    g <- df %>%filter(bandwidth !=0) %>%
        ggplot(mapping = aes(x = position,y = score,colour = bandwidth)) + 
        geom_line(show.legend = FALSE) + 
        geom_vline(xintercept = opt_position,colour = "red",size = 3,alpha = 0.2) +
        geom_hline(yintercept = opt_score,colour = "red",size = 3,alpha = 0.2)
    generalplot(g,name)
    g2 <- df %>% filter(bandwidth !=0) %>%
        filter(opt_position -100 < position,position <opt_position+100) %>%
        ggplot(mapping = aes(x = position,y = score,colour = bandwidth)) + 
        geom_line(show.legend = FALSE) + 
        geom_vline(xintercept = opt_position,colour = "red",size = 3,alpha = 0.2) +
        geom_hline(yintercept = opt_score,colour = "red",size = 3,alpha = 0.2)
    generalplot(g2,paste0(name,"neighbor"))
    df %>% filter(bandwidth != 0) %>% nest(-bandwidth) %>%
        apply(MARGIN=1,
              FUN = function(ls)plot_each_bandwidth(id,
                                                    opt_score,
                                                    opt_position,
                                                    ls$bandwidth,
                                                    ls$data))
    0
}

rawdata %>% mutate(kl_divergence = 1000*kl_divergence) %>% nest(-id) %>% apply(MARGIN = 1,FUN = function(ls){plot_path(ls$id,ls$data)})

sample1 <- rawdata %>% nest(-id) %>% slice(1) %>% unnest() %>% 
    mutate(kl_divergence = 1000*kl_divergence)  %>% select(-id)

sample2 <- rawdata %>% nest(-id) %>% slice(2) %>% unnest() %>% 
    mutate(kl_divergence = 1000*kl_divergence) %>% select(-id)
    
optimalloc<- list(one=sample1 %>% filter(bandwidth == 0) %>% pull(position),
                  two=sample2 %>% filter(bandwidth == 0) %>% pull(position))

