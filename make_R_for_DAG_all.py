import sys 

input_path = sys.argv[1] 
input_file = sys.argv[2]
input_table = sys.argv[3]
output_path = sys.argv[4]

output_file = input_file.split('/')[-1].rstrip('.input.tsv')
output=open(output_path + '/' + output_file, 'w')

vec=open(input_file, 'r')
vec_line=''.join(vec.readlines())

w_line = """library(readxl)
library(ggplot2)
library(tidyverse)  # ggplot, dplyr, %>%, and friends
library(ggdag)  # Make DAGs with ggplot
library(dagitty)  # Do basic DAG math
library(broom)
setwd('""" + input_path + """/')
filepath<-'""" + input_path + """/'
filename<-'"""+ input_table + """'
data <- read.table(paste0(filepath,filename),header = TRUE)
data <- data.frame(data)
gg_list = data$geo
gem_list = c()
for(yhname in gg_list){
     yh_g_list = strsplit(yhname,fixed=TRUE ,split='.')[[1]]
     gem_list = unique(c( gem_list , yh_g_list ))}



label <- data$label
names(label) <- data$name
node_dag <- dagify(""" + vec_line + """,coords = data,labels = label)

node_dag_tidy <- node_dag %>% tidy_dagitty()
name_split <- data.frame(t(data.frame(strsplit(node_dag_tidy$data$name, split="_"))))
node_dag_tidy <- node_dag %>% tidy_dagitty() %>%
    mutate(protein=name_split$X1) %>%
    mutate(type=name_split$X2) %>%
    mutate(geo=name_split$X7)

##make node_dag_tidy##
p_color<-dplyr::case_when(
    node_dag_tidy$data$protein == "surface" ~ "#B6C90E",
    node_dag_tidy$data$protein == 'envelope' ~ "#909656",
    node_dag_tidy$data$protein == 'nucleocapsid' ~ "#FCD92B",
    node_dag_tidy$data$protein == 'membrane' ~ "#6A8DFD",
    node_dag_tidy$data$protein == "ORF1ab" ~ "orange",
    node_dag_tidy$data$protein == 'ORF3a' ~ "#998C73",
    node_dag_tidy$data$protein == 'ORF6' ~ "#48CF8B",
    node_dag_tidy$data$protein == 'ORF7a' ~ "#CC9833",
    node_dag_tidy$data$protein == 'ORF7b' ~ "#CC6125",
    node_dag_tidy$data$protein == 'ORF8' ~ "#7BCAD1",
    node_dag_tidy$data$protein =='ORF10' ~ "#D0FF94",
    TRUE ~ "black"
)
p_color<-data.frame(p_color)


node_dag_tidy <-node_dag_tidy %>% mutate(p_color)
#c("Asia", "Africa", "Europe", "Oceania", "NorthAmerica", "SouthAmerica", "NA")

###make plot Automatically
geom_list<-c("Asia", "Africa", "Europe", "Oceania", "NorthAmerica", "SouthAmerica") ##, "NA")

for (country in gem_list){
    print(country)
    p_color<-dplyr::case_when(
        node_dag_tidy$data$protein == "surface" ~ "#B6C90E",
        node_dag_tidy$data$protein == 'envelope' ~ "#909656",
        node_dag_tidy$data$protein == 'nucleocapsid' ~ "#FCD92B",
        node_dag_tidy$data$protein == 'membrane' ~ "#6A8DFD",
        node_dag_tidy$data$protein == "ORF1ab" ~ "orange",
        node_dag_tidy$data$protein == 'ORF3a' ~ "#998C73",
        node_dag_tidy$data$protein == 'ORF6' ~ "#48CF8B",
        node_dag_tidy$data$protein == 'ORF7a' ~ "#CC9833",
        node_dag_tidy$data$protein == 'ORF7b' ~ "#CC6125",
        node_dag_tidy$data$protein == 'ORF8' ~ "#7BCAD1",
        node_dag_tidy$data$protein =='ORF10' ~ "#D0FF94",
        TRUE ~ "black"
    )
    p_color<-data.frame(p_color)

    node_dag_tidy <- node_dag %>% tidy_dagitty()
    name_split <- data.frame(t(data.frame(strsplit(node_dag_tidy$data$name, split="_"))))
    node_dag_tidy <- node_dag %>% tidy_dagitty() %>%
        mutate(protein=name_split$X1) %>%
        mutate(type=name_split$X2) %>%
        mutate(geo=name_split$X7) %>%
        mutate(p_color)

    detail<-ifelse(grepl(country,node_dag_tidy$data$geo)==TRUE,"YES","YES")
    detail<-data.frame(detail)
    node_dag_tidy <- node_dag_tidy %>% mutate(detail)

    TF<-grepl(country,node_dag_tidy$data$geo)
  ##  for ( i in c(1:length(TF))){
  ##      if(TF[i] == FALSE){
  ##          node_dag_tidy$data$p_color[i]<-"gray"
	    
  ##      }
  ##  } #
    print(country)
    print(node_dag_tidy)
    name=paste0(country,":","Mutation Type Change")
    ##make ggplot##
    ggplot(node_dag_tidy, aes(x = node_dag_tidy$data$x, y = node_dag_tidy$data$y, xend = xend, yend = yend)) +

        geom_point(aes(color=node_dag_tidy$data$p_color),color=node_dag_tidy$data$p_color,shape=ifelse(node_dag_tidy$data$protein=="Merge",15,ifelse(node_dag_tidy$data$type=="nonsyn",21,16)),size=0.5,alpha=0.7,stroke=2) +
        coord_cartesian(xlim = c(0,19)) +
        geom_text(aes(label=node_dag_tidy$data$label),size = 2, color="#3E4032",vjust = 3, hjust=0.5) +
        geom_dag_edges_arc(curvature = 0.000001,arrow = grid::arrow(length = grid::unit(2, "pt"), type = "closed"), na.rm = TRUE,position = "identity",n = 50,lineend = "square",linejoin = "round", edge_width=0.05,edge_color="#EBE4DF") +


        labs(x="Month",y="",title = (name)) +
        theme(legend.position = "bottom",legend.key.size = unit(1, 'mm'),panel.background = element_rect(fill = 'white'),legend.title=element_blank(),panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"),plot.title = element_text(hjust = 0.5),axis.text.y = element_blank()) +
        scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15, 16,17,18,19),labels = paste0(c(201912,202001,202002,202003,202004,202005,202006,202007,202008,202009,202010,202011,202012,202101,202102,202103, 202104,202105,202106,202107)))

        ##ggplot save##
        output=gsub(".table.tsv",".pdf",filename)
        ggsave(paste0(filepath,country,"_",output),width = 30,height = 20 ,units=c("cm"), limitsize = FALSE)
         }"""



output.write(w_line)


output.close()
vec.close()
