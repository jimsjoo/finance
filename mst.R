#require( gdata )
require( igraph )
#require( corrplot )

##### a function that standardized the given data
normalize <- function(x) { 
  x <- sweep(x, 2, apply(x, 2, min)) 
  sweep(x, 2, apply(x, 2, max), "/") 
} 

##### a function that calculates Euclidea's distance
distance_euclidian<-function(x){
  return( round( sqrt( 2 * ( 1-x ) ), 2) )
}
  
##### read data
#pl <- "C:/strawberry/perl/bin/perl.exe"
#df <- read.xls ( "ksp200.xlsx", sheet=1, header=TRUE, perl=pl )
df <- read.csv("mst_fut_stk.csv",header=TRUE)

##### set the number of variables
obs<-60

##### normalization 
df.norm <- normalize(df)

##### get pearson correlation matrix
df.corr <- cor( df.norm[,1:obs], method="pearson" )
#df.corr <- cor( df[,1:obs], method="pearson" )

##### get Euclidean distance
df.dist <- distance_euclidian(df.corr)

##### convert the correlation matrix to the dataframe without duplicate values and 1
data <- as.data.frame( as.table( df.dist ) )
combinations <- combn( colnames( df.dist ) , 2 , FUN=function( x ) { paste( x , collapse = "_" ) } )
data <- data[ data$Var1 != data$Var2 , ]
data <- data[ paste( data$Var1 , data$Var2 , sep="_" ) %in% combinations , ]

##### name the columns
names( data ) <- c( "from", "to", "weight" )

#### prepare igraph object, vertices and edges
g<-graph.full(n=obs)
V(g)$name<-colnames(df.corr)

g<-graph.edgelist(cbind(data$from, data$to))
E(g)$weight<-data$weight

##### get Minimum Spanning Trees
mstPrim <- minimum.spanning.tree(g, algorithm = "prim")

#### draw trees
V(g)$name<-colnames(df.corr)
par(mar=c(0,0,0,0)) 
tkplot(mstPrim,layout=layout.kamada.kawai,edge.width=2,
       vertex.size=10,
       vertex.color="cyan",       
       vertex.label=V(g)$name,
       vertex.label.color="blue",
       vertex.label.cex=1.0,
       vertex.label.font=0.0,
       vertex.frame.color="white",
       edge.color="black",
       edge.arrow.size=0.0,
       edge.label=E(mstPrim)$weight,
       edge.label.cex=1.0,
       edge.label.color="red"
)

##### get the list of pairs to be combined with edges
pairs_linked <- cbind(V(g)$name[get.edgelist(mstPrim)[,1]],V(g)$name[get.edgelist(mstPrim)[,2]])
vertex_degree <- cbind(V(g)$name, degree(mstPrim))

for(v in 1:NROW(vertex_degree)){    
   vn <- vertex_degree[v, 1]
   dg <- vertex_degree[v, 2]
   adj <- ''
   irows1 <- which( pairs_linked[,1]==vn )     
   for(i in irows1){          
     adj <- paste( adj, pairs_linked[i,2], sep=" " )
   }
   irows2 <- which( pairs_linked[,2]==vn )     
   for(i in irows2){          
     adj <- paste( adj, pairs_linked[i,1], sep=" " )
   }
   cat(vn, dg, '\t', adj, '\n')
}

#### END OF CODE ####
