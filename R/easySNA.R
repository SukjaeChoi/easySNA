
#' fromtoSNA
#'
#' Analyse SNA for from to data
#' @param docs dataframe. It has from, to, cnt column.
#' @return matrix
#' @export
#' @examples
#' fromtoSNA(docs = docs)

fromtoSNA <- function(docs)
{
  levels_all <- union(levels(docs[,1]), levels(docs[,2]))
  docs$from <- factor(docs[,1], levels=levels_all)   # 발신 목록을 팩터화
  docs$to <- factor(docs[,2], levels=levels_all)     # 수신 목록을 팩터화
  
  docs$from2 <- as.integer(docs$from)
  docs$to2 <- as.integer(docs$to)
  docs$cnt <- as.integer(docs[,3])
  fromto.mat <- cbind(from=docs[,"from2"], to=docs[,"to2"], cnt=docs[,3])
  
  if(!require(igraph)) {install.packages("igraph"); library(igraph)} 
  
  fromto.edge <- graph.edgelist(fromto.mat[, 1:2])   # edge 리스트 추출
  E(fromto.edge)$weight <- fromto.mat[, 3]           # edge마다 횟수 정보로 가중치 부여
  fromto.diag <- rep(0, nrow(docs)) + 5   	         # 자기자신 노드의 기본값(5)
  fromto.name <- levels_all   	                     # 노드의 이름
  
  plot.igraph(fromto.edge, vertex.size=10, vertex.shape="circle", 
              vertex.label=fromto.name, vertex.label.font=1, 
              vertex.label.cex=1+sqrt(fromto.diag)/15, 
              edge.width=2+E(fromto.edge)$weight/2, 
              edge.arrow.width=E(fromto.edge)$weight/50)
  
  return(fromto.edge)
}



#' actionSNA
#'
#' Analyse SNA by action
#' @param docs dataframe
#' @return matrix
#' @export
#' @examples
#' actionSNA(docs = docs)

actionSNA <- function(docs)
{
  docs[is.na(docs)] <- 0
  rownames(docs) <- docs[, 1]  	  # 첫번째 컬럼으로 행의 이름을 붙임
  docs <- docs[-1]             		# 첫번째 컬럼은 필요 없으므로 제외
  docs.mat1 <- as.matrix(docs)  	# 행렬로 변환
  
  n <- nrow(docs.mat1)
  m <- ncol(docs.mat1)
  n.n.mat <- matrix(0, n, n)                       	# n*n의 빈 행렬 생성
  colnames(n.n.mat) <- dimnames(docs.mat1)[[1]]
  
  mat.doc <- cbind(n.n.mat, docs.mat1) ; mat.doc
  doc.mat <- cbind(t(docs.mat1), matrix(0, m, m)) ; doc.mat
  docs.mat2 <- rbind(mat.doc, doc.mat)
  
  if(!require(sna)) {install.packages("sna"); library(sna)} 
  vertex.col <- c(rep("blue", n), rep("green", m))  # 사람이름 수만큼 blue, 책이름 수만큼 green
  vertex.cex <- c(rep(2, n), rep(2, m))             # 사람이름 수와 책이름 수만큼 숫자 2를 생성   	
  
  gplot(docs.mat2, mode = "circle", displaylabels = T, boxed.labels = F, 
        vertex.col = vertex.col, vertex.cex = vertex.cex, 
        label.col = vertex.col, label.cex = 1.2, usearrows = F)
  
  docs.mat3 <- docs.mat1 %*% t(docs.mat1)
  diag(docs.mat3) <- rep(0, dim(docs.mat3)[1])	
  docs.mat3[docs.mat3 > 0] <- 1  		                # 0 이상은 모두 1로 통일
  
  degree(docs.mat3, cmode = "indegree")   	# 입선
  degree(docs.mat3, cmode = "outdegree")  	# 출선
  degrees <- degree(docs.mat3)            	# 입선과 출선을 합한 것
  
  names <- rownames(docs.mat3)
  name.degree <- data.frame(name=names, degree=degrees)
  
  docs.gdist <- geodist(docs.mat3)$gdist
  dimnames(docs.gdist) <- list(as.vector(name.degree$name), as.vector(name.degree$name))
  
  return(docs.gdist)
}
