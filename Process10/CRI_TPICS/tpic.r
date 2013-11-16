Arguments <- commandArgs(trailingOnly = T)
L = as.numeric(Arguments[1])# Average length of DNA fragments
print(paste("L is ",L))
Genome = Arguments[2]               ##
sgr.file<- "Coverage"
bed.file  <- "TPIC"
zeta.file <-"ZetaRegions.txt"
Sh.param<-7 # 
tree.p.value<-0.01 
num.sims<-30000
min.site.length<-36 #minimum length (in bp) of a peak ... Maybe read length would be nice
if(Genome %in% "hg18" || Genome %in% "GRCh37"){
	index.set<-c(1:22, "X","Y") #set of chromosomes
}
if(Genome %in% "mm8" || Genome %in% "mm9"){
	index.set<-c(1:24, "X","Y") #set of chromosomes
}

#########################
findSh<-function(theta, H,term){
	sum<-1
	p<-1-1/theta +exp(-theta)/theta
	for(k in 2:(H+1)){
		j<-k-1
		p<-1-(k/theta)*p
	}
	sum<-sum+ 1 / (1 - p)
	new.term<-1/(1-p)
	i<-H+2
	if(term > new.term){
		term<-new.term
		while(log(new.term)>-300){
			new.term<-new.term*theta/i
			sum<-sum+new.term
			i<-i+1
		}   
		return(c(sum,term))
	}else{
		return(0)
	}
}
#########################
randomTree<-function(theta, H=0){
prev.state<-H
curr.state<-H+1
lattice.nodes<-c(prev.state,curr.state)
lattice.mins<-c(prev.state)
p<-list()
p[[1]]<-1-1/theta +exp(-theta)/theta
if(H>0){
for(k in 2:(H+1)){
	j<-k-1
	p[[k]]<-1-(k/theta)*p[[j]]
}}
while(curr.state>H){
	if(length(p)<curr.state){
		t<-curr.state-1
		p[[curr.state]]<-1-(curr.state/theta)*p[[t]]
	}
	prob.up<-p[[curr.state]]
	if(prob.up<0){
		prob.up<-0
	}else if(prob.up>1){
		prob.up<-1
	}
	prob.down<-1-prob.up
	go.up<-curr.state+1
	go.down<-curr.state-1
	next.state<-sample(c(go.up,go.down),1, prob=c(prob.up,prob.down))
	if(next.state==curr.state+1){
		lattice.nodes<-c(lattice.nodes, next.state)
	}
	if(curr.state<next.state && curr.state<prev.state){
		lattice.mins<-c(lattice.mins, curr.state)
	}else if(curr.state+1==next.state && curr.state==prev.state+1){
		lattice.mins<-c(lattice.mins, curr.state)
	}
	prev.state<-curr.state
	curr.state<-next.state
}
lattice.mins<-c(lattice.mins, lattice.nodes[length(lattice.nodes)])
lattice.path<-cbind(lattice.nodes,lattice.mins)
Max.m<-MaximalM(lattice.path,0)
return(Max.m)
}

#####################################
simulateT<-function(theta, H,num){
matchings<-matrix(nrow=1, ncol=num)
	for(t in 1:num.sims){
		match.sim<-randomTree(theta,H)
		matchings[t]<-length(match.sim)
	}

	sim.list<-sort(matchings)


return(sim.list)
}

#########################
IdentifyPeaks<-function(){
#system("perl zeta.pl")
smaller.regions<-read.table(file=zeta.file)
OuTF <- file(paste(bed.file, "bed", sep="."), "w")
counts.tot<-as.matrix(smaller.regions[,2:5])
input.thetas<-unique(counts.tot[,4])
data.thetas<-c()
for(j in input.thetas){
	intvls<-which(counts.tot[,4]==j)
	intvl.length<-sum(as.numeric(counts.tot[intvls,2])-as.numeric(counts.tot[intvls,1]))
	intvl.theta<-sum(counts.tot[intvls,3])*L/intvl.length
	data.thetas<-c(data.thetas,intvl.theta)
}
rm(counts.tot)
sim.list<-list()
H<-200
Sh<-findSh(200,H,1000)
while(Sh[1]> Sh.param){
	H<-H+1
	Sh<-findSh(200,H, Sh[2])
}
Tall.Mchs<-simulateT(200,H,num.sims)
Zeros<-which(data.thetas==0)
data.thetas[Zeros]<-.01	
for(i in 1:length(data.thetas)){
	H<-max(ceiling(data.thetas[i]),2)
	Sh<-findSh(data.thetas[i],H,1000)
	while(Sh[1]> Sh.param){
		H<-H+1
		Sh<-findSh(data.thetas[i],H,Sh[2])
	}
	if(data.thetas[i]<200){
		Mchs<-simulateT(data.thetas[i],H,num.sims)
	}else{
		Mchs<-Tall.Mchs
	}
	sim.list[[input.thetas[i]]]<-c((H+1),Mchs)
}
	print("Simulation of random trees is finished!", quote=FALSE)

for(chromo in index.set){
	chrom<-paste("chr",chromo,sep="")
	regions.chr<-smaller.regions[which(smaller.regions[,1]==chrom),]
	regions.chr<-as.matrix(regions.chr[,c(2,3,5)])

	data<-read.table(file=paste(sgr.file, chrom, "sgr", sep="."))
	data<-as.matrix(data[,2:3])

	for(k in 1:nrow(regions.chr)){
		end.i<-which(data[,1]<=regions.chr[k,2])
		if(length(end.i)>0){
			end.i<-end.i[length(end.i)]
			input.theta<-regions.chr[k,3]
			###
			H.mil<-sim.list[[input.theta]][1]
			in.order<-sim.list[[input.theta]][-1]	
			new.data<-data[1:end.i,,drop=FALSE]
			##############
			indices<-which(new.data[,2]>=H.mil)
			shape<-matrix(ncol=5)
			start<-1
			gap<-1
			end<-start+gap
			while(end<=length(indices)){
				if(indices[end]==indices[start]+gap){
					gap<-gap+1
					end<-start+gap
				}else{
					pot.peak<-new.data[indices[start]:indices[(end-1)],, drop=FALSE]
					last.end<-new.data[(indices[(end-1)]+1),1]-1
					crit<-criticalPoint(pot.peak,last.end)
					if(length(crit)>6){
						crit<-crit[,c(1,4,5)]
					}
					if(max(pot.peak[,2])-H.mil > 2*max(in.order)){
						tree.stat<-(max(pot.peak[,2])-H.mil)/2
						p.value<-0
					}else{
						tree.stat<-MaximalM(crit, H.mil)
						tree.stat<-length(tree.stat)
						count.sim<-which(in.order>=tree.stat)
						p.value<-length(count.sim)/num.sims
					}
					if(length(crit)==3){
						if(is.na(shape[1,1])){
							shape[1,]<-c(crit[1], crit[3], crit[2],tree.stat,p.value)
						}else{ shape<-rbind(shape, c(crit[1], crit[3], crit[2],tree.stat,p.value))}
					}else{
						if(is.na(shape[1,1])){
							shape[1,]<-c(crit[1,1], crit[nrow(crit),3], max(crit[,2]),tree.stat,p.value)
						}else{ shape<-rbind(shape, c(crit[1,1], crit[nrow(crit),3], max(crit[,2]),tree.stat,p.value))}
					}
					start<-end
					gap<-1
					end<-start+gap
				}
			}
			if(length(indices)-start>1){	
				pot.peak<-new.data[indices[start]:indices[length(indices)],, drop=FALSE]
				if(indices[length(indices)]==end.i){
					last.end<-data[(end.i+1),1]-1
				}else{
					last.end<-indices[length(indices)]+1
					last.end<-new.data[last.end,1]-1
				}
				crit<-criticalPoint(pot.peak,last.end)
				if(length(crit)>6){crit<-crit[,c(1,4,5)]}
				if(max(pot.peak[,2])-H.mil > 2*max(in.order)){
					tree.stat<-(max(pot.peak[,2])-H.mil)/2
					p.value<-0
				}else{
					tree.stat<-MaximalM(crit, H.mil)
					tree.stat<-length(tree.stat)
					count.sim<-which(in.order>=tree.stat)
					p.value<-length(count.sim)/num.sims
				}
				if(length(crit)==3){
					if(is.na(shape[1,1])){
						shape[1,]<-c(crit[1], crit[3], crit[2],tree.stat,p.value)
					}else{ shape<-rbind(shape, c(crit[1], crit[3], crit[2],tree.stat,p.value))}
				}else{
					if(is.na(shape[1,1])){
						shape[1,]<-c(crit[1,1], crit[nrow(crit),3], max(crit[,2]),tree.stat,p.value)
					}else{ shape<-rbind(shape, c(crit[1,1], crit[nrow(crit),3], max(crit[,2]),tree.stat,p.value))}
				}
			}
			data<-data[-c(1:end.i),, drop=FALSE]
			if(is.na(shape[1,1])){	
			}else{
				for(p in 1:nrow(shape)){
					if(shape[p,2]-shape[p,1]>=min.site.length && p==1){
						new.tag<-paste(chrom, formatC(shape[p,1], format="f", digits=0),formatC(shape[p,2], format="f", digits=0),shape[p,3], shape[p,4], shape[p,5], 1 ,"\n", sep=" ")
						cat(new.tag, file = OuTF)
					}else if(shape[p,2]-shape[p,1]>=min.site.length){
						new.tag<-paste(chrom, formatC(shape[p,1], format="f", digits=0),formatC(shape[p,2], format="f", digits=0),shape[p,3], shape[p,4], shape[p,5], 0 ,"\n", sep=" ")
						cat(new.tag, file = OuTF)
					}
				}
			}
		}
	}			
print(paste("Trees from chromosome", chrom, "have been created and scored.", sep=" "),quote=FALSE)
}
close(OuTF)
unlink(OuTF)

possible.k<-c()
found.trees<-read.table(file=paste(bed.file, "bed", sep="."))
found.probs<-as.matrix(found.trees[,6])
found.probs<-sort(found.probs)
ratio.p<-tree.p.value/length(found.probs)
for(k in 1:length(found.probs)){
	if(found.probs[k]<=k*ratio.p){
		possible.k<-c(possible.k,k)
	}
}
sig.p.value<-max(possible.k)*ratio.p
sig.found.p<-which(as.matrix(found.trees[,6])<=sig.p.value)
found.trees<-found.trees[sig.found.p, , drop=FALSE]
OuTMerged<-file(paste(bed.file, "final.bed",sep="."), "w")

for(chromo in index.set){
        chrom<-paste("chr",chromo,sep="")
	starts.chrom<-as.numeric(found.trees[which(found.trees[,1]==chrom), 2])
	ends.chrom<-as.numeric(found.trees[which(found.trees[,1]==chrom), 3])
	new.region<-as.numeric(found.trees[which(found.trees[,1]==chrom), 7])
	if(length(starts.chrom)>0){
		start.chr<-starts.chrom[1]
		end.chr<-ends.chrom[1]
		if(length(starts.chrom)>1){
			for(j in 2:length(starts.chrom)){
				new.region.ind<-new.region[j]
				new.end<-ends.chrom[j]
				new.start<-starts.chrom[j]
				if(new.region[j]==1 && new.start - end.chr < L){
					end.chr<-new.end
				}else{
					next.tag<-paste(chrom, formatC(start.chr, format="f", digits=0), formatC(end.chr, format="f", digits=0), "\n", sep=" ")
					cat(next.tag, file = OuTMerged)
					start.chr<-new.start
					end.chr<-new.end
				}
			}
		}
		next.tag<-paste(chrom, formatC(start.chr, format="f", digits=0), formatC(end.chr, format="f", digits=0), "\n", sep=" ")
		cat(next.tag, file = OuTMerged)
	}
}
close(OuTMerged)
unlink(OuTMerged)
print(paste("Trees with p-value less than", sig.p.value, "are statistically significant",sep=" "),quote=FALSE)
}

#########################
# Function returns a n x 5 matrix, where n is the number
# of critical points.  The first column is the nucleotide position
# (first column of data), fourth column is the count 
# (second column of data), and the fifth column gives the ending position for the node 
#--all nucleotides between 1st col and 5th col have same value (in 4th col).  The two middle columns are zero.
criticalPoint<-function(data,last.end=0)
{

	F1<-data[1,]
	end.vect<-c()
if(length(data)<=2){
  return(c(data,last.end))
}
# Removes consecutive numbers that are equal.
	for(i in 2:nrow(data)){
		k<-i-1;
		if(data[k,2] != data[i,2] ){
			F1<-rbind(F1,data[i,])
			end.pos<-data[i,1]-1
			end.vect<-c(end.vect,end.pos) 
		}
		
	}
	end.vect<-c(end.vect, last.end)
	if(length(F1)>2){
		F1<-cbind(F1,end.vect)
	}else{ F1<-c(F1,end.vect)}
	if(length(F1)<=6){
		return(F1)
	}else{
		num2<-nrow(F1)-1
		F2<-c(F1[1,1],0,0, F1[1,2], F1[1,3]) 
		for(i in 2:num2){
			j<-i+1
			k<-i-1
			if (F1[k,2] > F1[i,2] && F1[i,2] < F1[j,2]){
				v<-c(F1[i,1], 0, 0, F1[i,2], F1[i,3])
				F2<-rbind(F2,v)
			}
			else if (F1[k,2] < F1[i,2] && F1[i,2] > F1[j,2]){ 
				v<-c(F1[i,1], 0, 0, F1[i,2], F1[i,3])
				F2<-rbind(F2,v)
			}
		}
		i<-num2+1
		v<-c(F1[i,1], 0, 0, F1[i,2], F1[i,3])
		F2<-rbind(F2,v)
		return(F2)		
	}
}
#########################
# This function takes in data, a 2 column (nucl,  height) matrix and returns
# a newick string for the rooted binary tree associated to the data.
getNewick<-function(data){
	S<-c()
	Fe<-criticalPoint(data)
	if(length(Fe)>10){
		F2<-Fe[,1:4]
	}else{
		return()
	}
	M<-max(F2[,4])  # M is maximum of "count" values for nucleotides
	i=M;
# while loop below starts with highest "count" values as these 
# excursion sets appear first
	while(i>=0){
		j=2
		while(j<nrow(F2)){ 
			u<-j-1
			p<-j+1
			if((F2[j,4]==i) && (F2[u,4]>=i) && (F2[p,4]>i)){ # will merge two nodes at height i
				dist1<-F2[u,4]-F2[j,4]
				dist2<-F2[p,4]-F2[j,4]
				v<-c(F2[j,1], F2[u,1], F2[p,1], F2[j,4]) # node indexed by j-1 is merging with that indexed by j+1 at spot j
				if(F2[u,2]==0 && F2[p,2]==0){ # tests to see if node indexed by j-1 and j+1 are both leaves
					str1<-paste('(',F2[u,1],':', dist1,',', F2[p,1],':', dist2,')', sep="")
					str2<-paste(F2[j,1],':', F2[j,4], sep="")
					newstr<-paste(str1,str2, sep="");
					S<-c(S,newstr)
				}else if(F2[u,2]>0){ # tests to see if node j-1 is NOT a leaf 
# (i.e. nucleotide F2(j-1,2) was an earlier node merged here with nucleotide F2(j-1,3))
					test1<-paste(")",F2[u,1], ":",sep="")
					test2<-paste(")",F2[p,1], ":",sep="")        
					f<-grep(test1,S)  
					chop.newick<-unlist(strsplit(S[f], as.character(F2[u,1])))
					s<-paste("(",chop.newick[1],as.character(F2[u,1]), sep="")   # chops off the end of the Newick string
# constructed previously, since node j-1 is 
# being merged with another (i.e. it is not adjacent 
# to the root
					S<-S[-f]	  # removes old Newick string
					
					if(F2[p,2]>0){ # tests to see if additionally node j+1 is NOT a leaf
						seh=grep(test2,S)	
						chop.newick<-unlist(strsplit(S[seh], as.character(F2[p,1])))
						s2<-paste(chop.newick[1],as.character(F2[p,1]), sep="")  
						S<-S[-seh]
						
					}
					else{ s2<-as.character(F2[p,1])
					}
					P<-paste(s, ':', dist1,',', s2,':', dist2,')' , sep="")
					P2<-paste(P, F2[j,1],':', F2[j,4] , sep="")
					S<-c(S,P2)
				}
				else{ # node j-1 is a leaf but node j+1 is not a leaf
					test2<-paste(")",F2[p,1], ":",sep="")
					f<-grep(test2,S)
					chop.newick<-unlist(strsplit(S[f], as.character(F2[p,1])))
					s<-paste(chop.newick[1],as.character(F2[p,1]) , sep="")  
					S<-S[-f]
					P<-paste("(", F2[u,1], ":", dist1,",", s, ":", dist2,")", sep="")
					P2<-paste(P, F2[j,1],":", F2[j,4], sep="")
					S<-c(S,P2)
					
				}
				F2[u,]<-v
				F2<-F2[-c(j,p), ,drop=FALSE]

			} 
			else{ j<-j+1
			}
		}
		i<-i-1
		
	}
	crit=paste(S[length(S)], ';', sep="");
	return(crit)
}

##########################
MaximalM<-function(crit, simulated){
if(simulated>0){
if(length(crit)>3){ 
	if(crit[1,2]<crit[2,2]){
		crit<-crit[-1,,drop=FALSE]
	}}
if(length(crit)>3){ 
	t<-nrow(crit)
	if(crit[t,2]<crit[(t-1),2]){
		crit<-crit[-t,,drop=FALSE]
	}}
if(length(crit)<=3){
	num.verts<-crit[2]-simulated +1
	matching.size<-floor(num.verts/2)
	return(c(1:matching.size))
}else{
new.crit<-crit[,2]
lowest.v<-c(simulated:(new.crit[1]-1), new.crit[2])
graph.mat<-cbind(simulated:new.crit[1], lowest.v)
p<-2
while(p<length(new.crit)-1){
	if(new.crit[p+1]-new.crit[p] >1){
		lowest.v<-c((new.crit[p]+1):(new.crit[p+1]-1), new.crit[p+2])
		new.addition<-cbind((new.crit[p]+1):new.crit[p+1], lowest.v)
	}else if(new.crit[p+1]==new.crit[p]+1){
		new.addition<-c(new.crit[p+1], new.crit[p+2])
	}
	graph.mat<-rbind(graph.mat,new.addition)
	p<-p+2
}
if(p<-length(new.crit)-1){
	new.addition<-cbind((new.crit[p]+1):new.crit[p+1], (new.crit[p]+1):new.crit[p+1])
	graph.mat<-rbind(graph.mat,new.addition)
}
}
}else if(simulated==0){
	graph.mat<-crit
}
Matching<-list()
edg.vect<-c()
for(j in 2:nrow(graph.mat)){
	sub.index<-which(graph.mat[,1]==(graph.mat[j,1]-1))
	k<-j-1
	while(k>=1&& graph.mat[k,2]>=(graph.mat[j,1]-1)){
		if(length(which(sub.index==k))>0){
			edg.vect<-rbind(edg.vect, c(k,j))
			k<-k-1
		}else{ k<-k-1}
	}
}
if(is.null(edg.vect)){

}else{
while(nrow(edg.vect)>0){
	curr.vert<-1
	while(curr.vert<=nrow(edg.vect)){
		left.adj.vert<-which(edg.vect[,1]==edg.vect[curr.vert,1])
		right.adj.vert<-which(edg.vect[,2]==edg.vect[curr.vert,1])
	if(length(left.adj.vert)+length(right.adj.vert)==1){
		t<-length(Matching)+1
		Matching[[t]]<-edg.vect[curr.vert,]
		adjacent.vert<-edg.vect[curr.vert,2]
		edg.vect<-edg.vect[-curr.vert,, drop=FALSE]
		left.adj.vert<-which(edg.vect[,1]==adjacent.vert)
		if(length(left.adj.vert)>0){
			edg.vect<-edg.vect[-left.adj.vert,, drop=FALSE]
		}
		right.adj.vert<-which(edg.vect[,2]==adjacent.vert)
		if(length(right.adj.vert)>0){
			edg.vect<-edg.vect[-right.adj.vert,, drop=FALSE]
		}		
	}else{
		left.adj.vert<-which(edg.vect[,1]==edg.vect[curr.vert,2])
		right.adj.vert<-which(edg.vect[,2]==edg.vect[curr.vert,2])
		if(length(left.adj.vert)+length(right.adj.vert)==1){
			t<-length(Matching)+1
			Matching[[t]]<-edg.vect[curr.vert,]
			adjacent.vert<-edg.vect[curr.vert,1]
			edg.vect<-edg.vect[-curr.vert,, drop=FALSE]
			left.adj.vert<-which(edg.vect[,1]==adjacent.vert)
			if(length(left.adj.vert)>0){
				edg.vect<-edg.vect[-left.adj.vert,, drop=FALSE]
			}
			right.adj.vert<-which(edg.vect[,2]==adjacent.vert)
			if(length(right.adj.vert)>0){
				edg.vect<-edg.vect[-right.adj.vert,, drop=FALSE]
			}	
		}else{curr.vert<-curr.vert+1} 
	
	}
	}
}}	
return(Matching)
}

IdentifyPeaks()
