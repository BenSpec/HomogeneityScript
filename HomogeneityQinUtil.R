#Basic Lexicon: some, all, no, notAll, the, notThe, with the the same in the and notThe
#worlds =  onlySome, all, no
#FlatPriors on Lexica

 Some<-c(0,1,1)
 All <-c(0,0,1)
 No<-c(1,0,0)
 NotAll<-c(1,1,0)
 nworlds<-3
cost<-c(1,1,1.5,2,0,1)
#cost<-c(,5,6,7,0,2)
 ProbaWorlds<-c(1/3, 1/3, 1/3)
 temp<-5
 iter<-5

#Refined lexica

LexSome<-list(Some = Some)
LexAll<-list(All = All)
LexNo<-list(No = No)
LexNotAll<-list(NotAll = NotAll)
LexThe<-list(c(0,1,1), c(0,0,1))
Lex<-list(Some = LexSome, All = LexAll, No = LexNo, NotAll = LexNotAll, The = LexThe)
Lexica<-expand.grid(Lex)
Lex1<-unlist(Lexica[1,])
Lex11<-matrix(Lex1, nrow=3, ncol=5)
Lex11<-cbind(Lex11, c(1,0,0))
Lex2<-unlist(Lexica[2,])
Lex22<-matrix(Lex2, nrow=3, ncol=5)
Lex22<-cbind(Lex22, c(1,1,0))
LEX<-list(Lex11, Lex22)
for (i in 1:length(LEX))
   {colnames(LEX[[i]])<-c("some","all", "no", "NotAll", "the", "NotThe")
   	row.names(LEX[[i]])<-c("Zero","JustSome","All")}
ProbaLEX<-c(rep.int(1/length(LEX), length(LEX)))


Worlds<-c(1,2,3)
Qsome<-function(x)
{if (x==1)
	return (c(1))
 if ((x==2) || (x==3))
	return (c(2,3))
}

Qall<-function(x)
{if (x==3)
	return (c(3))
if ((x==1)||(x==2))
	return (c(1,2))
 
}

Qtotal<-function(x)
{return (c(x))}

Qs<-list(Qsome = Qsome, Qall= Qall, Qtotal = Qtotal)
ProbaQ<-c(1/3, 1/3, 1/3)


#Literal listener
L0<-list()
for (i in 1:length(LEX))
{L_0aux<-LEX[[i]]*ProbaWorlds
L_0<-t(t(L_0aux)/rowSums(t(L_0aux)))
L0<-append(L0, list(L_0))
}

#utility of a message in a world, for each Q and Lexicon
UsuperList<-list()
for (l in (1:length(LEX)))
{
{UList<-list()
for (q in (1:length(Qs)))	
{
listL_Q<-list()
for (w in (1:length(Worlds)))
{X<-rbind(L0[[l]][Qs[[q]](w),], rep(c(0),length(cost)))
listL_Q<-append(listL_Q, list(colSums(X)))}
L_Q<-do.call(rbind, listL_Q)
U<-t(log(t(L_Q)) - cost)
UList<-append(UList, list(U))}}
UsuperList<-append(UsuperList, list(UList))}

WeighedU<-function(l,q)
{return (ProbaLEX[l]*UsuperList[[l]][[q]])}

WeighedUlistbyQ<-function(q)
{WeighedUlist<-list()
	for (l in 1:length(LEX))
	{X<-WeighedU(l,q)
	 WeighedUlist<-append(WeighedUlist, list(X))}
	return(WeighedUlist)
	
	
	}

U1_q<-function(q)
{return(Reduce('+',WeighedUlistbyQ(q)))}

U1<-lapply(1:length(Qs), U1_q)


SoftMax<-function(x)
{return(exp(temp*x))}

S_aux <- lapply(U1, SoftMax)


NormS <-function(x)
{return(x/rowSums(x))}


S1<-lapply(S_aux, NormS)


#Listener
#L(w,Q|u) \propto P(w)*P(Q)*S(u|w,Q)
#L(w|u) = Sum_Q(P(Q)L(w,Q|u))
#U(u|Q,w) = log(L(Q(w), Q|u)) - c(u)
S<-S1
for (k in 1:iter)
{
L_aux<-list()
for (u in 1:length(cost))
{X<-matrix(,nrow=nworlds, ncol = length(Qs))
  for (w in 1:nworlds)	
	{for (q in 1:length(Qs))
	
		{X[w,q]<-ProbaWorlds[w]*ProbaQ[q]*S[[q]][w,u]}}
L_aux<-append(L_aux, list(X))	
}



L_u<-lapply(L_aux, function(x){x/sum(x)})

L_fun<-function(u,q,w)
{return(L_u[[u]][w,q])}


L<-list()
for (q in 1:length(Qs))
{M<-outer(1:nworlds, 1:length(cost), FUN = Vectorize(function(w,u){L_fun(u,q,w)}))
L<-append(L, list(M))	
	
	 	}

L_qWq_u<-function(u,q,w)
{return(sum(L[[q]][Qs[[q]](w),u]))}

# L_qWq_u<-L
# for (q in 1:length(Qs))
# {for (w in 1:nworlds)
# {for (u in 1:length(cost))
 # L_qWq_u[[q]][w,u]<-sum(L[[q]][Qs[[q]](w),u])	
	
	   # }
# }	

U_fun<-function(u,q,w)
{return(log(L_qWq_u(u,q,w))-cost[u])}

S_fun_aux<-function(u,q,w)
{return(SoftMax(U_fun(u,q,w)))}


S_aux<-list()
for (q in 1:length(Qs))
{
 M<-outer(1:nworlds, 1:length(cost), FUN = Vectorize(function(w,u){S_fun_aux(u,q,w)})) 	
S_aux<-append(S_aux, list(M))	}


S<-lapply(S_aux, NormS)

}


Lbis<-list()
for (u in 1:length(cost))
{M<-outer(1:nworlds, 1:length(Qs), FUN = Vectorize(function(w,q){L_fun(u,q,w)}))
Lbis<-append(Lbis, list(M))}

names(Lbis)<-c("some","all", "no", "NotAll", "the", "NotThe")

for (i in 1:length(Lbis))
   {colnames(Lbis[[i]])<-names(Qs)
   	row.names(Lbis[[i]])<-c("Zero","JustSome","All")}


L_q<-sapply(Lbis, FUN = function(x){colSums(x)})
L_w<-sapply(Lbis, FUN = function(x){rowSums(x)})



names(S)<-names(Qs)

for (i in 1:length(S))
   {colnames(S[[i]])<-c("some","all", "no", "NotAll", "the", "NotThe")
   	row.names(S[[i]])<-c("Zero","JustSome","All")}

