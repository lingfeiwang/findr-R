# Copyright 2016, 2017 Lingfei Wang
# 
# This file is part of Findr.
# 
# Findr is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Findr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with Findr.  If not, see <http://www.gnu.org/licenses/>.
# 
findr.gtype='integer'
findr.ftype='double'
as.gtype=as.integer
as.ftype=as.double

findr.libload=function(loglv,rs,nth)
{
	.C("external_R_lib_init",loglv,rs,nth)
}

findr.lib=function(loglv=6,rs=NULL,nth=NULL){
	if((length(loglv)>1)||(!is.null(rs)&&(length(rs)>1))||(!is.null(nth)&&(length(nth)>1)))
		stop('All parameters must be singlet.')
	loglvi=as.integer(loglv)
	if(is.null(rs))
	{rsi=as.integer(0)}
	else
		rsi=as.integer(rs)
	if(is.null(nth))
	{nthi=as.integer(0)}
	else
	{
		nthi=as.integer(nth)
		if(nthi<=0)
			stop('Number of threads must be positive')
	}
	if((loglvi!=loglv)||(!is.null(rs)&&(rsi!=rs))||(!is.null(nth)&&(nthi!=nth)))
		stop('All parameters must be integer')
	if((loglvi<0)||(rsi<0))
		stop('All parameters must be nonnegative')
	if((loglvi>12))
		stop('Maximum log level is 12')
	findr.libload(loglvi,rsi,nthi)
}

findr.pijs_gassist_pv=function(dg,dt,dt2,na=NULL) {
	#Validity check
	if((length(dim(dg))!=2)||(length(dim(dt))!=2)||(length(dim(dt2))!=2))
		stop('Wrong input data dimensions. Check length(dim(x))==2 for x=dg,dt,dt2.')
	if((typeof(dg)!=findr.gtype)||(typeof(dt)!=findr.ftype)||(typeof(dt2)!=findr.ftype))
		stop(paste('Wrong input data types. Check typeof(dg) is ',findr.gtype,', and typeof(dt or dt2) is ',findr.ftype,sep=''))
	ng=dim(dg)[1]
	nt=dim(dt2)[1]
	ns=dim(dg)[2]
	if((dim(dt)[1]!=ng)||(dim(dt)[2]!=ns)||(dim(dt2)[2]!=ns))
		stop('Wrong input data shapes. Check dg and dt have same dimension, dt2 and dt have same column count')
	#Set nv value
	if(!is.null(na))
	{
		if(length(na)>1)
			stop('na should have length 1.')
		nai=as.integer(na)
		if(nai!=na)
			stop('Data type of na is not integer.')
		if(nai<as.integer(max(dg)))
			stop('Input number of alleles smaller than maximum genotype value.')
		nvx=as.integer(nai+1)
	}
	else
		nvx=as.integer(max(dg)+1)
	
	if(min(dg)<0)
		stop('Genotype values should not be negative.')
	#Output buffer
	ans=list('p1'=array(as.ftype(0),ng),'p2'=matrix(as.ftype(0),ng,nt),'p3'=matrix(as.ftype(0),ng,nt),'p4'=matrix(as.ftype(0),ng,nt),'p5'=matrix(as.ftype(0),ng,nt),'ret'=as.integer(0))
	ans=.C("external_R_pijs_gassist_pv",ng,nt,ns,
		 as.integer(dg),as.ftype(dt),as.ftype(dt2),
		 'p1'=ans$p1,'p2'=ans$p2,'p3'=ans$p3,'p4'=ans$p4,'p5'=ans$p5,
		 nvx,'ret'=ans$ret)
	if(ans$ret!=0)
		stop('C library execution failed.')
	ans=list('p1'=ans$p1,'p2'=ans$p2,'p3'=ans$p3,'p4'=ans$p4,'p5'=ans$p5)
	rownames(ans$p1)=rownames(dt)
	rownames(ans$p2)=rownames(dt)
	rownames(ans$p3)=rownames(dt)
	rownames(ans$p4)=rownames(dt)
	rownames(ans$p5)=rownames(dt)
	colnames(ans$p2)=rownames(dt2)
	colnames(ans$p3)=rownames(dt2)
	colnames(ans$p4)=rownames(dt2)
	colnames(ans$p5)=rownames(dt2)
	ans
}

findr.pijs_gassist=function(dg,dt,dt2,na=NULL,nodiag=FALSE) {
	#Validity check
	if(typeof(nodiag)!='logical')
		stop('Data type of nodiag is not logical.')
	if(length(nodiag)>1)
		stop('nodiag should have length 1.')
	if((length(dim(dg))!=2)||(length(dim(dt))!=2)||(length(dim(dt2))!=2))
		stop('Wrong input data dimensions. Check length(dim(x))==2 for x=dg,dt,dt2.')
	if((typeof(dg)!=findr.gtype)||(typeof(dt)!=findr.ftype)||(typeof(dt2)!=findr.ftype))
		stop(paste('Wrong input data types. Check typeof(g) is ',findr.gtype,', and typeof(dt or dt2) is ',findr.ftype,sep=''))
	ng=dim(dg)[1]
	nt=dim(dt2)[1]
	ns=dim(dg)[2]
	if((dim(dt)[1]!=ng)||(dim(dt)[2]!=ns)||(dim(dt2)[2]!=ns))
		stop('Wrong input data shapes. Check dg and dt have same dimension, dt2 and dt have same column count')
	#Set nv value
	if(!is.null(na))
	{
		if(length(na)>1)
			stop('na should have length 1.')
		nai=as.integer(na)
		if(nai!=na)
			stop('Data type of na is not integer.')
		if(nai<as.integer(max(dg)))
			stop('Input number of alleles smaller than maximum genotype value.')
		nvx=as.integer(nai+1)
	}
	else
		nvx=as.integer(max(dg)+1)
	
	if(min(dg)<0)
		stop('Genotype values should not be negative.')
	#set nodiag value
	if(nodiag)
		nd=as.integer(1)
	else
		nd=as.integer(0)
	#Output buffer
	ans=list('p1'=array(as.ftype(0),ng),'p2'=matrix(as.ftype(0),ng,nt),'p3'=matrix(as.ftype(0),ng,nt),'p4'=matrix(as.ftype(0),ng,nt),'p5'=matrix(as.ftype(0),ng,nt),'ret'=as.integer(0))
	ans=.C("external_R_pijs_gassist",ng,nt,ns,
		 as.integer(dg),as.ftype(dt),as.ftype(dt2),
		 'p1'=ans$p1,'p2'=ans$p2,'p3'=ans$p3,'p4'=ans$p4,'p5'=ans$p5,
		 nvx,nd,'ret'=ans$ret)
	if(ans$ret!=0)
		stop('C library execution failed.')
	ans=list('p1'=ans$p1,'p2'=ans$p2,'p3'=ans$p3,'p4'=ans$p4,'p5'=ans$p5)
	rownames(ans$p1)=rownames(dt)
	rownames(ans$p2)=rownames(dt)
	rownames(ans$p3)=rownames(dt)
	rownames(ans$p4)=rownames(dt)
	rownames(ans$p5)=rownames(dt)
	colnames(ans$p2)=rownames(dt2)
	colnames(ans$p3)=rownames(dt2)
	colnames(ans$p4)=rownames(dt2)
	colnames(ans$p5)=rownames(dt2)
	ans
}

findr.pij_gassist_any=function(dg,dt,dt2,name,na=NULL,nodiag=FALSE) {
	#Validity check
	if(typeof(nodiag)!='logical')
		stop('Data type of nodiag is not logical.')
	if(length(nodiag)>1)
		stop('nodiag should have length 1.')
	if((length(dim(dg))!=2)||(length(dim(dt))!=2)||(length(dim(dt2))!=2))
		stop('Wrong input data dimensions. Check length(dim(x))==2 for x=dg,dt,dt2.')
	if((typeof(dg)!=findr.gtype)||(typeof(dt)!=findr.ftype)||(typeof(dt2)!=findr.ftype))
		stop(paste('Wrong input data types. Check typeof(g) is ',findr.gtype,', and typeof(dt or dt2) is ',findr.ftype,sep=''))
	ng=dim(dg)[1]
	nt=dim(dt2)[1]
	ns=dim(dg)[2]
	if((dim(dt)[1]!=ng)||(dim(dt)[2]!=ns)||(dim(dt2)[2]!=ns))
		stop('Wrong input data shapes. Check dg and dt have same dimension, dt2 and dt have same column count')
	#Set nv value
	if(!is.null(na))
	{
		if(length(na)>1)
			stop('na should have length 1.')
		nai=as.integer(na)
		if(nai!=na)
			stop('Data type of na is not integer.')
		if(nai<as.integer(max(dg)))
			stop('Input number of alleles smaller than maximum genotype value.')
		nvx=as.integer(nai+1)
	}
	else
		nvx=as.integer(max(dg)+1)
	
	if(min(dg)<0)
		stop('Genotype values should not be negative.')
	#set nodiag value
	if(nodiag)
		nd=as.integer(1)
	else
		nd=as.integer(0)
	#Output buffer
	ans=list('p'=matrix(as.ftype(0),ng,nt),'ret'=as.integer(0))
	ans=.C(name,ng,nt,ns,
		 as.integer(dg),as.ftype(dt),as.ftype(dt2),
		 'p'=ans$p,nvx,nd,'ret'=ans$ret)
	if(ans$ret!=0)
		stop('C library execution failed.')
	ans=ans$p
	rownames(ans)=rownames(dt)
	colnames(ans)=rownames(dt2)
	ans
}

findr.pij_gassist=function(dg,dt,dt2,na=NULL,nodiag=FALSE) {
	findr.pij_gassist_any(dg,dt,dt2,"external_R_pij_gassist",na=na,nodiag=nodiag)
}

findr.pij_gassist_trad=function(dg,dt,dt2,na=NULL,nodiag=FALSE) {
	findr.pij_gassist_any(dg,dt,dt2,"external_R_pij_gassist_trad",na=na,nodiag=nodiag)
}

findr.pijs_cassist_pv=function(dc,dt,dt2) {
	#Validity check
	if((length(dim(dc))!=2)||(length(dim(dt))!=2)||(length(dim(dt2))!=2))
		stop('Wrong input data dimensions. Check length(dim(x))==2 for x=dc,dt,dt2.')
	if((typeof(dc)!=findr.ftype)||(typeof(dt)!=findr.ftype)||(typeof(dt2)!=findr.ftype))
		stop(paste('Wrong input data types. Check typeof(dc, dt, or dt2) is ',findr.ftype,sep=''))
	ng=dim(dc)[1]
	nt=dim(dt2)[1]
	ns=dim(dc)[2]
	if((dim(dt)[1]!=ng)||(dim(dt)[2]!=ns)||(dim(dt2)[2]!=ns))
		stop('Wrong input data shapes. Check dc and dt have same dimension, dt2 and dt have same column count')
	
	#Output buffer
	ans=list('p1'=array(as.ftype(0),ng),'p2'=matrix(as.ftype(0),ng,nt),'p3'=matrix(as.ftype(0),ng,nt),'p4'=matrix(as.ftype(0),ng,nt),'p5'=matrix(as.ftype(0),ng,nt),'ret'=as.integer(0))
	ans=.C("external_R_pijs_cassist_pv",ng,nt,ns,
		 as.ftype(dc),as.ftype(dt),as.ftype(dt2),
		 'p1'=ans$p1,'p2'=ans$p2,'p3'=ans$p3,'p4'=ans$p4,'p5'=ans$p5,
		 'ret'=ans$ret)
	if(ans$ret!=0)
		stop('C library execution failed.')
	ans=list('p1'=ans$p1,'p2'=ans$p2,'p3'=ans$p3,'p4'=ans$p4,'p5'=ans$p5)
	rownames(ans$p1)=rownames(dt)
	rownames(ans$p2)=rownames(dt)
	rownames(ans$p3)=rownames(dt)
	rownames(ans$p4)=rownames(dt)
	rownames(ans$p5)=rownames(dt)
	colnames(ans$p2)=rownames(dt2)
	colnames(ans$p3)=rownames(dt2)
	colnames(ans$p4)=rownames(dt2)
	colnames(ans$p5)=rownames(dt2)
	ans
}

findr.pijs_cassist=function(dc,dt,dt2,nodiag=FALSE) {
	#Validity check
	if(typeof(nodiag)!='logical')
		stop('Data type of nodiag is not logical.')
	if(length(nodiag)>1)
		stop('nodiag should have length 1.')
	#Validity check
	if((length(dim(dc))!=2)||(length(dim(dt))!=2)||(length(dim(dt2))!=2))
		stop('Wrong input data dimensions. Check length(dim(x))==2 for x=dc,dt,dt2.')
	if((typeof(dc)!=findr.ftype)||(typeof(dt)!=findr.ftype)||(typeof(dt2)!=findr.ftype))
		stop(paste('Wrong input data types. Check typeof(dc, dt, or dt2) is ',findr.ftype,sep=''))
	ng=dim(dc)[1]
	nt=dim(dt2)[1]
	ns=dim(dc)[2]
	if((dim(dt)[1]!=ng)||(dim(dt)[2]!=ns)||(dim(dt2)[2]!=ns))
		stop('Wrong input data shapes. Check dc and dt have same dimension, dt2 and dt have same column count')
	
	#set nodiag value
	if(nodiag)
		nd=as.integer(1)
	else
		nd=as.integer(0)
	#Output buffer
	ans=list('p1'=array(as.ftype(0),ng),'p2'=matrix(as.ftype(0),ng,nt),'p3'=matrix(as.ftype(0),ng,nt),'p4'=matrix(as.ftype(0),ng,nt),'p5'=matrix(as.ftype(0),ng,nt),'ret'=as.integer(0))
	ans=.C("external_R_pijs_cassist",ng,nt,ns,
		 as.ftype(dc),as.ftype(dt),as.ftype(dt2),
		 'p1'=ans$p1,'p2'=ans$p2,'p3'=ans$p3,'p4'=ans$p4,'p5'=ans$p5,
		 nd,'ret'=ans$ret)
	if(ans$ret!=0)
		stop('C library execution failed.')
	ans=list('p1'=ans$p1,'p2'=ans$p2,'p3'=ans$p3,'p4'=ans$p4,'p5'=ans$p5)
	rownames(ans$p1)=rownames(dt)
	rownames(ans$p2)=rownames(dt)
	rownames(ans$p3)=rownames(dt)
	rownames(ans$p4)=rownames(dt)
	rownames(ans$p5)=rownames(dt)
	colnames(ans$p2)=rownames(dt2)
	colnames(ans$p3)=rownames(dt2)
	colnames(ans$p4)=rownames(dt2)
	colnames(ans$p5)=rownames(dt2)
	ans
}

findr.pij_cassist_any=function(dc,dt,dt2,name,nodiag=FALSE) {
	#Validity check
	if(typeof(nodiag)!='logical')
		stop('Data type of nodiag is not logical.')
	if(length(nodiag)>1)
		stop('nodiag should have length 1.')
	#Validity check
	if((length(dim(dc))!=2)||(length(dim(dt))!=2)||(length(dim(dt2))!=2))
		stop('Wrong input data dimensions. Check length(dim(x))==2 for x=dc,dt,dt2.')
	if((typeof(dc)!=findr.ftype)||(typeof(dt)!=findr.ftype)||(typeof(dt2)!=findr.ftype))
		stop(paste('Wrong input data types. Check typeof(dc, dt, or dt2) is ',findr.ftype,sep=''))
	ng=dim(dc)[1]
	nt=dim(dt2)[1]
	ns=dim(dc)[2]
	if((dim(dt)[1]!=ng)||(dim(dt)[2]!=ns)||(dim(dt2)[2]!=ns))
		stop('Wrong input data shapes. Check dc and dt have same dimension, dt2 and dt have same column count')
		
	#set nodiag value
	if(nodiag)
		nd=as.integer(1)
	else
		nd=as.integer(0)
	#Output buffer
	ans=list('p'=matrix(as.ftype(0),ng,nt),'ret'=as.integer(0))
	ans=.C(name,ng,nt,ns,
		 as.ftype(dc),as.ftype(dt),as.ftype(dt2),
		 'p'=ans$p,nd,'ret'=ans$ret)
	if(ans$ret!=0)
		stop('C library execution failed.')
	ans=ans$p
	rownames(ans)=rownames(dt)
	colnames(ans)=rownames(dt2)
	ans
}

findr.pij_cassist=function(dc,dt,dt2,nodiag=FALSE) {
	findr.pij_cassist_any(dc,dt,dt2,"external_R_pij_cassist",nodiag=nodiag)
}

findr.pij_cassist_trad=function(dc,dt,dt2,nodiag=FALSE) {
	findr.pij_cassist_any(dc,dt,dt2,"external_R_pij_cassist_trad",nodiag=nodiag)
}

findr.pij_rank_pv=function(dt,dt2) {
	#Validity check
	if((length(dim(dt))!=2)||(length(dim(dt2))!=2))
		stop('Wrong input data dimensions. Check length(dim(x))==2 for x=dt,dt2.')
	if((typeof(dt)!=findr.ftype)||(typeof(dt2)!=findr.ftype))
		stop(paste('Wrong input data types. Check typeof(dt or dt2) is ',findr.ftype,sep=''))
	ng=dim(dt)[1]
	nt=dim(dt2)[1]
	ns=dim(dt)[2]
	if(dim(dt2)[2]!=ns)
		stop('Wrong input data shapes. Check t and t2 have same column count')

	#Output buffer
	ans=list('p'=matrix(as.ftype(0),ng,nt),
					 'ret'=as.integer(0))
	ans=.C("external_R_pij_rank_pv",ng,nt,ns,
		 as.ftype(dt),as.ftype(dt2),
		 'p'=ans$p,'ret'=ans$ret)
	if(ans$ret!=0)
		stop('C library execution failed.')
	rownames(ans$p)=rownames(dt)
	colnames(ans$p)=rownames(dt2)
	ans=ans$p
	ans
}

findr.pij_rank=function(dt,dt2,nodiag=FALSE) {
	#Validity check
	if(typeof(nodiag)!='logical')
		stop('Data type of nodiag is not logical.')
	if((length(dim(dt))!=2)||(length(dim(dt2))!=2))
		stop('Wrong input data dimensions. Check length(dim(x))==2 for x=dt,dt2.')
	if((typeof(dt)!=findr.ftype)||(typeof(dt2)!=findr.ftype))
		stop(paste('Wrong input data types. Check typeof(dt or dt2) is ',findr.ftype,sep=''))
	ng=dim(dt)[1]
	nt=dim(dt2)[1]
	ns=dim(dt)[2]
	if(dim(dt2)[2]!=ns)
		stop('Wrong input data shapes. Check t and t2 have same column count')

	#set nodiag value
	if(nodiag)
		nd=as.integer(1)
	else
		nd=as.integer(0)
	#Output buffer
	ans=list('p'=matrix(as.ftype(0),ng,nt),
					 'ret'=as.integer(0))
	ans=.C("external_R_pij_rank",ng,nt,ns,
		 as.ftype(dt),as.ftype(dt2),
		 'p'=ans$p,nd,'ret'=ans$ret)
	if(ans$ret!=0)
		stop('C library execution failed.')
	rownames(ans$p)=rownames(dt)
	colnames(ans$p)=rownames(dt2)
	ans=ans$p
	ans
}

findr.netr_one_greedy=function(dp,namax=NULL,nimax=NULL,nomax=NULL) {
	#Validity check
	if((length(dim(dp))!=2))
		stop('Wrong input data dimensions. Check length(dim(dp))==2.')
	if((typeof(dp)!=findr.ftype))
		stop(paste('Wrong input data types. Check typeof(dp) is ',findr.ftype,sep=''))
	nt=dim(dp)[1]
	if(dim(dp)[2]!=nt)
		stop('Wrong input data shapes. Check dp is a square matrix.')
	if(!is.null(rownames(dp))&!is.null(colnames(dp))&!all(rownames(dp)==colnames(dp)))
		stop('Different row and column names for dp. Ensure dp is a square matrix.')
	
	if(!is.null(namax))
	{
		if(length(namax)>1)
			stop('Length of namax must be 1.')
		nav=as.integer(namax)
		if(nav!=namax)
			stop('Data type of namax is not integer.')
		if(nav<=0)
			stop('namax must be positive')
	}
	else
	{nav=0}
	
	if(!is.null(nimax))
	{
		if(length(nimax)>1)
			stop('Length of namax must be 1.')
		niv=as.integer(nimax)
		if(niv!=nimax)
			stop('Data type of nimax is not integer.')
		if(niv<=0)
			stop('namax must be positive')
	}
	else
	{niv=0}

	if(!is.null(nomax))
	{
		if(length(nomax)>1)
			stop('Length of namax must be 1.')
		nov=as.integer(nomax)
		if(nov!=nomax)
			stop('Data type of nomax is not integer.')
		if(nov<=0)
			stop('namax must be positive')
	}
	else
	{nov=0}
	
	#Output buffer
	ans=list('net'=matrix(FALSE,nt,nt),'ret'=as.integer(0))
	ans=.C("external_R_netr_one_greedy",nt,as.ftype(dp),nav,niv,nov,'net'=ans$net,'ret'=ans$ret)
	if(ans$ret!=0)
		stop('C library execution failed.')
	rownames(ans$net)=rownames(dp)
	colnames(ans$net)=rownames(dp)
	ans=ans$net
	ans
}



























