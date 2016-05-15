# Copyright 2016 Lingfei Wang
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

findr.lib=function(loglv=6,rs=0,nth=0){
	loglvi=as.integer(loglv)
	rsi=as.integer(rs)
	nthi=as.integer(nth)
	if((loglvi!=loglv)||(rsi!=rs)||(nthi!=nth))
		stop('All parameters must be integer')
	if((loglvi<0)||(rsi<0)||(nthi<0))
		stop('All parameters must be nonnegative')
	findr.libload(loglvi,rsi,nthi)
}

findr.pijs_gassist_a=function(dg,dt,dt2,na=NULL,nodiag=FALSE) {
	#Validity check
	if(typeof(nodiag)!='logical')
		stop('Data type of nodiag is not logical.')
	if(!is.null(na))
	{
		nai=as.integer(na)
		if(nai!=na)
			stop('Data type of na is not integer.')
	}
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
	if(is.null(na)||(na<=0))
		nvx=as.integer(max(dg)+1)
	else
		nvx=nai+1
	if((min(dg)<0)||(max(dg)>=nvx))
		stop('Genotype values do not fall into nv range.')
	#set nodiag value
	if(nodiag)
		nd=as.integer(1)
	else
		nd=as.integer(0)
	#Output buffer
	ans=list('p1'=array(as.ftype(0),ng),'p2b'=matrix(as.ftype(0),ng,nt),'p2c'=matrix(as.ftype(0),ng,nt),'p3'=matrix(as.ftype(0),ng,nt),'ret'=as.integer(0))
	ans=.C("external_R_pijs_gassist_a",ng,nt,ns,
		 as.integer(dg),as.ftype(dt),as.ftype(dt2),
		 'p1'=ans$p1,'p2b'=ans$p2b,'p2c'=ans$p2c,'p3'=ans$p3,
		 nvx,nd,'ret'=ans$ret)
	if(ans$ret!=0)
		stop('C library execution failed.')
	ans=list('p1'=ans$p1,'p2b'=ans$p2b,'p2c'=ans$p2c,'p3'=ans$p3)
	rownames(ans$p1)=rownames(dt)
	rownames(ans$p2b)=rownames(dt)
	rownames(ans$p2c)=rownames(dt)
	rownames(ans$p3)=rownames(dt)
	colnames(ans$p2b)=rownames(dt2)
	colnames(ans$p2c)=rownames(dt2)
	colnames(ans$p3)=rownames(dt2)
	ans
}


findr.pijs_gassist_tot=function(dg,dt,dt2,na=NULL,nodiag=FALSE) {
	#Validity check
	if(typeof(nodiag)!='logical')
		stop('Data type of nodiag is not logical.')
	if(!is.null(na))
	{
		nai=as.integer(na)
		if(nai!=na)
			stop('Data type of na is not integer.')
	}
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
	if(is.null(na)||(na<=0))
		nvx=as.integer(max(dg)+1)
	else
		nvx=nai+1
	if((min(dg)<0)||(max(dg)>=nvx))
		stop('Genotype values do not fall into nv range.')
	#set nodiag value
	if(nodiag)
		nd=as.integer(1)
	else
		nd=as.integer(0)
	#Output buffer
	ans=list('p1'=array(as.ftype(0),ng),'p2b'=matrix(as.ftype(0),ng,nt),'p2c'=matrix(as.ftype(0),ng,nt),'p3'=matrix(as.ftype(0),ng,nt),'ret'=as.integer(0))
	ans=.C("external_R_pijs_gassist_tot",ng,nt,ns,
		 as.integer(dg),as.ftype(dt),as.ftype(dt2),
		 'p1'=ans$p1,'p2b'=ans$p2b,'p2c'=ans$p2c,'p3'=ans$p3,
		 nvx,nd,'ret'=ans$ret)
	if(ans$ret!=0)
		stop('C library execution failed.')
	ans=list('p1'=ans$p1,'p2b'=ans$p2b,'p2c'=ans$p2c,'p3'=ans$p3)
	rownames(ans$p1)=rownames(dt)
	rownames(ans$p2b)=rownames(dt)
	rownames(ans$p2c)=rownames(dt)
	rownames(ans$p3)=rownames(dt)
	colnames(ans$p2b)=rownames(dt2)
	colnames(ans$p2c)=rownames(dt2)
	colnames(ans$p3)=rownames(dt2)
	ans
}

findr.pij_rank_a=function(dt,dt2,nodiag=FALSE) {
	#Validity check
	if(typeof(nodiag)!='logical')
		stop('Data type of nodiag is not logical.')
	if((length(dim(dt))!=2)||(length(dim(dt2))!=2))
		stop('Wrong input data dimensions. Check length(dim(x))==2 for x=dt,dt2.')
	if((typeof(dt)!=findr.ftype)||(typeof(dt2)!=findr.ftype))
		stop(paste('Wrong input data types. Check typeof(dt or t2) dis ',findr.ftype,sep=''))
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
	ans=.C("external_R_pij_rank_a",ng,nt,ns,
		 as.ftype(dt),as.ftype(dt2),
		 'p'=ans$p,nd,'ret'=ans$ret)
	if(ans$ret!=0)
		stop('C library execution failed.')
	rownames(ans$p)=rownames(dt)
	colnames(ans$p)=rownames(dt2)
	ans
}



























