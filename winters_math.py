# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 08:21:55 2018

@author: bwinters
"""

#==============================================================================
# This is a matrix class with row and column opporations    
#==============================================================================

class matrix:
    def __init__(self,List,rows, cols):
        self.list=List
        self.olist=List.copy()
        self.rows=rows
        self.cols=cols
        if rows*cols!=len(List):
            raise ValueError('rows times columns does not equal number of elements')
        return None
    
    @classmethod
    def from_file(cls, filename):
        rows=0
        List=[]
        mfile=open(filename,'r')
        for line in mfile:
            rows=rows+1
            line=line[0:-1]
            row=line.split(',')
            row=[int(item) for item in row]
            List=List+row
        cols=len(List)//rows
        return cls(List,rows, cols)
    
    @classmethod
    def copy(cls,M):
        List=M.list.copy()
        N=matrix(List,M.rows,M.cols)
        return N
        
        
    def get_entry(self,i,j):
        k=i*self.cols+j
        return self.list[k]
        
    def set_entry(self,i,j,entry):
        if i>self.rows:
            raise ValueError('not in row range')
        if j>self.cols:
            raise ValueError('not in column range')
        k=i*self.cols+j
        self.list[k]=entry
        return None
    def invert(self):
        M=matrix(self.list.copy(),self.rows,self.cols)
        if M.rows!=M.cols:
            raise ValueError('Must be a square matrix')
        I=M.id_pre()
        for start in range(M.rows):
            i,j=M.find_pivot(start)
            if j>start:
                raise ValueError('Matrix is not invertable')
            M.swap_rows(start,i)
            I.swap_rows(start,i)
            entry=M.get_entry(start,start)
            if entry!=0:
                M.row_mult(start,1/entry)
                I.row_mult(start,1/entry)
            else:
                raise ValueError('Matrix is not invertable')
            for k in range(start+1,M.rows):
                entry=M.get_entry(k,start)
                M.add_rows(start,k,-entry)
                I.add_rows(start,k,-entry)
            for pivot in range(M.rows-1,0,-1):
                for r in range(pivot):
                    entry=M.get_entry(r,pivot)
                    M.add_rows(pivot,r,-entry)
                    I.add_rows(pivot,r,-entry)

        return I
        
    def swap_cols(self,i,j):
        R=self.rows
        for k in range(R):
            a=self.get_entry(k,i)
            b=self.get_entry(k,j)
            self.set_entry(k,i,b)
            self.set_entry(k,j,a)
        return None
        
    def swap_rows(self,i,j):
        C=self.cols
        for k in range(C):
            a=self.get_entry(i,k)
            b=self.get_entry(j,k)
            self.set_entry(i,k,b)
            self.set_entry(j,k,a)
        return None
        
    def add_cols(self,i,j,scalar):
        R=self.rows
        for k in range(R):
            a=self.get_entry(k,i)
            b=self.get_entry(k,j)
            self.set_entry(k,j,scalar*a+b)
        return None
    def row_mult(self,r,scalar):
        if r not in range(self.rows):
            raise ValueError('Row out of range')
        for i in range(self.cols):
            entry=self.get_entry(r,i)
            self.set_entry(r,i,scalar*entry)
        return None

    def col_mult(self,c,scalar):
        if c not in range(self.cols):
            raise ValueError('Column out of range')
        for i in range(self.rows):
            entry=self.get_entry(i,c)
            self.set_entry(i,c,scalar*entry)
        return None        
    def add_rows(self,i,j,scalar):
        C=self.cols
        for k in range(C):
            a=self.get_entry(i,k)
            b=self.get_entry(j,k)
            self.set_entry(j,k,scalar*a+b)
        return None        

    def combine_rows(self,i,j,a,b,c,d):
        """replaces rows i,j with linear combinations
        of row i and row j"""
        if i not in range(self.rows):
            raise ValueError('row number out of range')
        if j not in range(self.rows):
            raise ValueError('row number out of range')
        C=self.cols
        for k in range(C):
            entry_i=a*self.get_entry(i,k)+b*self.get_entry(j,k)
            entry_j=c*self.get_entry(i,k)+d*self.get_entry(j,k)
            self.set_entry(i,k,entry_i)
            self.set_entry(j,k,entry_j)
        return None
        
    def combine_columns(self,i,j,a,b,c,d):
        """Replaces column i,j with a linear combinations 
        of column i and column j"""
        if i not in range(self.cols):
            raise ValueError('column number out of range')
        if j not in range(self.cols):
            raise ValueError('column number out of range')
        R=self.rows
        for k in range(R):
            entry_i=a*self.get_entry(k,i)+b*self.get_entry(k,j)
            entry_j=c*self.get_entry(k,i)+d*self.get_entry(k,j)
            self.set_entry(k,i,entry_i)
            self.set_entry(k,j,entry_j)
        return None

        
    def show(self):
        R=self.rows
        C=self.cols
        out=[]
        for i in range(R):
            out.append(self.list[i*C:(i+1)*C])
        return out
    
    def multiply(self,M,N):
        R=M.rows
        C1=M.cols
        R2=N.rows
        C=N.cols
        List=(R*C)*[0]
        if C1!=R2:
            raise ValueError('cols of first must equal rows of second')
        for i in range(R):
            for j in range(C):
                S=0
                for k in range(C1):
                   S=S+M.get_entry(i,k)*N.get_entry(k,j) 
                List[i*C+j]=S
        P=matrix(List,R,C)
        return P
        
    def id_post(self):
        C=self.cols
        List=[]
        for i in range(C):
            for j in range(C):
                if i==j:
                    List.append(1)
                else:
                    List.append(0)
        P=matrix(List,C,C)
        return P
    def id_pre(self):
        R=self.rows
        List=[]
        for i in range(R):
            for j in range(R):
                if i==j:
                    List.append(1)
                else:
                    List.append(0)
        P=matrix(List,R,R)
        return P

    def sn_compact(self):
        M=matrix(self.list.copy(),self.rows,self.cols)
        crank,rrank,brank,z_col=M.ranks()
        T=[]
        ones=0
        entry=M.get_entry(ones,ones)
        while entry==1:
            ones+=1
            if ones<M.rows and ones<M.cols:
                entry=M.get_entry(ones,ones)
            else:
                return (ones, T,crank,rrank,brank,ones+len(T),z_col)
        count=ones
        while entry!=0:
            count+=1
            T.append(entry)
            if count<M.rows and count<M.cols:
                entry=M.get_entry(count,count)
            else:
                break
        return (ones, T,crank,rrank,brank,ones+len(T),z_col)
        
    def smith_normal(self):
        top=0
        self.U=self.id_pre()
        self.V=self.id_post()
        K,L=self.find_pivot(top)
        if K==-1:
            return None
        while K>=0:
            if top%20==0:
                print('Now doing row ',top, ' of ',self.rows )
            self.swap_rows(K,top)
            self.U.swap_rows(K,top)
            self.swap_cols(L,top)
            self.V.swap_cols(L,top)
            while not self.pivot_zeroed(top):
                self.improve_pivot_down(top)
                self.improve_pivot_right(top)
            entry=self.get_entry(top,top)
            if entry<0:
                self.row_mult(top,-1)
                self.U.row_mult(top,-1)
                #self.V.col_mult(top,-1)
            top+=1
            if top<self.rows:
                K,L=self.find_pivot(top)
            else:
                return None   
        return None
    
    def find_pivot(self,i):
        K=i
        L=i
        entry=self.get_entry(K,L)        
        if entry!=0:
            return (K,L)
        while entry==0 and L<self.cols:
            while entry==0 and K<self.rows:
                K+=1
                if K<self.rows:
                    entry=self.get_entry(K,L)
                    if entry!=0:
                        return (K,L)
            K=i
            L+=1
            if L<self.rows:
                entry=self.get_entry(K,L)
                if entry!=0:
                    return (K,L)
        return (-1,-1)
    def pivot_zeroed(self,i):
        for k in range(i+1,self.rows):
            if self.get_entry(k,i)!=0:
                return False
        for k in range(i+1, self.cols):
            if self.get_entry(i,k)!=0:
                return False
        return True
    def improve_pivot_down(self,i):
        if i not in range(self.rows):
            raise ValueError('i is not in row range')
        for j in range(i+1,self.rows):
            a=self.get_entry(i,i)
            b=self.get_entry(j,i)
            if not self.divides(a,b):
                if self.divides(b,a):
                    self.swap_rows(i,j)
                    self.U.swap_rows(i,j)
                else:
                    self.refine_row(i,j)
            a=self.get_entry(i,i)
            b=self.get_entry(j,i)
            scalar=(-1)*(b//a)
            self.add_rows(i,j,scalar)            
            self.U.add_rows(i,j,scalar)
        return None
    def improve_pivot_right(self,i):
        for j in range(i+1,self.cols):
            a=self.get_entry(i,i)
            b=self.get_entry(i,j)
            if not self.divides(a,b):
                if self.divides(b,a):
                    self.swap_cols(i,j)
                    self.V.swap_cols(i,j)
                else:
                    self.refine_col(i,j)
            a=self.get_entry(i,i)
            b=self.get_entry(i,j)
            scalar=(-1)*(b//a)
            self.add_cols(i,j,scalar)            
            self.V.add_cols(i,j,scalar)
        return None       
        
        
    def refine_row(self,i,j):
        a=self.get_entry(i,i)
        b=self.get_entry(j,i)
        beta,sigma,tau=self.gcd(a,b)
        alpha=a//beta
        gamma=b//beta
        self.combine_rows(i,j,sigma,tau,-gamma,alpha)
        self.U.combine_rows(i,j,sigma,tau,-gamma,alpha)
        return None
    
        
    def refine_col(self,i,j):
        a=self.get_entry(i,i)
        b=self.get_entry(i,j)
        beta,sigma,tau=self.gcd(a,b)
        alpha=a//beta
        gamma=b//beta
        self.combine_columns(i,j,sigma,tau,-gamma,alpha)
        self.V.combine_columns(i,j,sigma,tau,-gamma,alpha)
        return None
        
    def refine_col_mat(self,i,j,sigma,tau):
        P=self.id_post()
        P.set_entry(i,i,sigma)
        P.set_entry(j,i,tau)
        return P
        
    def sweep_down(self,i):
        R=self.rows
        for k in range(i+1,R):
            if self.get_entry(k,i):
                self.add_rows_2(i,k)
        return None
    def sweep_right(self,i):
        C=self.cols
        for k in range(i+1,C):
            if self.get_entry(i,k)!=0:
                self.add_cols_2(i,k)
        return None
        
    def zero_cols(self):
        count=0
        C=self.cols
        R=self.rows
        for k in range(C):
            flag=True
            for i in range(R):
                entry=self.get_entry(i,k)
                if entry!=0:
                    flag=False
                    break
            if flag:
                count+=1
        return count
           
        
    def ranks(self):
        #This returns the rank of Zp and of B_{p-1}
        self.smith_normal()
        crank=self.cols
        rrank=self.rows
        brank=self.last_nonzero_diag()+1
        z_col=self.zero_cols()
        return (crank,rrank,brank,z_col)
        
    def last_nonzero_diag(self):
        i=0
        entry=self.get_entry(i,i)
        if entry==0:
            return -1
        while entry!=0:
            i+=1
            if i<self.rows and i<self.cols:
                entry=self.get_entry(i,i)
                if entry==0:
                    return i-1
            else:
                return i-1
        return i-1
    def div(self,p,q):
        a=p//q
        r=p%q
        return(a,r)
    def divides(self,p,q):
        r=abs(q)%abs(p)
        if r==0:
            return True
        else:
            return False
        return None            
    def gcd(self,m,n):
        if self.divides(m,n):
            return (abs(m),m//abs(m),0)
        if self.divides(n,m):
            return (abs(n),0,n//abs(n))
        p,q=m,n
        a,r=self.div(p,q)
        A=[a]
        while r!=0:
            p,q=q,r
            a,r=self.div(p,q)
            A.append(a)
        A.pop()
        s0=1
        t0=-A.pop()
        while A:
            s1=t0
            t1=s0-t0*A.pop()
            s0,t0=s1,t1
        return q,s0,t0    
        

        
#==============================================================================
# These deal with the prime order function
#==============================================================================

primes_cache={}

def prime_order(n):
    """This will give the order of a prime in the increasing list
    of primes."""
    file=open('primes_less_than_100000.csv','r')
    primes=[]
    if n in primes_cache:
        return primes_cache[n]
    for line in file:
        row=line.split(',')
        row=[int(item) for item in row]
        primes=primes+row
    for i in range(len(primes)):
        primes_cache[i]=primes[i]
    return primes[n]

def vertex_order(n,vertices):
    """This gives the order of a vertex in the increasing list of
    vertices."""
    if n not in vertices:
        raise ValueError('n is not a vertex')
        return None
    m=vertices.index(n)
    return m
    
