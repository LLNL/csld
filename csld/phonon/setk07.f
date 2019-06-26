      subroutine setk07(na,nb,nc,a,b,c,ptk,nptk,idef,ntet,nkmax,ntmax)          
c     set the k-points in one-heigth the reciprocal cell for a                  
c     simple orthorhombic lattice with parameters a, b, c                       
c     symmetry is d2h                                                           
      implicit real*8(a-h,o-z)                                                  
      real*4 avol                                                               
      dimension ptk(4,nkmax),idef(5,ntmax)                                      
      equivalence (ivol,avol)                                                   
      pi = 3.141592653589793238d0                                               
      if(na.le.0.or.nb.le.0.or.nc.le.0) goto 97                                 
      if(a.le.0.0d0 .or. b.le.0.0d0 .or. c.le.0.0d0) goto 98                    
      nptk = (na+1)*(nb+1)*(nc+1)                                               
      if(nptk.gt.nkmax) stop '*** <SETK07> nptk exceeds nkmax ***'              
      ntet = 6*na*nb*nc                                                         
      if(ntet.gt.ntmax) stop '*** <SETK07> ntet exceeds ntmax ***'              
c *** set the k-points                                                          
      ak = pi/a/na                                                              
      bk = pi/b/nb                                                              
      ck = pi/c/nc                                                              
      write(6,100) nptk,ntet,na*ak,nb*bk,nc*ck                                  
      nx1=na+1                                                                  
      ny1=nb+1                                                                  
      nz1=nc+1                                                                  
      w = 1.0d0/(na*nb*nc)                                                      
      nptk=0                                                                    
      do 5 i=1,nx1                                                              
      do 5 j=1,ny1                                                              
      do 5 k=1,nz1                                                              
c        nptk = (i-1)*ny1*nz1 + (j-1)*nz1 + k                                   
      wk = w                                                                    
      if(i.eq.1) wk = wk/2.0d0                                                  
      if(j.eq.1) wk = wk/2.0d0                                                  
      if(k.eq.1) wk = wk/2.0d0                                                  
      if(i.eq.nx1) wk = wk/2.0d0                                                
      if(j.eq.ny1) wk = wk/2.0d0                                                
      if(k.eq.nz1) wk = wk/2.0d0                                                
      nptk=nptk+1                                                               
      ptk(1,nptk)=(i-1)*ak                                                      
      ptk(2,nptk)=(j-1)*bk                                                      
      ptk(3,nptk)=(k-1)*ck                                                      
      ptk(4,nptk) = wk                                                          
    5 continue                                                                  
c *** define the tetrahedra                                                     
      ix=ny1*nz1                                                                
      ntet=0                                                                    
      i7=0                                                                      
      do 10 i=1,na                                                              
      do 10 j=1,nb                                                              
      i7=(i-1)*ix+(j-1)*nz1                                                     
      do 10 k=1,nc                                                              
      i7=i7+1                                                                   
      i6=i7+ix                                                                  
      i2=i6+nz1                                                                 
      i1=i2+1                                                                   
      ntet=ntet+1                                                               
      idef(1,ntet)=i7                                                           
      idef(2,ntet)=i6                                                           
      idef(3,ntet)=i2                                                           
      idef(4,ntet)=i1                                                           
      i8=i7+1                                                                   
      i5=i6+1                                                                   
      ntet=ntet+1                                                               
      idef(1,ntet)=i7                                                           
      idef(2,ntet)=i6                                                           
      idef(3,ntet)=i5                                                           
      idef(4,ntet)=i1                                                           
      ntet=ntet+1                                                               
      idef(1,ntet)=i7                                                           
      idef(2,ntet)=i8                                                           
      idef(3,ntet)=i5                                                           
      idef(4,ntet)=i1                                                           
      i3=i7+nz1                                                                 
      i4=i3+1                                                                   
      ntet=ntet+1                                                               
      idef(1,ntet)=i7                                                           
      idef(2,ntet)=i3                                                           
      idef(3,ntet)=i2                                                           
      idef(4,ntet)=i1                                                           
      ntet=ntet+1                                                               
      idef(1,ntet)=i7                                                           
      idef(2,ntet)=i3                                                           
      idef(3,ntet)=i4                                                           
      idef(4,ntet)=i1                                                           
      ntet=ntet+1                                                               
      idef(1,ntet)=i7                                                           
      idef(2,ntet)=i8                                                           
      idef(3,ntet)=i4                                                           
      idef(4,ntet)=i1                                                           
   10 continue                                                                  
      avol=1.d0/dfloat(ntet)                                                    
      do 15 it=1,ntet                                                           
   15 idef(5,it)=ivol                                                           
      return                                                                    
   97 write(6,101)                                                              
      goto 99                                                                   
   98 write(6,102)                                                              
   99 stop                                                                      
  100 format(' sampling the 8th part of a rectangular parallelepiped'/          
     .1x,i5,' k-points',i7,' tetrahedra'/                                       
     .' kxmax =',d11.4,'  kymax =',d11.4,'  kzmax =',d11.4)                     
  101 format(' *** <SETK07> na, nb or nc is not positive ***')                  
  102 format(' *** <SETK07> a, b or c is not positive ***')                     
      end
