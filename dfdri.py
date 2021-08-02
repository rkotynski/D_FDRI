#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Created on 27 Jul 2021

This D-FDRI code in Python enables to generate a measurement matrix with sampling patterns and a reconstruction matrix for singlepixel imaging (SPI)
The description of D-FDRI is included in the following paper. If you find the code useful, please cite this work:

Anna Pastuszczak, Rafał Stojek, Piotr Wróbel, and Rafał Kotyński, "Differential real-time single-pixel imaging with Fourier domain regularization - applications to VIS-IR imaging and polarization imaging,"
Optics Express, https://doi.org/10.1364/OE.433199 (open accesss)


@author: anna.pastuszczak@fuw.edu.pl, rafal.stojek@fuw.edu.pl, piotr.wrobel@fuw.edu.pl, rafal.kotynski@fuw.edu.pl


  GPL LICENSE INFORMATION

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


import numpy as np
from scipy import fft,random,linalg,io

class DFDRI:
    def __init__(self,μ=0.5,ϵ=1e-7,p=1,m=7,dim=(256,256),CR=0.03,verbose=True):
        self.μ=μ
        self.ϵ=ϵ
        self.p=p
        self.m=m
        self.CR=CR
        self.dim=dim
        self.verbose=verbose
        self.M_dct=None
        self.Mbin=None
        self.Pg=None
        if self.verbose:
            print('Differential Fourier Domain Regularized Inversion (D-FDRI)')
            print('This is a python module for compressive differential Single Pixel Imaging (SPI).')
            print('It calculates binary sampling patterns (a binary measurement matrix) and reconstructs images from compressive measurements at a very high speed.')
            print('For details please see https://doi.org/10.1364/OE.433199 (open access)\nIf you find this code useful, please cite the paper in your work.\n\n')

            print('Default parameters:')
            print(f'dim={dim}\t-image resolution')
            print(f'p={p}\t\t\t\t-order of the finite difference operator')
            print(f'p+m={p+m}\t\t\t-number of binary patterns used to measure the 0th spatial freq.')
            print(f'μ={μ}, ϵ={ϵ}\t-FDRI parameters')
            print(f'CR={CR}\t\t\t-compression ratio')
        
    def differential_operator(self,p=None):    
        '''
        Returns a function which is a differential operator of order p (p=1 or p=2)
        The operator takes a matrix as its first argument
        If the second argument ax=1, then it operates on the matrix from the right side
        '''
        if p is None:
            p=self.p
        if p==2:
           return lambda Y,right=False: np.hstack((np.zeros((Y.shape[0],1)) ,Y , np.zeros((Y.shape[0],1))))-0.5* np.hstack((np.zeros((Y.shape[0],2)) ,Y )) -0.5* np.hstack((Y,np.zeros((Y.shape[0],2))  ))  if right else   Y[1:-1,:]-0.5*Y[:-2,:]-0.5*Y[2:,:] 
        elif p==1:
           return lambda Y,right=False: np.hstack((np.zeros((Y.shape[0],1)),Y ))- np.hstack((Y,np.zeros((Y.shape[0],1)) ))   if right else  np.diff(Y,axis=0)
        else:
            print('Unsupported order p')
            
    
    def auxiliary_matrix_a(self,p=None, m=None,use_precalculated=True,maxiter=5000):
        if p is None:
            p=p=self.p
        if m is None:
            m=m=self.m        
        assert((m&1)and (m>1))   # m must be odd  
        if use_precalculated:
            A_precalculated={(2,7):[[0, 0, 0, 0, 1, 1, 1], [1, 0, 0, 1, 1, 0, 0], [1, 1, 0, 0, 1, 0, 0],
                                    [1, 1, 1, 1, 0, 0, 0], [1, 0, 0, 1, 0, 1, 0], [0, 1, 0, 1, 1, 0, 1],
                                    [0, 0, 1, 0, 1, 1, 1], [0, 0, 1, 0, 0, 1, 1], [0, 1, 1, 1, 0, 0, 0]],
                             (2,5):[[1, 0, 1, 1, 0], [1, 1, 0, 1, 0], [1, 0, 0, 1, 0], [0, 0, 1, 1, 0],
                                    [1, 1, 1, 0, 0], [0, 1, 0, 0, 1], [1, 0, 0, 1, 1]],
                             (2,3):[[0, 1, 1],[0, 0, 1], [0, 1, 0], [1, 0, 0], [1, 1, 0]],
                             (1,7):[[0, 1, 0, 0, 1, 1, 1], [1, 1, 0, 0, 0, 0, 1], [0, 0, 1, 1, 1, 0, 0],
                                     [1, 0, 1, 1, 0, 0, 1],[0, 1, 1, 0, 0, 1, 0], [0, 1, 0, 0, 1, 1, 0],
                                     [1, 1, 1, 0, 0, 0, 1], [0, 0, 1, 0, 1, 1, 1]],
                             (1,5):[[1, 1, 0, 1, 0], [1, 0, 1, 1, 0],[1, 1, 0, 0, 0],[0, 1, 0, 1, 1],
                                    [1, 1, 1, 0, 0],[1, 0, 1, 0, 1]],
                             (1,3):[[0, 1, 1], [1, 1, 0], [1, 0, 0], [0, 1, 0]]}
            if (p,m) in A_precalculated.keys():
                A=np.array(A_precalculated[p,m])
                if self.verbose:
                    print(f'Using the precalculated differential DC-decomposition matrix A\n{A}')
                return A
        if self.verbose:
            print(f'Preparing the differential DC-decomposition matrix A (m={m},p={p})')
        DIFF=self.differential_operator(p)              
    
        A=[]
        for i in range(2**m):                
            b=np.binary_repr(i,m)
            v=np.array([int(b[j]) for j in range(m)])
            if v.sum()==m//2 or v.sum()==(m//2)+1:
                A.append(v)
        A0=np.array(A)
        ntst=maxiter
        Abest=None
        stdbest=np.inf        
        while(True):            
            prm=np.random.permutation(A0.shape[0])[:m+p]
            M=DIFF(A0[prm,:])
            r=np.linalg.matrix_rank(M)
            if r==m:
                A=A0[prm,:]
                ntst-=1
                coef=np.linalg.inv(DIFF(A).T).sum(axis=1)
                stdnew=(DIFF(coef.reshape((1,-1)),ax=1)**2).sum()
                if stdnew<stdbest: # select the matrix giving the most uniform distribution 
                    stdbest=stdnew
                    Abest=np.array(A)
                    if self.verbose:
                        print(f'iter={maxiter-ntst}, crit={stdbest}')
                if ntst<=0:
                    A=Abest
                    if self.verbose:
                        print(f'Finished calculating the auxiliary matrix\n A={A}')
                    return A
    
    
    def dct_sampling_functions(self,dim=None,rows=None,CR=0.03):
        '''
        Create a matrix M with rows containg low-frequency continuous-valued 2d DCT functions
        The DCT functions are stored in rows of matrix M, and their locations in the 2D DCT basis are
        indicated by the logical matrix SM    
        The compression ratio CR may be passed to the function instead of the number of rows
        If binarize is True then a binarization with a randomly varied threshold is applied to M
        Rhe range of threshold values is governed by the value of m
        '''
        if dim is None:
            dim=self.dim
        cols=np.prod(dim)
        if rows is None:
            rows=int(round(cols*CR))
        InvTransform=lambda x:fft.idctn(x,norm='ortho')
        (x,y)=np.meshgrid(range(dim[1]),range(dim[0]))
        Avg=1/(x+y+1e-7) # selection function for the 2d DCT basis
        Avg[0,0]=0 # the 0th spatial frequency is excluded
        I=np.argsort(-Avg.reshape(-1))[:rows]
        SM=np.zeros(Avg.size,dtype=bool) # selection matrix
        SM[I]=True
        P=np.flatnonzero(SM)
        M=np.zeros((rows,cols))
        I=np.zeros(dim)  
        for r in range(rows):        
            I.reshape(-1)[P[r]]=1        
            it=InvTransform(I).reshape((1,-1))
            I.reshape(-1)[P[r]]=0
            M[r,:]=it
        self.M_dct=M
        return M,SM.reshape(dim)
    
    
    def binary_measurement_matrix(self,M_dct=None,A=None):
        if A is None:
            A=self.auxiliary_matrix_a()
        if M_dct is None:
            if self.M_dct is None:
                self.dct_sampling_functions()
            M_dct=self.M_dct
        m=A.shape[1]
        p=A.shape[0]-m
        k_dct=M_dct.shape[0]
        k=m+p+k_dct
        n=M_dct.shape[1]
        M=np.zeros((k,n),dtype=np.single,order='F')
        I=np.random.randint(0,m,n)
        for i in range(m):
            M[:m+p,I==i]=A[:,i].reshape((-1,1))    
        for r in range(k_dct):        
                argsrt=np.argsort(M_dct[r,:].ravel()+1e-7*np.random.rand(n))
                M[m+p+r,argsrt[round(n/2 * ( 2*np.random.rand()+m-1)/m) :] ]=1
        self.Mbin=M
        return M
    
        
    def d_fdri(self,Mbin=None, p=None,dim=None,μ=None, ϵ=None, tol=1e-7):
        if Mbin is None:
            if self.Mbin is None:
                self.binary_measurement_matrix()
            Mbin=self.Mbin
        if p is None:
            p=self.p
        if dim is None:
            dim=self.dim
        if μ is None:
           μ=self.μ
        if ϵ is None:
           ϵ=self.ϵ  
        DIFF=self.differential_operator(p)
        M=DIFF(Mbin.reshape(-1,np.prod(dim)))
        Ny,Nx=dim
        w=lambda N:(2*np.pi/N)*np.hstack((np.arange(N//2),np.arange(-N//2,0)))   
        (wx,wy)=np.meshgrid(w(Nx),w(Ny))
        D=1/np.sqrt((1-μ)**2 * (np.sin(wx)**2+np.sin(wy)**2) +ϵ +μ**2*(wx**2+wy**2)/(2*np.pi**2))
        row_fft2=lambda X: fft.fftn(X.reshape((-1,Ny,Nx)),axes=(-2,-1)).reshape((-1,Ny*Nx))
        row_ifft2=lambda X: fft.ifftn(X.reshape((-1,Ny,Nx)),axes=(-2,-1)).reshape((-1,Ny*Nx))
        col_fft2=lambda X: fft.fftn(X.T.reshape((-1,Ny,Nx)),axes=(-2,-1)).reshape((-1,Ny*Nx)).T
        col_ifft2=lambda X: fft.ifftn(X.T.reshape((-1,Ny,Nx)),axes=(-2,-1)).reshape((-1,Ny*Nx)).T 
        FILT_R=lambda X: row_fft2(row_ifft2(X)*D.reshape(-1)); # F*D*F'*X
        FILT_L=lambda X: col_fft2(D.reshape(-1,1)*col_ifft2(X)); # F*D*F'*X
        a=FILT_R(M.reshape((-1,Nx*Ny))).real
        # use svd to calculate the pseudoinverse             
        U, S, V = linalg.svd(a,full_matrices=False); #a = U*S*V'             
        inv_S=1/S
        I=np.abs(S)<tol*np.abs(S).max()
        if np.any(I):
            print('Warning: zero svd values: ',I.sum())
        inv_S[I]=0
        P=np.array(FILT_L((V.T.conj()*inv_S.reshape(-1))@U.T.conj()).real)
        self.Pg=np.array(DIFF(P,right=True),dtype=np.single,order='F')
        return self.Pg
        
    
    def reconstruct(self,P=None,channel=1):
        if P is None:
            P=self.Pg
        def relu(x):
            x[x<0]=0
            return x
        if channel==1:
            P0=np.array(P,dtype=np.single,order='F')
        else:
            P0=np.array(-P,dtype=np.single,order='F')
        return lambda y: relu(P0@y.reshape((-1,1)))
    
    def psnr(self,orig, tstimg): 
        MSE = np.mean((orig - tstimg) ** 2) 
        #vmax = max(orig.max(),tstimg.max())
        vmax=orig.max()
        PSNR = 10 * np.log10(vmax**2 / MSE) 
        return PSNR


def main():
    DFDRI(verbose=True)

if __name__ == "__main__":
    main()

