import numpy as np
from wyniki import *
#only arrays are storred in the file so using '*' is justified

#experimental bounds
ones=np.eye(3)

exp1=np.array([[1.3e-3, 0, 0],\
               [6.8e-4, 2.2e-4, 0],\
               [2.7e-3, 1.2e-3, 2.8e-3]],\
              dtype=np.float64)
exp2=np.array([[2.4e-2, 0, 0],\
               [2.5e-2, 2.2e-2, 0],\
               [6.9e-2, 1.2e-2, 1.0e-1]],\
              dtype=np.float64)
exp3=np.array([[1.0e-2, 0, 0],\
               [1.7e-2, 1.4e-2, 0],\
               [4.5e-2, 5.3e-2, 1.0e-1]],\
              dtype=np.float64)
exp=np.array([exp1, exp2, exp3])

exp_mi1=np.triu(np.subtract(ones, exp1))
exp_mi2=np.triu(np.subtract(ones, exp2))
exp_mi3=np.triu(np.subtract(ones, exp3))
exp_mi=np.array([exp_mi1, exp_mi2, exp_mi3])

exp_ma1=np.add(exp_mi1, exp1)
exp_ma2=np.add(exp_mi2, exp2)
exp_ma3=np.add(exp_mi3, exp3)
exp_ma=np.array([exp_ma1, exp_ma2, exp_ma3])

#s_min and s_max: row-scenarios; columns-mass
s_mi=np.array([[s3_mi1, s3_mi2, s3_mi3],\
               [s2_mi1, s2_mi2, s2_mi3],\
               [s1_mi1, s1_mi2, s1_mi3]])

s_ma=np.array([[s3_ma1, s3_ma2, s3_ma3],\
               [s2_ma1, s2_ma2, s2_ma3],\
               [s1_ma1, s1_ma2, s1_ma3]])

#tables of indices
row=[0, 1, 2, 1, 2, 2]
col=[0, 1, 2, 0, 0, 1]
#row=[2]
#col=[2]
#tables of indices
srow=[0, 0, 0, 1, 1, 1, 2, 2, 2]
mcol=[0, 1, 2, 0, 1, 2, 0, 1, 2]
#srow=[0]
#mcol=[2]
#stepsize
stepsize=1e-07
for scen, mass in zip(srow,mcol):
    c=[0, 0, 0, 0, 0, 0]
    nc=[0, 0, 0, 0, 0, 0]
    indx=0
    print 'scenario:', scen+1
    print 'mass:', mass+1
    #comp=exp_mi[mass]!=s3_mi[mass]
    compmi=np.subtract(s_mi[scen,mass],exp_mi[mass])
    compma=np.subtract(exp_ma[mass],s_ma[1,mass])
    #print 'compmi:', compmi
    for i, j in zip(row, col):
        print i+1,j+1
        #print compmi[i,j]/stepsize
        for fix in xrange(int(compmi[i,j]/stepsize)):
            step=np.full((3,3),stepsize,dtype=np.float64)
            step=np.tril(step)
            step[i,j]=0
            newa=np.copy(s_mi[scen,mass])
            newa[i,j]=exp_mi[mass][i,j]
            if fix != 0:
                newa[i,j]+=stepsize
            #print 'fix', fix
            while True:
                #print newa
                old_step=np.copy(step)
                step=np.where(newa<s_ma[scen,mass],step,0)
                if np.array_equal(step,old_step)==False:
                    test=np.subtract(old_step,step)
                    newa=np.subtract(newa,test)
                #print step
                if (step==0).all():
                    #print 'break by newa'
                    break
                #if (newa>s_ma[scen,mass]).any():
                #    print 'break by newa'
                #    break
                u,s,vh=np.linalg.svd(newa)
                #print 's:',s
                if (s<1.0).any():
                    c[indx]+=1
                elif (s>=1.0).all():
                    print 'fix:',fix
                    print 'newa:\n',newa
                    print 'singular values:',s
                    nc[indx]+=1
                newa+=step
        indx+=1
    print 'number of contractions for a elements: ', c
    print 'number of non-contractions for a elements: ', nc
