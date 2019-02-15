from numpy import *
import sys


#permutation function
def diag_perm( mat, tab ):
    for ind in range(0, len(tab)): #diagonal permutation according to perm
        nind = perm[ind]
        mat[ind, ind] = mat[nind, nind]
    return mat

s = [1.0003,0.999,0.9964]
e = [0.9972, 0.9998, 0.9987]
print('Choosen s: ', s)
print('Choosen e: ', e)
n = len(e)
s = sort(s)
print('-2*10^-16*n*max(s) =',-2e-16*n*max(s))
a = diag(s)
nzeros = n - sum(sign(s))
print('nzeros: ', nzeros)
sum_correction = 0
print('sum_correction: ', sum_correction)
temp_e = [1 if x1!=x2 else 0 for (x1, x2) in zip(e, map(abs, e))] #e~=abs(e)
temp_s = [1 if x1!=x2 else 0 for (x1, x2) in zip(s, map(abs, s))] #s~=abs(s)
e_or_s = [1 if x1!=0 or x2!=0 else 0 for (x1, x2) in zip(temp_e, temp_s)] #or(e~=abs(e), s~=abs(s))
if max(e_or_s)==1:
    print('sorry e s must be nonneg')
    sys.exit()
temp = min(log(cumprod(sort(e)) / cumprod(s)))
print('temp: ', temp)
if temp < -2e-16*n*max(s):
    print('sorry majorization not satsified')
    sys.exit()
for i in range(0, n - 1):
    posn = i
    print('posn, petla i: ', posn)
    for j in range(i, n - 1):
        if a[j, j] <= e[i]:
            posn = j
        print('posn, petla j: ', posn)
    j = posn
    print("j = posn: ", j)
    correction = max([a[j, j]-e[i], 0])
    print('correction: ',correction)
    if correction > 0:
        print('j')
        print([j, a[j, j], e[i]])
        print(a)
    a[j, j] = min([a[j, j], e[i]])
    print('a:')
    print(a)
    sum_correction = correction + sum_correction
    print('sum_correction: ', sum_correction)
    correction = max([(e[i]-a[j+1, j+1]), 0])
    print('correction: ',correction)
    if correction > 0:
        print('j+1')
        print([a[j+1, j+1], e[i]])
    a[j+1, j+1] = max(a[j+1,j+1], e[i])
    print('a:')
    print(a)
    sum_correction = correction + sum_correction
    print('sum_correction: ', sum_correction)
    perm = list(range(0, i))
    perm.append(posn)
    perm.append(posn+1)
    print('perm: ', perm)
    for j in range(i, n):
        if (j != posn) and (j != posn+1):
            perm.append(j)
            print('appended perm: ', perm)
    diag_perm( a, perm)
    print('a:')
    print(a)
    if e[i] != 0:
        ee = e[i]
        print('ee: ', ee)
        s1 = a[i+1, i+1]
        print('s1: ', s1)
        s2 = a[i, i]
        print('s2: ', s2)
        U = eye(2)
        print('U:')
        print(U)
        if (s2**2 - s1**2) != 0:
            cu = sqrt((s1**2-ee**2)/(s1**2-s2**2))
            print('cu: ', cu)
            su = sqrt((ee**2-s2**2)/(s1**2-s2**2))
            print('su: ', su)
            U = [[-cu, su], [su, cu]]
            print('U:')
            print(U)
        a[i:i+2, :] = matmul(U,a[i:i+2, :])
        print('a:')
        print(a)
        if abs(a[i, i+1]) > abs(a[i, i]):
            a[:, [i,i+1]] = a[:, [i+1,i]]
            print('a:')
            print(a)
        tangent = a[i, i+1]/a[i,i]
        print('tangent: ', tangent)
        cosine = 1/sqrt(1+tangent**2)
        print('cosine: ', cosine)
        sine = cosine*tangent
        print('sine: ', sine)
        V = [[cosine, sine], [-sine, cosine]]
        print('V:')
        print(V)
        a[:, i:i+2] = matmul(a[:, i:i+2],transpose(V))
        print('a:')
        print(a)
        a[i, i+1] = 0
        print('a:')
        print(a)
        a[:,i] = sign(a[i,i])*a[:,i]
        print('a:')
        print(a)
        a[i+1,i+1] = abs(a[i+1,i+1])
        print('a:')
        print(a)
    else:
        if nzeros == 1:
            if prod(sign(e[i+1:n])) == 0:
                V = [[0, 1], [1, 0]]
                print('V:')
                print(V)
                a[:, i:i+1] = a[:, i:i+1]*V
                print('a:')
                print(a)
                a[i+1, i+1] = abs(a[i+1, i+1])
                print('a:')
                print(a)
                nzeros = nzeros+1
                print('nzeros: ', nzeros)
            else:
                f = prod(e[i+1:n])/prod(diag(a[i+2:n, i+2:n]))
                print('f:', f)
                a[i+1, i] = sqrt(a[i+1, i+1]**2 - f**2)
                print('a:')
                print(a)
                a[i+1, i+1] = f
                print('a:')
                print(a)
        nzeros = nzeros - 1
        print('nzeros: ', nzeros)
    
    indx = argsort(diag(a[i+1:n, i+1:n])).tolist()
    print('indx: ', indx)
    perm = range(0,i+1) + [x+i+1 for x in indx]
    print('perm: ', perm)
    a = diag_perm(a, perm)
    print('last a:')
    print(a)
    print('==================================================')

if sum_correction > 2e-16*n*max(s):
    print('warning')
    print(sum_correction)
