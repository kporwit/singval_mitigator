import copy
import numpy as np
import sys

#diagonal permutation according to tab
def diag_perm(mat, tab):
    mat_copy = copy.deepcopy(mat)
    for ind in range(0, len(tab)):
        nind = tab[ind]
        mat[ind, ind] = mat_copy[nind, nind]
    return mat

#lower triangular matrix creation and comparison
def create_alpha(singval_tab, eigenval_tab, verbosity):
    if verbosity > 0:
        print 'Choosen s: ' + str(singval_tab)
        print 'Choosen e: ' + str(eigenval_tab)
    n = len(eigenval_tab)
    singval_tab = np.sort(singval_tab)
    if verbosity > 1:
        print '-2*10^-16*n*max(singval_tab) =' + str(-2**(-16)*n*max(singval_tab))
    a = np.diag(singval_tab)
    nzeros = n - sum(np.sign(singval_tab))
    if verbosity > 1:
        print 'nzeros: ' + str(nzeros)
    sum_correction = 0
    if verbosity > 1:
        print 'sum_correction: ' + str(sum_correction)
    temp_e = [1 if x1!=x2 else 0 for (x1, x2) in zip(eigenval_tab, map(abs,
                                                                       eigenval_tab))] #e~=abs(e)
    temp_s = [1 if x1!=x2 else 0 for (x1, x2) in zip(singval_tab, map(abs,
                                                                      singval_tab))] #s~=abs(s)
    e_or_s = [1 if x1!=0 or x2!=0 else 0 for (x1, x2) in zip(temp_e, temp_s)] #or(e~=abs(e), s~=abs(s))
    if max(e_or_s)==1:
        print 'sorry e s must be nonneg'
        sys.exit()
    temp = min(np.log(np.cumprod(np.sort(np.absolute(eigenval_tab))) / np.cumprod(singval_tab)))
    if verbosity > 1:
        print 'temp: ' + str(temp)
    if temp < -2e-16*n*max(singval_tab):
        if verbosity > 0:
            print 'sorry majorization not satsified'
            print 'temp comparison: ' + str(-2e-16*n*max(singval_tab))
        b = np.ones(3)
        return b
    for i in range(0, n - 1):
        posn = i
        if verbosity > 1:
            print 'posn, loop i: ' + str(posn)
        for j in range(i, n - 1):
            if a[j, j] <= eigenval_tab[i]:
                posn = j
            if verbosity > 1:
                print 'posn, loop j: ', str(posn)
        j = posn
        if verbosity > 1:
            print "j = posn: " + str(j)
        correction = max([a[j, j]-eigenval_tab[i], 0])
        if verbosity > 1:
            print 'correction: ' + str(correction)
        if correction > 0:
            if verbosity > 0:
                print 'j'
                print str([j, a[j, j], eigenval_tab[i]])
                print str(a)
            b = np.ones(3)
            return b
        a[j, j] = min([a[j, j], eigenval_tab[i]])
        if verbosity > 1:
            print 'a:'
            print str(a)
        sum_correction = correction + sum_correction
        if verbosity > 1:
            print 'sum_correction: ' + str(sum_correction)
        correction = max([(eigenval_tab[i]-a[j+1, j+1]), 0])
        if verbosity > 1:
            print 'correction: ' + str(correction)
        if correction > 0:
            if verbosity > 0:
                print 'j+1'
                print str([a[j+1, j+1], eigenval_tab[i]])
            b = np.ones(3)
            return b
        a[j+1, j+1] = max(a[j+1,j+1], eigenval_tab[i])
        if verbosity > 1:
            print 'a:'
            print str(a)
        sum_correction = correction + sum_correction
        if verbosity > 1:
            print 'sum_correction: ' + str(sum_correction)
        perm = list(range(0, i))
        perm.append(posn)
        perm.append(posn+1)
        if verbosity > 1:
            print 'perm: ' + str(perm)
        for j in range(i, n):
            if (j != posn) and (j != posn+1):
                perm.append(j)
                if verbosity > 1:
                    print 'appended perm: ' + str(perm)
        if verbosity > 1:
            print 'perm: ' + str(perm)
        diag_perm(a, perm)
        if verbosity > 1:
            print('a:')
            print str(a)
        if eigenval_tab[i] != 0:
            ee = eigenval_tab[i]
            s1 = a[i+1, i+1]
            s2 = a[i, i]
            U = np.eye(2)
            if verbosity > 1:
                print 'ee: ' + str(ee)
                print 's1: ' + str(s1)
                print 's2: ' + str(s2)
                print 'U:'
                print str(U)
            if (s2**2 - s1**2) != 0:
                cu = np.sqrt((s1**2-ee**2)/(s1**2-s2**2))
                su = np.sqrt((ee**2-s2**2)/(s1**2-s2**2))
                U = [[-cu, su], [su, cu]]
                if verbosity > 1:
                    print 'cu: ' + str(cu)
                    print 'su: ' + str(su)
                    print 'U:'
                    print str(U)
            a[i:i+2, :] = np.matmul(U,a[i:i+2, :])
            if verbosity > 1:
                print 'a:'
                print str(a)
            if abs(a[i, i+1]) > abs(a[i, i]):
                a[:, [i,i+1]] = a[:, [i+1,i]]
                if verbosity > 1:
                    print 'a:'
                    print str(a)
            tangent = a[i, i+1]/a[i,i]
            cosine = 1/np.sqrt(1+tangent**2)
            sine = cosine*tangent
            V = [[cosine, sine], [-sine, cosine]]
            if verbosity > 1:
                print 'tangent: ' + str(tangent)
                print 'cosine: ' + str(cosine)
                print 'sine: ' + str(sine)
                print 'V:'
                print str(V)
            a[:, i:i+2] = np.matmul(a[:, i:i+2], np.transpose(V))
            if verbosity > 1:
                print 'a:'
                print str(a)
            a[i, i+1] = 0
            if verbosity > 1:
                print 'a:'
                print str(a)
            a[:,i] = np.sign(a[i, i])*a[:, i]
            if verbosity > 1:
                print 'a:'
                print str(a)
            a[i+1, i+1] = abs(a[i+1, i+1])
            if verbosity > 1:
                print 'a:'
                print str(a)
        else:
            if nzeros == 1:
                if prod(sign(eigenval_tab[i+1:n])) == 0:
                    V = [[0, 1], [1, 0]]
                    if verbosity > 1:
                        print 'V:'
                        print str(V)
                    a[:, i:i+1] = a[:, i:i+1]*V
                    if verbosity > 1:
                        print 'a:'
                        print str(a)
                    a[i+1, i+1] = abs(a[i+1, i+1])
                    if verbosity > 1:
                        print 'a:'
                        print str(a)
                    nzeros = nzeros+1
                    if verbosity > 1:
                        print 'nzeros: ' + str(nzeros)
                else:
                    f = prod(eigenval_tab[i+1:n])/prod(diag(a[i+2:n, i+2:n]))
                    a[i+1, i] = sqrt(a[i+1, i+1]**2 - f**2)
                    if verbosity > 1:
                        print 'f:' + str(f)
                        print 'a:'
                        print str(a)
                    a[i+1, i+1] = f
                    if verbosity > 1:
                        print 'a:'
                        print str(a)
            nzeros = nzeros - 1
            if verbosity > 1:
                print 'nzeros: ' + str(nzeros)
        indx = np.argsort(np.diag(a[i+1:n, i+1:n])).tolist()
        perm = range(0,i+1) + [x+i+1 for x in indx]
        a = diag_perm(a, perm)
        if verbosity > 0:
            print 'indx: ' + str(indx)
            print 'perm: ' + str(perm)
            print 'last a:'
            print str(a)
            print('==================================================')
    if sum_correction > 2e-16*n*max(singval_tab):
        if verbosity > 0:
            print 'warning'
            print str(sum_correction)
        b = np.ones(3)
        return b
    return a
