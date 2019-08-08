


def LPTClusterEquivalentByTranslation(LPT10,LPT20,returnTranslated):
    from copy import deepcopy
#perform round operation:
    lpt1=deepcopy(LPT10)
    LPT1=deepcopy(LPT10)
    lpt2=deepcopy(LPT20)
    LPT2=deepcopy(LPT20)
#    print "lp1"
#    print lpt1
#    print "lpt2"
#    print lpt2
#    print "\n"
    npt=len(LPT1)
    trLPT2=LPT2
    isSame=False
    LPT1=roundlpts(LPT1)
    LPT1.sort()
    sortLPT1=LPT1
    LPT2=roundlpts(LPT2)
    LPT2.sort()
    sortLPT2=LPT2
    if len(LPT2)!=npt:
#        print('Clusters length different')
        return isSame
    else:
        i=0
        while i < npt:
            if LPT1[0][1] == LPT2[i][1]:
                trDirect=array(LPT1[0][0])-array(LPT2[i][0])
            #loop over atomic positions in LPT2 and lpt2 (return lpt2)
                for j in range(len(LPT2)):
                    LPT2[j][0]=(array(LPT2[j][0])+trDirect).tolist()
                    lpt2[j][0]=(array(lpt2[j][0])+trDirect).tolist()
            #loop end
                LPT2.sort()
                if sortLPT1==LPT2:
                    if returnTranslated:
                        #isSame=lpt2
                        isSame=roundlpts(lpt2)
                    else:
                        isSame=True
                else:
                    isSame=False
                break
            else:
                i+=1
            if i==npt:
                isSame=False
        return isSame
#test
if False:
    print('LPTCluster...test')
    t1=[[[0.0, -1.0, -1.0], 10], [[0.0, 0.0, 0.0], 25], [[-1.0, -1.0, 0.0], 7]]
#    t2=[[[0.0, -1.0, -1.0], 7], [[1.0, -1.0, -2.0], 10], [[1.0, 0.0, -1.0], 25]]
    t2=[[[1.0, -1.0, -2.0], 10],[[0.0, -1.0, -1.0], 7], [[1.0, 0.0, -1.0], 25]]
#    print(LPTClusterEquivalentByTranslation(t1,t2,False))
    print(LPTClusterEquivalentByTranslation(t1,t2,True))



def relativePosition(origMappedOne2One, final):
    resl=[]
    for i in final:
       # resl.append([origMappedOne2One.allindex(i)+1])
        tmp=(array(allindex(i,origMappedOne2One))+1).tolist()
        resl.append(tmp)
    #print tmp
    for i in range(len(resl)):
        flag=len(resl[i])>1
        resl[i]=resl[i][0]
        if flag:
            for j in range(i+1,len(resl)):
                resl[j]=allremove(resl[i], resl[j])
                #resl[j].remove(resl[i])
    return resl


def FCTrans(npt, DIM, gamma, pi):
    dimTensor=int(math.pow(DIM, npt))
    Gm=List2D(0.0, dimTensor, dimTensor)
    for i in range(dimTensor):
        idx1=list(array(IntegerDigits(i,DIM, npt))+1)
        for j in range(dimTensor):
            idx2=list(array(IntegerDigits(j,DIM, npt))+1)
            Gm[i][j]=1.0
            for k in range(npt):
                tmp1=idx1[k]
                tmp2=idx2[pi[k]-1]
                Gm[i][j]=Gm[i][j]*gamma[tmp1-1][tmp2-1]
    return Gm


def ListFlat(listin):
    tmp=[]
    for i in listin:
        tmp=tmp+i
    return tmp


def roundlpt(lpt):
    std=3 #for low symmetry system
   # std=4 #for high symmetry system
    tmp=list(map(float,lpt[0]))
    return [[round(tmp[0],std),round(tmp[1],std),round(tmp[2],std)],lpt[1]]

def roundlpts(lpts):
    for i in range(len(lpts)):
        lpts[i]=roundlpt(lpts[i])
    return lpts



def allindex(value, qlist):
    indices = []
    idx = -1
    while True:
        try:
            idx = qlist.index(value, idx+1)
            indices.append(idx)
        except ValueError:
            break
    return indices


def allremove(case, listin):
    indlist=allindex(case, listin)
    tmp=[]
    for i in range(len(listin)):
        if not i in indlist:
            tmp.append(listin[i])
    return tmp


def List2D(n,x,y):
    x=int(x)
    y=int(y)
    tmp=[]
    for i in range(x):
        tmp.append([n]*y)
    return tmp


def IntegerDigits(n, b, lg):
    tmp=base10ton(n, b)
    tmp=[int(i) for i in str(tmp)]
    if len(tmp)>lg:
        print('Wrong with base transformation')
        return [0]*lg
    elif len(tmp)==lg:
        return tmp
    else:
        return [0]*(lg-len(tmp))+tmp


def base10ton(num, base):
    """Change ``num'' to given base
    Upto base 36 is supported."""

    converted_string, modstring = "", ""
    currentnum = num
    if not 1 < base < 37:
        raise ValueError("base must be between 2 and 36")
    if not num:
        return '0'
    while currentnum:
        mod = currentnum % base
        currentnum = currentnum // base
        converted_string = chr(48 + mod + 7*(mod > 10)) + converted_string
    return converted_string


