##from pymprog import *
import numpy as np
import time
import math

M=10000

## 纸卷每个长16，现需要25个3m长，20个6m长，18个7m长的纸卷，如何切割时用料最少
## model 1
## min x1+x2+x3
## 5x1>=25
## 2x2>=20
## 2x3>=15

g_A=np.mat([[5,0,0,1,-1,0,0,0,0],
          [0,2,0,0,0,1,-1,0,0],
          [0,0,2,0,0,0,0,1,-1]])

g_b=np.mat([[25],[20],[18]])

g_c=np.mat([[1],[1],[1],[M],[M],[M],[M],[M],[M]])

g_xb=[3,5,7]
g_xn=[0,1,2,4,6,8]

g_xb.sort()
g_xn.sort()

g_sp=[]

def calc_sigma(_A,_BI,_c,_xn,_xb):
    ctn=_c[_xn].T
    ctb=_c[_xb].T
    N=_A[:,_xn]
    sigma=ctn-ctb.dot(_BI).dot(N)
    return sigma

def calc_theta(_BI,_xin,_b,_A):
    bi=_BI.dot(_b)
    ai=_BI.dot(_A[:,_xin])
    theta=[]
    for i in range(0,bi.shape[0]):
        if(ai[i]>0):
            theta.append(float(bi[i]/ai[i]))
        else:
            theta.append(M)
    theta=np.array(theta)
    return theta

def calc_obj(_xb,_xn,_A,_b,_c):
    sol=[]
    for i in range(0,_A.shape[1]):
        sol.append(0)
    _B=_A[:,_xb]
    _BI=_B.I
    _bi=_BI.dot(_b)
    index=0
    for i in _xb:
        sol[i]=float(_bi[index])
        index=index+1
    sol=np.mat(sol)
    print(sol,"sol")
    obj=sol.dot(_c)
    return obj

def calc_shadow_price(_xb,_A,_c):
    _B=_A[:,_xb]
    _BI=_B.I
    _sp=_c[_xb].T.dot(_BI)
##    print(_sp,"影子价格")
    return _sp

def simplex_main(A,xb,xn,b,c):
    stop=0
    while(stop==0):
##        print(xb,"xb",xn,"xn")
        B=A[:,xb]
        print(B,"B")
        BI=B.I
        sigma=calc_sigma(A,BI,c,xn,xb)

##        print(sigma,"检验数")
    
##        print(BI.dot(b),"b")

        #存在使目标函数优化的非基变量
        ##min =>  "<" ;  max=>  ">"
        if((sigma<0).any()==True):
            stop=0
        else:#最优
            stop=1
            continue
        
        #入基变量
        ##min =>  "argmin" ;  max=>  "argmax"
        xin=xn[np.argmin(sigma)]

        #计算出基变量的theta值
        theta=calc_theta(BI,xin,b,A)
##        print(theta,"theta")
        if((theta>0).any()==True):
            stop=0
        else:#无界
            stop=2
            continue
        #出基变量
        xout=xb[np.argmin(theta)]

        #更新基可行解
        xb.remove(xout)
        xb.append(xin)
        xn.remove(xin)
        xn.append(xout)
        xb.sort()
        xn.sort()
##        print(xin,"xin",xout,"xout")
        print("------------------------------------")
    if(stop==1):
        print(calc_obj(xb,xn,A,b,c),"obj")

    elif(stop==2):
        print("模型无界")


def dp_for_subproblem(_sp):
    print(_sp.item(1),"shou")
    max_item=[6,3,3]
    _obj=[]
    _index=[]
    ##初始化数组
    for i in range(0,3):
        _1=[]
        _0=[]
        for j in range(0,17):
            _1.append(1)
            _0.append(0)
        _obj.append(_1)
        _index.append(_0)
    print(_obj[2][15])
    ##计算3m长对应状态
    for i in range(0,17):
        _obj[0][i]=1-_sp.item(0)*math.floor(i/3)

    ##计算6m长对应状态
    for i in range(0,17):
        index=-1
        for j in range(0,i+1):
            obj=_obj[0][j]-_sp.item(1)*math.floor((i-j)/6)
            if (obj<_obj[1][i]):
                _obj[1][i]=obj
                index=j
        _index[1][i]=index

    ##计算7m长对应状态
    index=-1
    for i in range(0,17):
       obj=_obj[1][i]-_sp.item(2)*math.floor((16-i)/7)
       if(obj<_obj[2][16]):
           _obj[2][16]=obj
           index=i
    _index[2][16]=index
    
    if(_obj[2][16]>=-0.0000000001):
        return 0,0
    else:
        sol3=math.floor((16-_index[2][16])/7)
        sol2=math.floor((_index[2][16]-_index[1][_index[2][16]])/6)
        sol1=math.floor(_index[1][_index[2][16]]/3)
        sol=np.mat([[sol1],[sol2],[sol3]])
        return 1,sol
    
if __name__=='__main__':


    
    start_time=time.time()
    algorithm_stop=0
    while(algorithm_stop==0):
        simplex_main(g_A,g_xb,g_xn,g_b,g_c)
        g_sp=calc_shadow_price(g_xb,g_A,g_c)
        print("-------------------simplex completed!!!!!!!!!!!!!!!!!")
        best,sub_sol=dp_for_subproblem(g_sp)
        if(best==0):
            algorithm_stop=1
            break
        else:
            g_xn.append(g_A.shape[1])
            g_A=np.c_[g_A,sub_sol]
            ro=np.mat([[1]])
            g_c=np.r_[g_c,ro]
   
    
##    simplex_main(g_A,g_xb,g_xn,g_b,g_c)

##    print(g_A[:,[9,10,11]])
    end_time=time.time()
    print("Execution Time: ", end_time - start_time)
##    print(g_xb,"g_xb",g_xn,"g_xn")
##    g_sp=calc_shadow_price(g_xb,g_A,g_c)
    
    

    


