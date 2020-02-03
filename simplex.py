import numpy as np

## max 50x1+100x2
## x1+x2<=300
## 2x1+x2<=400
## x2<=250

M=10000

##系数矩阵-标准化后
A=np.mat([[1,1,1,0,0],[2,1,0,1,0],[0,1,0,0,1]])

b=np.mat([[300],[400],[250]])

c=np.mat([[50],[100],[0],[0],[0]])

row=A.shape[0]
column=A.shape[1]

xb=[2,3,4]
xn=[0,1]

xb.sort()
xn.sort()

def calc_sigma(_BI):
    ctn=c[xn].T
    ctb=c[xb].T
    N=A[:,xn]
    sigma=ctn-ctb.dot(_BI).dot(N)
    return sigma

def calc_theta(_BI,_xin):
    bi=_BI.dot(b)
    ai=_BI.dot(A[:,_xin])
    theta=[]
    for i in range(0,bi.shape[0]):
        if(ai[i]!=0):
            theta.append(float(bi[i]/ai[i]))
        else:
            theta.append(M)
    theta=np.array(theta)
    return theta

def calc_obj(_xb,_xn):
    sol=[]
    for i in range(0,column):
        sol.append(0)
    _B=A[:,_xb]
    _BI=_B.I
    _bi=_BI.dot(b)
    index=0
    for i in _xb:
        sol[i]=float(_bi[index])
        index=index+1
    sol=np.mat(sol)
    obj=sol.dot(c)

    return obj


if __name__=='__main__':
    stop=0
    while(stop==0):
        print(xb,"xb",xn,"xn")
        B=A[:,xb]
        print(B,"B")
        BI=B.I
        sigma=calc_sigma(BI)

        print(sigma,"检验数")
    
        print(BI.dot(b),"b")

        #存在使目标函数优化的非基变量
        if((sigma>0).any()==True):
            stop=0
        else:#最优
            stop=1
            continue
        
        #入基变量
        xin=xn[np.argmax(sigma)]

        #计算出基变量的theta值
        theta=calc_theta(BI,xin)
        print(theta,"theta")
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
        print(xin,"xin",xout,"xout")
        print("------------------------------------")
    print(calc_obj(xb,xn),"obj")

    


