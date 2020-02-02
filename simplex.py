import numpy as np

## max 5x1+6x2+7x4
## x1+x2+2x3<=1
## x1+2x2+x3<=2

####确定数据
##系数矩阵-标准化后
A=[[1,1,1,0],[1,2,0,1]]
A=np.mat(A)

b=[[1],[2]]
b=np.mat(b)

c=[[5],[6],[0],[0]]
c=np.mat(c)

m=A.shape[0]

xb=[2,3]
xn=[0,1]

xb.sort()
xn.sort()
max_a=max(xb)
print(max_a)

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
    theta=np.array(theta)
    return theta

stop=0
while(stop==0):
    print(stop)
    B=A[:,xb]
    BI=B.I
    sigma=calc_sigma(BI)

    print(sigma)
    print(xb,"xb",xn,"xn")
    print(BI.dot(b),"b")

    #存在使目标函数优化的非基变量
    if((sigma>0).any()==True):
        stop=0
    else:
        stop=1
        continue
        
    #入基变量
    xin=xn[np.argmax(sigma)]

    #计算出基变量的theta值
    theta=calc_theta(BI,xin)
    if((theta>0).any()==True):
        stop=0
    else:
        stop=1
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



    
##    print(xout,"xout")
##
##        
##    
##
##    print(theta,"theta")
##    print(np.argmax(sigma))
##
##    print(xin)
##    
##    print((sigma>0).all())
##    stop=1




##C_A = c[[0,2]]    #先取出想要的行数据
##C_A = C_A[:,[2,3]] #再取出要求的列数据



B=[[0.5,0,0],[0.5,1,0],[0,0,1]]

B=np.mat(B)

c=B[:,[1,2]]
##c=np.mat(c)

b=[[1],[0],[0]]

##print(B.I.dot(b))
print(c)
##print(A.I)
##print(A)
##print(B,B.I)



