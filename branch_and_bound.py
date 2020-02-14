from pymprog import *
import numpy as np
import time
import math
import copy

M=10000

## 纸卷每个长16，现需要25个3m长，20个6m长，18个7m长的纸卷，如何切割时用料最少
## model 1
## min x1+x2+x3
## 5x1>=25
## 2x2>=20
## 2x3>=15

##g_A=np.mat([[5,0,0,1,-1,0,0,0,0],
##            [0,2,0,0,0,1,-1,0,0],
##            [0,0,2,0,0,0,0,1,-1]])
##
##g_b=np.mat([[25],[20],[18]])
##
##g_c=np.mat([[1],[1],[1],[M],[M],[M],[M],[M],[M]])
##
##g_xb=[3,5,7]
##g_xn=[0,1,2,4,6,8]
##
##g_xb.sort()
##g_xn.sort()
##
##g_sp=[]




## model 2
## max 0.043x1+0.027x2+0.025x3+0.022x4+0.045x5
## x1+x2+x3+x4+x5<=1000
## x2+x3+x4>=400
## 2x1+2x2+x3+x4+5x5<=1.4(x1+x2+x3+x4+x5)
## 9x1+15x2+4x3+3x4+2x5<=5(x1+x2+x3+x4+x5)

g_A=np.mat([[1,1,1,1,1,1,0,0,0,0],
          [0,1,1,1,0,0,1,-1,0,0],
          [0.6,0.6,-0.4,-0.4,3.6,0,0,0,1,0],
          [4,10,-1,-2,-3,0,0,0,0,1]])

g_b=np.mat([[1000],[400],[0],[0] ])

g_c=np.mat([[-0.043],[-0.027],[-0.025],[-0.022],[-0.045],[0],[0],[0],[0],[0]])

g_xb=[5,6,8,9]
g_xn=[0,1,2,3,4,7]

g_xi=[6] 

g_x_check=[0,1,2,3,4]

x_bounds = []
x_bounds.append([0, M])
x_bounds.append([0, M])
x_bounds.append([0, M])
x_bounds.append([0, M])
x_bounds.append([0, M])

g_Node_list=[]

g_stop=0

g_node_counter=0

class Node:
    def __init__(self, x_bounds=[], index=0):
        self._x_bounds = x_bounds
        self._index = index
        print("create Node:", index)
        print('')

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

def generate_matrix(_x_bounds,_g_A,_g_xb,_g_xn,_g_b,_g_c,_g_xi):
    _A=copy.deepcopy(_g_A)
    _b=copy.deepcopy(_g_b)
    _xb=copy.deepcopy(_g_xb)
    _xn=copy.deepcopy(_g_xn)
    _c=copy.deepcopy(_g_c)
    _xi=copy.deepcopy(_g_xi)
    for i in range(0,5):
        if(_x_bounds[i][0]!=0):
            ro=np.mat(np.zeros(_A.shape[1]))
            ro.itemset(i,1)
            _A=np.r_[_A,ro]
            
            _col=np.mat(np.zeros(_A.shape[0])).T
            _col.itemset((_A.shape[0]-1),1)
            _A=np.c_[_A,_col]
            _col.itemset((_A.shape[0]-1),-1)
            _A=np.c_[_A,_col]
             
            _b=np.r_[_b,np.mat([_x_bounds[i][0]])]
            _xb.append(_A.shape[1]-2)
            _xn.append(_A.shape[1]-1)
            _c=np.r_[_c,np.mat([[M],[0]])]
            _xi.append(_A.shape[1]-2)
            print(_c,"cccccc")
        if(_x_bounds[i][1]!=M):
            ro=np.mat(np.zeros(_A.shape[1]))
            ro.itemset(i,1)
            _A=np.r_[_A,ro]
            
            _col=np.mat(np.zeros(_A.shape[0])).T
            _col.itemset((_A.shape[0]-1),1)
            _A=np.c_[_A,_col]
            
            _b=np.r_[_b,np.mat([_x_bounds[i][1]])]
            _xb.append(_A.shape[1]-1)
            _c=np.r_[_c,np.mat([[0]])]
    
    return _A,_xb,_xn,_b,_c,_xi

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
##    print(sol,"sol")
    obj=sol.dot(_c)
    return obj,sol

def calc_shadow_price(_xb,_A,_c):
    _B=_A[:,_xb]
    _BI=_B.I
    _sp=_c[_xb].T.dot(_BI)
##    print(_sp,"影子价格")
    return _sp

def simplex_main(_g_A,_g_xb,_g_xn,_g_b,_g_c,_g_xi):
    A=copy.deepcopy(_g_A)
    xb=copy.deepcopy(_g_xb)
    xn=copy.deepcopy(_g_xn)
    b=copy.deepcopy(_g_b)
    c=copy.deepcopy(_g_c)
    xi=copy.deepcopy(_g_xi)
    g_stop=0
    while(g_stop==0):
##        print(xb,"xb",xn,"xn")
        B=A[:,xb]
##        print(B,"B")
##        print(c,"c")
        BI=B.I
        sigma=calc_sigma(A,BI,c,xn,xb)

##        print(sigma,"检验数")

##        print(BI,b)
##        print(BI.dot(b),"b")

        #存在使目标函数优化的非基变量
        ##min =>  "<" ;  max=>  ">"
        if((sigma<0).any()==True):
            g_stop=0
        else:#最优
            g_stop=1
            continue
        
        #入基变量
        ##min =>  "argmin" ;  max=>  "argmax"
        xin=xn[np.argmin(sigma)]

        #计算出基变量的theta值
        theta=calc_theta(BI,xin,b,A)
##        print(theta,"theta")
        if((theta>0).any()==True):
            g_stop=0
        else:#无界
            g_stop=2
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

    
    if(g_stop==1):
        obj,sol=calc_obj(xb,xn,A,b,c)
##    print(xb,xi)
##    print(sol,"which sol")
    for i in xb:
        if i in xi:
            if(sol.item(i)>0 and sol.item(i)-sol.item(i+1)>0):
                g_stop=3
                break
##        print(obj,"obj")
##        print(sol,"sol")
    if(g_stop==2):
        print("模型无界")
        obj=[],sol=0
    elif(g_stop==3):
        print("模型无可行解")
        obj,sol=calc_obj(xb,xn,A,b,c)
        print(sol,"sol无可行解")
        print(xb)
    print(obj)
    return g_stop,obj,sol,xb

def check_integer(_g_xb,_g_sol):
    integer=True
    for i in _g_xb:
        if(math.pow(_g_sol.item(i)-round(_g_sol.item(i)),2)>0.000000001 and i in g_x_check):
            integer=False
            break
    return integer

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

def generate_node(_g_xb,_g_sol,_Node,g_node_counter):
    _x_bounds=copy.deepcopy(_Node._x_bounds)
    for i in _g_xb:
        if(math.pow(_g_sol.item(i)-round(_g_sol.item(i)),2)>0.000000001):
            _x=i
            _x_value=_g_sol.item(i)
            break
    
    _x_bounds1=copy.deepcopy(_x_bounds)
    _x_bounds2=copy.deepcopy(_x_bounds)
    _x_bounds1[_x][1]=math.floor(_x_value)
    _x_bounds2[_x][0]=math.ceil(_x_value)

    _node1=Node(_x_bounds1,g_node_counter)
    _node2=Node(_x_bounds2,g_node_counter+1)
    g_node_counter=g_node_counter+2
    return _node1,_node2

if __name__=='__main__':
    node=Node(x_bounds,0)
    g_Node_list.append(node)
    bnb_stop=0
    lower_bound=-M
    upper_bound=M
    lower_bound_list=[-M]
    upper_bound_list=[M]
    
    A,xb,xn,b,c,xi =generate_matrix(x_bounds,g_A,g_xb,g_xn,g_b,g_c,g_xi)

    stop,obj,sol,xb=simplex_main(A,xb,xn,b,c,xi)
    print(sol,"SOL")
    print(obj,"OBJ")
    lower_bound_list.append(obj.item(0))
    print("-------------------simplex completed!!!!!!!!!!!!!!!!!")

    current_node=g_Node_list[0]
    node1,node2=generate_node(xb,sol,current_node,g_node_counter)

    g_Node_list.append(node1)
    g_Node_list.append(node2)

    g_Node_list.pop(0)
    
    while(len(g_Node_list)>0):
        current_node=g_Node_list[0]
        print(current_node._index)
        print(current_node._x_bounds)
        
        print(lower_bound_list)
        print(g_xn,"g xn")
        A,xb,xn,b,c,xi =generate_matrix(current_node._x_bounds,g_A,g_xb,g_xn,g_b,g_c,g_xi)
        print(xn)
##        print(A)
        print(xi)
        stop,obj,sol,xb=simplex_main(A,xb,xn,b,c,xi)
        print(sol,"sol")
        print(xb,"xb")

        if(stop!=1):
            g_Node_list.pop(0)
            continue
##        break
        if(check_integer(xb,sol)==False):##非整数解
            if(obj>upper_bound_list[-1]):
                g_Node_list.pop(0)
            else:
                node1,node2=generate_node(xb,sol,current_node,g_node_counter)
                g_Node_list.append(node2)
                g_Node_list.append(node1)
                g_Node_list.pop(0)
                if(obj>lower_bound_list[-1]):
                    lower_bound_list.append(obj.item(0))
        else:##整数解
            g_Node_list.pop(0)
            if(obj<upper_bound_list[-1]):
                upper_bound_list.append(obj.item(0))

        for i in g_Node_list:
            print(i._index,i._x_bounds)
        
        print(lower_bound_list,"LB")
        print("----")
        print(upper_bound_list,"UB")
    

    start_time=time.time()
##    simplex_main(g_A,g_xb,g_xn,g_b,g_c)
##    print(g_xb)
    
##    algorithm_stop=0
##    while(algorithm_stop==0):
##        simplex_main(g_A,g_xb,g_xn,g_b,g_c)
##        g_sp=calc_shadow_price(g_xb,g_A,g_c)
##        
##        best,sub_sol=dp_for_subproblem(g_sp)
##        if(best==0):
##            algorithm_stop=1
##            break
##        else:
##            g_xn.append(g_A.shape[1])
##            g_A=np.c_[g_A,sub_sol]
##            ro=np.mat([[1]])
##            g_c=np.r_[g_c,ro]
   
    
##    simplex_main(g_A,g_xb,g_xn,g_b,g_c)

##    print(g_A[:,[9,10,11]])
    end_time=time.time()
    print("Execution Time: ", end_time - start_time)
##    print(g_xb,"g_xb",g_xn,"g_xn")
##    g_sp=calc_shadow_price(g_xb,g_A,g_c)
    
    

    


