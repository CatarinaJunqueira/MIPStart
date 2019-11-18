# -*- coding: utf-8 -*-
"""
Created on Thu May  9 16:12:39 2019

@author: catar
"""

from __future__ import print_function
from scipy.io import loadmat
from connect_ifttt import email_alert
import numpy as np
import cplex
import socket

def Write_MST_Gurobi(casename):
    
    global debug
    debug = False
    data = loadmat(casename)
    S1 = data['Seq_Retirada']
    C = data['C'][0][0]
    R = data['R'][0][0]
    Seq_Navio_Inv = data['Seq_Navio_Inv'].tolist()[0]
    Seq_Navio_Id_Inv = data['Seq_Navio_Id_Inv'].tolist()[0]
    Patios = data['patio'].tolist()
    q_o = data['q_o'].tolist()
    q_d = data['q_d'].tolist()
    q_r = data['q_r'].tolist()
    q_c = data['q_c'].tolist()
    w_o = data['w_o'].tolist()
    w_d = data['w_d'].tolist()
    w_a = data['w_a'].tolist()
    w_r = data['w_r'].tolist()
    w_c = data['w_c'].tolist()
    phi = data['phi'].tolist()
    Npatios = len(Patios)
    P=Npatios+1 # numero de portos
    
    for o in range(Npatios):
        for d in range(P):
            if phi[o][d].shape[1] != 0 :
                phi[o][d] = phi[o][d].tolist()[0]         
            else:
                phi[o][d] = []
    
    omega=[ [] for i in range(Npatios) ] # omega = conjunto dos indices dos conteineres em cada patio
    S = []
    for i in range(Npatios):
        Patios[i]=Patios[i][0]
        omega[i]=np.extract(Patios[i]!= 0 , Patios[i]).tolist()
        S.append(S1[0][i].tolist()[0])            
        
    N=[ 0 for i in range(Npatios) ] # N = quantidade de conteineres em cada patio
    for i in range(Npatios):
        N[i]=np.count_nonzero(Patios[i])
    
    T=N
    
    H=[] # H = numero de linhas de cada patio
    for i in range(Npatios):
        H.append(Patios[i].shape[0])
    
    W= [] # W = numero de colunas de cada patio
    for i in range(Npatios):
        W.append(Patios[i].shape[1])
    
    print('parametros criados')
    
    model = cplex.Cplex()
    start_time = model.get_time()
    model.objective.set_sense(model.objective.sense.minimize)    
    startVar=[]
    startVal=[]    
        
    #------------------------------------------------------------#
    #--------------------  Variaveis  ---------------------------#
    #------------------------------------------------------------#
    nvar = 0 
    model,nvar,startVar,startVal = variavel_v(model,S,N,T,nvar,omega,startVar,startVal)
    model,nvar,startVar,startVal = variavel_q(model,N,R,C,nvar,q_o,q_d,q_r,q_c,startVar,startVal)
    model,nvar,startVar,startVal = variavel_u(model,N,R,C,nvar,Seq_Navio_Inv,startVar,startVal)
    model,nvar,startVar,startVal = variavel_w(model,N,R,C,nvar,w_o,w_d,w_a,w_r,w_c,startVar,startVal)
    model,nvar,startVar,startVal = variavel_z(model,omega,N,T,R,C,S,Seq_Navio_Id_Inv,nvar,startVar,startVal)
    model,nvar,startVar,startVal = variavel_y(model,omega,Patios,S,N,H,W,T,nvar,startVar,startVal)
    model,nvar,startVar,startVal = variavel_b(model,omega,Patios,S,N,H,W,T,nvar,startVar,startVal)
    model,nvar,startVar,startVal = variavel_x(model,omega,N,H,W,T,nvar,startVar,startVal)
    print('variaveis criadas')
    
    solucao_inicial_gurobi_mst = casename + '.mst'    
    out_file = open(solucao_inicial_gurobi_mst,'w+') 
    out_file.write("# MIP start \n")   

    for Var,Val in zip(startVar,startVal):

        out_file.write(str(Var)+" "+str(int(Val)) + "\n")
    
    out_file.close()

            
#------------------------------------------------------------#
#--------------------  Variaveis  ---------------------------#
#------------------------------------------------------------#
def variavel_v(model,S,N,T,nvar,omega,startVar,startVal):
    global vind
    vind = dict()
    indx = nvar
    vnames = []
    val =[]
    for o in range(len(N)):  
        for n in omega[o]:
            flag = False
            for t in range(1,T[o]+1):
              #  if debug:
                vnames.append('v_'+str(n)+'_'+str(t))
                vind[(n, t)] = indx
                indx += 1
                val.append(0.0)            
                if S[o][t-1] == n or flag == True :
                    if flag == False:
                        flag = True
                        continue
                    else:                           
                        val[-1]= 1.0
                        flag = True 
                        
                        
    nv = len(vnames)
    #nv= indx-nvar
    nvar += nv
    lb = [0.0]*nv
    ub = [1.0]*nv
    ctypes =[model.variables.type.binary]*nv
    model.variables.add(obj=[], lb=lb, ub=ub, types=ctypes,names=vnames)                    
    #startVar+=list(range(nvar-nv,nvar))
    startVar+=vnames
    startVal+=val            
  #  model.MIP_starts.add(cplex.SparsePair(ind = vnames, val = val), model.MIP_starts.effort_level.repair, "first")               
 
    return model,nvar,startVar,startVal

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------# 
    
def variavel_q(model,N,R,C,nvar,q_o,q_d,q_r,q_c,startVar,startVal):
    
    qqnames = []
    global qind
    qind = dict()
    indx = nvar
    Y=np.size(q_o)
    
    for o in range(Y):
        qqnames.append('q_'+str(q_o[0][o])+'_'+str(q_d[0][o])+'_'+str(q_r[0][o])+'_'+str(q_c[0][o]))
    
    qnames = []
    val =[]
    for o in range(1,len(N)+1):  
        for d in range(o+1,len(N)+2):
              for r in range(1,R+1):
                  for c in range(1,C+1):
                      qind[(o,d, r, c)] = indx
                      indx += 1
                      if 'q_'+str(o)+'_'+str(d)+'_'+str(r)+'_'+str(c) in qqnames :
                          qnames.append('q_'+str(o)+'_'+str(d)+'_'+str(r)+'_'+str(c))
                          val.append(1.0)
                      else:
                          qnames.append('q_'+str(o)+'_'+str(d)+'_'+str(r)+'_'+str(c))
                          val.append(0.0)
    

    nq = len(qnames)   
    nvar += nq   
    lb = [0.0]*nq
    ub = [1.0]*nq
    ctypes =[model.variables.type.binary]*nq
    model.variables.add(obj=[], lb=lb, ub=ub, types=ctypes,names=qnames)
    #startVar+=list(range(nvar-nq,nvar)) 
    startVar+=qnames
    startVal+=val                       
   # model.MIP_starts.add(cplex.SparsePair(ind = qnames, val = val), model.MIP_starts.effort_level.repair, "first")    
    
    return model,nvar,startVar,startVal

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------# 

def variavel_u(model,N,R,C,nvar,Seq_Navio_Inv,startVar,startVal):
    
    unames = []
    val =[]
    global uind
    uind = dict()
    indx = nvar
    for o in range(1,len(N)+1):  
        for r in range(1,R+1):
            for c in range(1,C+1):
                unames.append('u_'+str(o)+'_'+str(r)+'_'+str(c))
                uind[(o, r, c)] = indx
                indx += 1
                if Seq_Navio_Inv[o-1][r-1,c-1] == 0:
                    val.append(0.0)
                else:
                    val.append(1.0)
    
    
    nu = len(unames)
    nvar += nu
    lb = [0.0]*nu
    ub = [1.0]*nu
    ctypes =[model.variables.type.binary]*nu
    model.variables.add(obj=[], lb=lb, ub=ub, types=ctypes,names=unames)
    #startVar+=list(range(nvar-nu,nvar))
    startVar+=unames
    startVal+=val                     
    #model.MIP_starts.add(cplex.SparsePair(ind = unames, val = val), model.MIP_starts.effort_level.repair, "first")     
    
    return model,nvar,startVar,startVal

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------# 

def variavel_w(model,N,R,C,nvar,w_o,w_d,w_a,w_r,w_c,startVar,startVal):
    wnames = []
    global wind
    wind = dict()
    indx = nvar                        
    obj = []
    val=[]
    wwnames = []
    Y=np.size(w_o)
    
    for o in range(Y):
        wwnames.append('w_'+str(w_o[0][o])+'_'+str(w_d[0][o])+'_'+str(w_a[0][o])+'_'+str(w_r[0][o])+'_'+str(w_c[0][o]))

    for o in range(1,len(N)+1):  
        for d in range(o+1,len(N)+2):
            for a in range(o+1,d+1):
                for r in range(1,R+1):
                    for c in range(1,C+1):
                        wind[(o, d,a, r, c)] = indx
                        indx += 1                        
                        if 'w_'+str(o)+'_'+str(d)+'_'+str(a)+'_'+str(r)+'_'+str(c) in wwnames :
                          wnames.append('w_'+str(o)+'_'+str(d)+'_'+str(a)+'_'+str(r)+'_'+str(c))
                          val.append(1.0)
                        else:
                          wnames.append('w_'+str(o)+'_'+str(d)+'_'+str(a)+'_'+str(r)+'_'+str(c))
                          val.append(0.0)
                        
                        if a == d:
                            obj.append(0.0)
                        else:
                            obj.append(1.0)
    
    nw = len(wnames)    
    nvar += nw    
    lb = [0.0]*nw
    ub = [1.0]*nw
    ctypes =[model.variables.type.binary]*nw
    model.variables.add(obj=obj, lb=lb, ub=ub, types=ctypes,names=wnames)
    #startVar+=list(range(nvar-nw,nvar))
    startVar+=wnames
    startVal+=val      
   # model.MIP_starts.add(cplex.SparsePair(ind = wnames, val = val), model.MIP_starts.effort_level.repair, "first") 
    
    return model,nvar,startVar,startVal

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#   
def variavel_z(model,omega,N,T,R,C,S,Seq_Navio_Id_Inv,nvar,startVar,startVal):
    znames = []
    val=[]
    global zind
    zind = dict()
    indx = nvar
    for o in range(len(N)):  
        for n in omega[o]:
            flag = 0
            for t in range(1,T[o]+1):
                for r in range(1,R+1):
                    for c in range(1,C+1):
                        znames.append('z_'+str(n)+'_'+str(t)+'_'+str(r)+'_'+str(c))
                        zind[( n, t,r,c)] = indx
                        indx += 1
                    #    if t==16 and r==1 and c==1:
                    #        print('z_'+str(n)+'_'+str(t)+'_'+str(r)+'_'+str(c))
                        if Seq_Navio_Id_Inv[o][r-1,c-1] == n and S[o][t-1] == n:                         
                           val.append(1.0)
                           flag=1
                           row=r
                           col=c
                           continue
                        else:
                           val.append(0.0) 
                        if flag == 1 and row == r and col == c:
                           val[-1]= 1.0
                           
    nz = len(znames)   
    nvar += nz   
    lb = [0.0]*nz
    ub = [1.0]*nz
    ctypes =[model.variables.type.binary]*nz
    model.variables.add(obj=[], lb=lb, ub=ub, types=ctypes,names=znames)
    #startVar+=list(range(nvar-nz,nvar))
    startVar+=znames
    startVal+=val 
   # model.MIP_starts.add(cplex.SparsePair(ind = znames, val = val), model.MIP_starts.effort_level.repair, "first")   
    
    return model,nvar,startVar,startVal

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#   
def variavel_y(model,omega,Patios,S,N,H,W,T,nvar,startVar,startVal):
    yynames = []
    val=[]
    global yind
    yind = dict()

    indx = nvar
    
    for o in range(len(N)): #para cada patio
        for t in range(1,T[o]): #para cada tempo que sai um conteiner n
            contador = H[o]+1
            for j in range(1,H[o]+1):  #percorrer as linhas do patio
                contador = contador - 1
                for i in range(1,W[o]+1): #percorrer as colunas  do patio                 
                    if Patios[o][j-1,i-1]==S[o][t-1]:
                       n=S[o][t-1]
                       yynames.append('y_'+str(i)+'_'+str(contador)+'_'+str(n)+'_'+str(t))
                       
                       continue    
    
    ynames = []
    for o in range(len(N)):  
        for i in range(1,W[o]+1): 
            for j in range(1,H[o]+1):
                for n in omega[o]:
                    for t in range(1,T[o]+1):
                        yind[(i,j,n, t)] = indx
                        indx += 1
                        if 'y_'+str(i)+'_'+str(j)+'_'+str(n)+'_'+str(t) in yynames:
                            ynames.append('y_'+str(i)+'_'+str(j)+'_'+str(n)+'_'+str(t))
                            val.append(1.0)
                        else:
                            ynames.append('y_'+str(i)+'_'+str(j)+'_'+str(n)+'_'+str(t))
                            val.append(0.0)
    
    ny = len(ynames)    
    nvar += ny    
    lb = [0.0]*ny
    ub = [1.0]*ny
    ctypes =[model.variables.type.binary]*ny
    model.variables.add(obj=[], lb=lb, ub=ub, types=ctypes,names=ynames)
    #startVar+=list(range(nvar-ny,nvar))
    startVar+=ynames
    startVal+=val 
   # model.MIP_starts.add(cplex.SparsePair(ind = ynames, val = val), model.MIP_starts.effort_level.repair, "first")   
    
    return model,nvar,startVar,startVal

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------# 
def variavel_b(model,omega,Patios,S,N,H,W,T,nvar,startVar,startVal):
    bnames = []
    val=[]
    global bind
    bind = dict()
    indx = nvar
    flag=0
    for o in range(len(N)):  
        for i in range(1,W[o]+1):
            contador = H[o]+1
            for j in range(1,H[o]+1):
                contador = contador - 1
                for n in omega[o]:
                    flag = 1
                    for t in range(1,T[o]+1):
                        bnames.append('b_'+str(i)+'_'+str(contador)+'_'+str(n)+'_'+str(t))
                        bind[(i, contador, n, t)] = indx
                        indx += 1
                        val.append(0.0) 
                        if Patios[o][j-1,i-1]==n and flag == 1: # se o conteiner n eh o que esta na posicao (i,j), entao b_ijnt=1
                           val[-1]= 1.0 
                           if Patios[o][j-1,i-1]==S[o][t-1]: # se chegamos no tempo t em que o conteiner n sai, nos proximos tempos b_ijnt deve ser igual a 0 
                              flag = 0

    nb = len(bnames)    
    nvar += nb    
    lb = [0.0]*nb
    ub = [1.0]*nb
    ctypes =[model.variables.type.binary]*nb
    model.variables.add(obj=[], lb=lb, ub=ub, types=ctypes,names=bnames) 
   # startVar+=list(range(nvar-nb,nvar))
    startVar+=bnames
    startVal+=val
   # model.MIP_starts.add(cplex.SparsePair(ind = bnames, val = val), model.MIP_starts.effort_level.repair, "first")   
    
    return model,nvar,startVar,startVal
#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------# 
def variavel_x(model,omega,N,H,W,T,nvar,startVar,startVal):
    xnames = []
    val=[]
    global xind
    xind = dict()
    
    indx = nvar
    for o in range(len(N)):  
        for i in range(1,W[o]+1): 
            for j in range(1,H[o]+1):
                for k in range(1,W[o]+1):
                    for l in range(1,H[o]+1):
                        for n in omega[o]:
                            for t in range(1,T[o]+1):
                                xnames.append('x_'+str(i)+'_'+str(j)+'_'+str(k)+'_'+str(l)+'_'+str(n)+'_'+str(t))
                                val.append(0.0)
                                xind[(i,j,k,l,n,t)] = indx
                                indx +=1
    
    nx = len(xnames)    
    nvar += nx    
    lb = [0.0]*nx
    ub = [1.0]*nx
    ctypes =[model.variables.type.binary]*nx
    obj = [1.0]*nx
    model.variables.add(obj=obj, lb=lb, ub=ub, types=ctypes,names=xnames)    
    #startVar+=list(range(nvar-nx,nvar))
    startVar+=xnames
    startVal+=val
    #model.MIP_starts.add(cplex.SparsePair(ind = xnames, val = val), model.MIP_starts.effort_level.repair, "first")   
    
    return model,nvar,startVar,startVal
#------------------------------------------------------------------------------------------------------------------#
# FIM
#------------------------------------------------------------------------------------------------------------------#

    
if __name__ == "__main__":
#    i=1
#    for n in range(2):
#        if n==3 or n==4 or n==7 or n==8 or n==10 or n==11 or n==12:
#           name_instance = 'SolucaoInicial_Aleatoria_Instancia_' + str(i) 
#           print(name_instance)
#           _ =  Write_MST_Gurobi(name_instance)
#           i = i+1
#        else:
#           name_instance = 'SolucaoInicial_Instancia_' + str(i) 
#           print(name_instance)
#           _ =  Write_MST_Gurobi(name_instance)
#           i = i+1    
     _ =  Write_MST_Gurobi('SolucaoInicial_InstanciaIntegradaTeste3.mat')
#     i=1
#     for n in range(12):
#         if n==3 or n==4 or n==7 or n==8 or n==10 or n==11 or n==12:
#             name_instance = 'SolucaoInicial_Aleatoria_Instancia_' + str(i) 
#             print(name_instance)
#             _ =  Write_MST_Gurobi(name_instance)
#             i = i+1
#         else:
#             name_instance = 'SolucaoInicial_Instancia_' + str(i) 
#             print(name_instance)
#             _ =  Write_MST_Gurobi(name_instance)
#             i = i+1
#             
#     i=1
#     for t in range(12):
#         if t==7:
#             name_instance = 'SolucaoInicial_Aleatoria_Instancia_' + str(i) + '_1' 
#             print(name_instance)
#             _ =  Write_MST_Gurobi(name_instance)
#             i = i+1
#         else:
#             name_instance = 'SolucaoInicial_Instancia_' + str(i) 
#             print(name_instance)
#             _ =  Write_MST_Gurobi(name_instance)
#             i = i+1             