# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 09:57:16 2019

@author: catar
"""
from __future__ import print_function
from scipy.io import loadmat
from connect_ifttt import email_alert
import numpy as np
import cplex
import socket

def SolucaoInicial(casename):
    
    global debug
    debug = False
    data = loadmat(casename)
    S1 = data['Seq_Retirada']
    C = data['C'][0][0]
    R = data['R'][0][0]
 #  Seq_Navio = data['Seq_Navio'].tolist()
    Seq_Navio_Inv = data['Seq_Navio_Inv'].tolist()[0]
    Seq_Navio_Id_Inv = data['Seq_Navio_Id_Inv'].tolist()[0]
  # Seq_Patio = data['Seq_Patio'].tolist()[0]
  # N_Mov=data['value'][0][0]
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
    TT =data['TT']
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
   # model,nvar,startVar,startVal = variavel_x(model,omega,N,H,W,T,nvar,startVar,startVal)
    print('variaveis criadas')
    
    #------------------------------------------------------------#
    #-------------------  Restricoes  ---------------------------#
    model = restricao_P0(model,omega,N,H,W,Patios)
    model = restricao_P1(model,omega,N,H,W,T)
    model = restricao_P2(model,omega,N,H,W,T)
    model = restricao_P3(model,omega,N,H,W,T)
    model = restricao_P6(model,omega,N,H,W,T)
    model = restricao_P7(model,omega,N,H,W,T)
   # model = restricao_P8(model,omega,N,H,W,T)
    model = restricao_P9(model,omega,N,H,W,T)
  #  model = restricao_P10(model,omega,N,H,W,T)
    model = restricao_PA(model,omega,N,H,W,T)
    print('restricoes do patio criadas')
    model = restricao_I1(model,omega,N,T)
    model = restricao_I2(model,omega,N,T)
    model = restricao_I3(model,omega,N,R,C,T)
    model = restricao_I4(model,omega,N,R,C,T,P) #problem
    model = restricao_I5(model,omega,N,R,C,T)
    model = restricao_I6(model,phi,R,C,T,P,N,startVar,startVal) #problem
    model = restricao_I7(model,omega,N,R,C,T)
    model = restricao_I8(model,omega,N,R,C,T,P) #problem
    model = restricao_I9(model,omega,N,R,C,P)
    model = restricao_I10(model,omega,N,R,C,TT)
    print('restricoes de integracao criadas')
    model = restricao_N1(model,P,R,C,TT)
    model = restricao_N2(model,R,C,P)
    model = restricao_N3(model,N,R,C)
    model = restricao_N4(model,P,R,C)
    print('restricoes criadas')
    
    #model.write(casename+".lp")
    
    variaveis = model.variables.get_num()
    print('Numero de Variaveis = ',variaveis)    
    restricoes = model.linear_constraints.get_num()
    print('Numero de Restricoes = ',restricoes)
  #  variaveis2 = len(startVar)
  #  print('Tamanho startVar = ', variaveis2)
  #  restricoes2 = len(startVal)
 #   print('Tamanho startVal = ', restricoes2)
    
    z = 'No modelo ha %s variaveis e %s restricoes' %(variaveis,restricoes)
    
    modelotxt = casename + '.txt'    
    out_file = open(modelotxt,'w+')  
    model.set_results_stream(out_file)    
    model.set_results_stream(modelotxt)
    
    model.MIP_starts.add(cplex.SparsePair(ind = startVar, val = startVal), model.MIP_starts.effort_level.repair, "first") 
    
   # model.MIP_starts.write("SolucaoInicial_Instancia_1.mst")
    model.parameters.timelimit.set(432000)  # limite de tempo em segundos (5 dias)
    model.parameters.threads.set(20)
    
    print('resolvendo o modelo')
    model.solve()
    
    out_file.seek(0)
    out_string =out_file.read()
    out_file.close()
    
    print(out_string)
   
    end_time = model.get_time()
    solvetime = end_time-start_time
    print('Duracao = ',solvetime)
    
    status = model.solution.get_status_string()
    fobj = model.solution.get_objective_value()
    Navio = obterNavio(model,R,C,Npatios)
    Pt = obterP(model,H[0],W[0],omega[0],T[0])
    
    print('\nSolution status = ',status)
    print('Valor da Funcao Objetivo: ',fobj )
    
    Y = 'Instancia: %s <br> Rodado em: %s <br> Solution status: %s  <br> Valor da Funcao Objetivo: %s  <br> Duracao: %s <br>' %(casename,socket.gethostname(),status,fobj,solvetime)
   
    email_alert(Y,z, out_string)
    
    excelfilename = 'Resultado_das_Instancias_MIP_Start'
    
    t1=0
    
  #  save_log_excel(model,casename,t1,solvetime,excelfilename)
            
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
#------------------------------------------------------------------------------------------------------------------# 

#------------------------------------------------------------------------------------------------------------------#
#-------------------  Restricoes  ---------------------------#
#------------------------------------------------------------------------------------------------------------------#
def restricao_P0(model,omega,N,H,W,Patios):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for i in range(W[o]): 
            for j in range(H[o]):
                for n in omega[o]:
                    if debug:
                        rest_names.append('restP0_'+str(i+1)+'_'+str(j+1)+'_'+str(n)+'_'+str(o+1))

                    if Patios[o][H[o]-j-1,i] == n :
                        rhs.append(1.0)
                    else :
                        rhs.append(0.0)

                    #add b coeficientes
                    cols = [bind[(i+1,j+1,n,1)]]
                    coefs = [1.0]
                    expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
#Constraint 1: In each time period, each block must either be within the stack or in the outside region:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P1(model,omega,N,H,W,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for n in omega[o]:
            for t in range(T[o]):
                if debug:
                    rest_names.append('restP1_'+str(n)+'_'+str(t+1)+'_'+str(o+1))
                rhs.append(1.0)
                cols = [vind[(n,t+1)]]
                coefs = [1.0]
                
                for i in range(W[o]):  
                    for j in range(H[o]):
                        #add b coeficientes
                        cols.append(bind[(i+1,j+1,n,t+1)])
                        coefs.append(1.0)
                expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
#Constraint 2: In each time period, each slot (i,j) must be occupied by at most one block:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P2(model,omega,N,H,W,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for i in range(W[o]):  
            for j in range(H[o]):
                for t in range(T[o]):
                    if debug:
                        rest_names.append('restP2_'+str(i+1)+'_'+str(j+1)+'_'+str(t+1)+'_'+str(o+1))
                    rhs.append(1.0)
                    cols = []
                    coefs = []
                    for n in omega[o]:          
                        #add b coeficientes
                        cols.append(bind[(i+1,j+1,n,t+1)])
                        coefs.append(1.0)
                        
                    expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model


#------------------------------------------------------------------------------------------------------------------#
#Constraint 3: garante que não hajam ‘buracos’ no pátio ao restringir que se há um contêiner posição $(i,j+1)$, 
#então a posição $(i,j)$ abaixo também deve estar ocupada:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P3(model,omega,N,H,W,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for i in range(W[o]):  
            for j in range(H[o]-1):
                for t in range(T[o]):
                    if debug:
                        rest_names.append('restP3_'+str(i+1)+'_'+str(j+1)+'_'+str(t+1)+'_'+str(o+1))
                    rhs.append(0.0)
                    cols = []
                    coefs = []
                    for n in omega[o]:          
                        #add b coeficientes
                        cols.append(bind[(i+1,j+2,n,t+1)])
                        coefs.append(1.0)
                        cols.append(bind[(i+1,j+1,n,t+1)])
                        coefs.append(-1.0)
                        
                    expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
#Constraint 4: restrição de equilíbrio de fluxo entre as variáveis de configuração e de movimento no pátio.  
#Vincula o layout no período t com o layout no período t + 1 através das retiradas e realocações executadas:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P6(model,omega,N,H,W,T):

    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):
        for i in range(W[o]):
            for j in range(H[o]):
                for n in omega[o]:
                    for t in range(1,T[o]):
                        if debug:
                            rest_names.append(
                                'restP6_' + str(i + 1) + '_' + str(j + 1) + '_' + str(n) + '_' + str(t + 1) + '_' + str(o + 1))
                        rhs.append(0.0)
                        cols = [bind[(i + 1, j + 1, n, t + 1)]]  # b_{i,j,n,t}
                        coefs = [1.0]
                        cols.append(bind[(i + 1, j + 1, n, t)])  # -b_{i,j,n,t-1}
                        coefs.append(-1.0)
                        cols.append(yind[(i + 1, j + 1, n, t)])  # y_{i,j,n,t-1}
                        coefs.append(1.0)
                        #for k in range(W[o]):
                         #   for l in range(H[o]):
                          #      if k != i or l != j:
                           #         cols.append(xind[(k + 1, l + 1, i + 1, j + 1, n,t)])  # -\sum_{k=1}^{W}\sum_{l=1}^{H}x_{k,l,i,j,n,t-1}
                            #        coefs.append(-1.0)
                             #       cols.append(xind[(i + 1, j + 1, k + 1, l + 1, n, t)])  # \sum_{k=1}^{W}\sum_{l=1}^{H}x_{i,j,k,l,n,t-1}
                              #      coefs.append(1.0)

                        expr.append(cplex.SparsePair(cols, coefs))


    senses=["E"]*len(expr)

    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)
    return model




#------------------------------------------------------------------------------------------------------------------#
#Constraint 5:  define a variável $v_{nt}$ e assegura que todos os contêineres sejam retirados do pátio:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P7(model,omega,N,H,W,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for n in omega[o]:
            for t in range(1,T[o]):
                if debug:
                    rest_names.append('restP7_'+str(n)+'_'+str(t+1)+'_'+str(o+1))
                rhs.append(0.0)
                cols = [vind[(n,t+1)]] #v_{nt}
                coefs = [1.0]                       
                for i in range(W[o]):          
                    for j in range(H[o]):
                        for tt in range(t): 
                            cols.append(yind[(i+1,j+1,n,tt+1)]) #-\sum_{i=1}^{W}\sum_{j=1}^{H}\sum_{t'=1}^{t-1}y_{ijnt'}
                            coefs.append(-1.0)
                        
                expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 6:  garante a política LIFO, ou seja, se no período t, o contêiner $n$ está abaixo do contêiner $q$  
# e o contêiner $n$ é remanejado, então no período $t + 1$ o contêiner $n$ não pode estar alocado em uma posição 
# abaixo do contêiner $q$:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P8(model,omega,N,H,W,T):
        
    rhs = []
    rest_names = []
    expr = []   

    for o in range(len(N)):
        M = N[o]*((H[o]-1)**2)

        for i in range(W[o]): 
            for k in range(W[o]):
                for j in range(H[o]-1):
                    for l in range(H[o]-1):
                        for t in range(T[o]):
                            if debug:
                                rest_names.append('restP8_'+str(i+1)+'_'+str(j+1)+'_'+str(k+1)+'_'+str(l+1)+'_'+str(t+1)+'_'+str(o+1))
                            rhs.append(M)
                            cols = [] 
                            coefs = []
                            
                            for n in omega[o]: 
                                cols.append(xind[(i+1,j+1,k+1,l+1,n,t+1)]) #\sum_{n=1}^{N}x_{i,j,k,l,n,t}
                                coefs.append(M)
                            
                                for jj in range(j+1,H[o]):
                                    for ll in range(l+1,H[o]):
                                        cols.append(xind[(i+1,jj+1,k+1,ll+1,n,t+1)])
                                        coefs.append(1.0)               
                        
                            expr.append(cplex.SparsePair(cols,coefs))   
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 7:  garante que sejam remanejados apenas os contêineres que estão acima, ou seja, na mesma coluna, 
# de um contêiner a ser retirado:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P9(model,omega,N,H,W,T):

        
    rhs = []
    rest_names = []
    expr = []
 
        
    for o in range(len(N)):
        M = (H[o] ** 2) * W[o]*(W[o]-1)
        for i in range(W[o]):
            for t in range(T[o]):
                for n in omega[o]:
                    if debug:
                       rest_names.append('restP9_'+str(n)+'_'+str(i+1)+'_'+str(t+1)+'_'+str(o+1))
                    rhs.append(M)
                    cols = [] 
                    coefs = []
                    
                    for j in range(H[o]):
                        cols.append(bind[(i+1,j+1,n,t+1)])
                        coefs.append(M)
                        #for k in range(W[o]):
                         #   for l in range(H[o]):
                          #      for ii in range(i):
                           #         cols.append(xind[(ii+1,j+1,k+1,l+1,n,t+1)])
                            #        coefs.append(1.0)
                             #   for iii in range(i+1,W[o]):
                              #      cols.append(xind[(iii+1,j+1,k+1,l+1,n,t+1)])
                               #     coefs.append(1.0)          
                        
                    expr.append(cplex.SparsePair(cols,coefs))   
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 8:  garante que nenhum contêiner pode ser remanejado para outra posição que esteja da mesma coluna na
# qual ele se encontra:
#------------------------------------------------------------------------------------------------------------------#
def restricao_P10(model,omega,N,H,W,T):
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for i in range(W[o]): 
            for j in range(H[o]):
                for l in range(H[o]):
                    for n in omega[o]: 
                        for t in range(T[o]):
                            if debug:
                                rest_names.append('restP10_'+str(i+1)+'_'+str(j+1)+'_'+str(l+1)+'_'+str(n)+'_'+str(t+1)+'_'+str(o+1))
                            rhs.append(0.0)
                            cols= [xind[(i+1,j+1,i+1,l+1,n,t+1)]]
                            coefs = [1.0]        
                        
                            expr.append(cplex.SparsePair(cols,coefs))   
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 9:  garante que um contêiner na posição $(i,j)$ só pode ser movido depois que o contêiner na posição 
# $(i,j+1)$ é movido. Se o contêiner na posição $(i,j+1)$ não é movido então temos que $b_{i(j+1)nt} = 1$ e 
# $x_{i(j+1)klnt} = 0$, e o lado esquerdo da equação se torna 0. Consequentemente o lado direito da equação também 
# deve ser 0. Dessa forma, nenhuma realocação ou remanejamento é permitido para o contêiner na posição $(i,j)$:
#------------------------------------------------------------------------------------------------------------------#
def restricao_PA(model,omega,N,H,W,T):
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for i in range(W[o]): 
            for j in range(H[o]-1):
                for t in range(T[o]):
                    if debug:
                        rest_names.append('restA_'+str(i+1)+'_'+str(j+1)+'_'+str(t+1)+'_'+str(o+1))
                    rhs.append(1.0)
                    cols = [] 
                    coefs = []
                    for n in omega[o]:
                        cols.append(bind[(i+1,j+2,n,t+1)]) #\sum_{n=1}^{N}b_{i,j+1,n,t}
                        coefs.append(1.0)
                        cols.append(yind[(i+1,j+1,n,t+1)]) #\sum_{n=1}^{N}y_{i,j,n,t}
                        coefs.append(1.0)
#                        for k in range(W[o]):
#                            for l in range(H[o]):
#                                cols.append(xind[(i+1,j+2,k+1,l+1,n,t+1)]) #\sum_{k=1}^{W}\sum_{l=1}^{H}\sum_{n=1}^{N}x_{i,j+1,k,l,n,t}
#                                coefs.append(-1.0)  
#                                cols.append(xind[(i+1,j+1,k+1,l+1,n,t+1)]) #\sum_{k=1}^{W}\sum_{l=1}^{H}\sum_{n=1}^{N}x_{i,j,k,l,n,t}
#                                coefs.append(1.0)

                    expr.append(cplex.SparsePair(cols,coefs))  

    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model


#------------------------------------------------------------------------------------------------------------------#
# Constraint 10: garante que em cada período de tempo um contêiner seja retirado do pátio:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I1(model,omega,N,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)): 
        for t in range(T[o]):
            if debug:
                rest_names.append('restI1_'+str(t+1)+'_'+str(o+1))
            rhs.append(t)
            cols = []
            coefs = []
            for n in omega[o]:          
                #add v coeficientes
                cols.append(vind[(n,t+1)])
                coefs.append(1.0)
                        
            expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model



#------------------------------------------------------------------------------------------------------------------#
# Constraint 11: define a variável $v_{nt}$. Quando um contêiner $n$ é retirado do pátio, a variável $v_{nt}$ se 
# torna 1 e se mantém igual a 1 nos períodos de tempo seguintes:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I2(model,omega,N,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):
        for n in omega[o]: 
            for t in range(T[o]-1):
                if debug:
                    rest_names.append('restI2_'+str(n)+'_'+str(t+1)+'_'+str(o+1))
                rhs.append(0.0)     
                    #add v coeficientes
                cols = [vind[(n,t+1)]]
                coefs = [1.0]
                cols.append(vind[(n,t+2)])
                coefs.append(-1.0)
                        
                expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model
#------------------------------------------------------------------------------------------------------------------#
# Constraint 12: garante que o contêiner $n$ seja carregado no navio no período de tempo $t$:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I3(model,omega,N,R,C,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):
        for n in omega[o]: 
            for t in range(T[o]-1):
                if debug:
                    rest_names.append('restI3_'+str(n)+'_'+str(t+1)+'_'+str(o+1))
                rhs.append(0.0)     
                    #add v coeficientes
                cols = [vind[(n,t+2)]]
                coefs = [-1.0]
                for r in range(R):
                    for c in range(C):
                        cols.append(zind[(n,t+1,r+1,c+1)])
                        coefs.append(1.0)
                    
                expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 13: assegura que uma posição $(r,c)$ no navio só pode ser ocupada por um contêiner, seja ele um 
# contêiner que foi carregado no porto atual (porto $o$), em algum porto anterior (porto $o-1$) ou um contêiner
# que já estava no navio e está sendo remanejado em $o$:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I4(model,omega,N,R,C,T,P):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(1,len(N)):
        for r in range(R):
            for c in range(C):
                for t in range(T[o]):
                    if debug:
                        rest_names.append('restI4_'+str(r+1)+'_'+str(c+1)+'_'+str(t+1)+'_'+str(o+1))
                    rhs.append(1.0)
                    cols = [] 
                    coefs = []
                    
                    for oo in range(o):
                        for d in range(o+1,P):
                            for a in range(o+1,d+1):
                                cols.append(wind[(oo+1,d+1,a+1,r+1,c+1)])
                                coefs.append(1.0)
                    
                    for d in range(o+1,P):
                        cols.append(qind[(o+1,d+1,r+1,c+1)])
                        coefs.append(1.0)                            
                    
                    for n in omega[o]: 
                        cols.append(zind[(n,t+1,r+1,c+1)])
                        coefs.append(1.0)
                    
                    expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 14: certifica que o contêiner $n$, depois de carregado, não mude de posição enquanto o navio estiver
# parado no mesmo porto:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I5(model,omega,N,R,C,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):
        for n in omega[o]: 
            for t in range(T[o]-1):
                for r in range(R):
                    for c in range(C):
                        if debug:
                            rest_names.append('restI5_'+str(n)+'_'+str(t+1)+'_'+str(r+1)+'_'+str(c+1)+'_'+str(o+1))
                        rhs.append(0.0)     
                       #add z coeficientes
                        cols=[zind[(n,t+1,r+1,c+1)]]
                        coefs = [1.0]
                        cols.append(zind[(n,t+2,r+1,c+1)])
                        coefs.append(-1.0)
                        
                        expr.append(cplex.SparsePair(cols,coefs))
    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 15:  garante que se há um contêiner na posição $(r,c)$ do navio, ele deve ser um contêiner que acabou 
# de ser embarcado, ou um contêiner de remanejamento:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I6(model,phi,R,C,T,P,N,startVar,startVal):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):    
          for r in range(R):
              for c in range(C):
                  for d in range(o+1,P):
                      rest_names.append('restI6_'+str(r+1)+'_'+str(c+1)+'_'+str(o+1)+'_'+str(d+1))
                      rhs.append(0.0)     
                      #add q coeficientes
                      cols = ['q_'+str(o+1)+'_'+str(d+1)+'_'+str(r+1)+'_'+str(c+1)]
                      coefs = [1.0]
                      
                      index = startVar.index('q_'+str(o+1)+'_'+str(d+1)+'_'+str(r+1)+'_'+str(c+1))
                      lhs = -startVal[index]
                      
                      for n in phi[o][d]:
                          cols.append('z_'+str(n)+'_'+str(T[o])+'_'+str(r+1)+'_'+str(c+1))
                          coefs.append(1.0)
                          
                          index = startVar.index('z_'+str(n)+'_'+str(T[o])+'_'+str(r+1)+'_'+str(c+1))
                          lhs -= startVal[index] 
                          
                      for a in range(o+1,d+1):
                          cols.append('w_'+str(o+1)+'_'+str(d+1)+'_'+str(a+1)+'_'+str(r+1)+'_'+str(c+1))
                          coefs.append(-1.0)                   
                          
                          index = startVar.index('w_'+str(o+1)+'_'+str(d+1)+'_'+str(a+1)+'_'+str(r+1)+'_'+str(c+1))
                          lhs += startVal[index] 
                          
                      if abs(lhs) > 1e-6 :
                          print(rest_names[-1])
                          
                          index = startVar.index('q_'+str(o+1)+'_'+str(d+1)+'_'+str(r+1)+'_'+str(c+1))
                          print('q_'+str(o+1)+'_'+str(d+1)+'_'+str(r+1)+'_'+str(c+1)+' = '+str(startVal[index]))
                          
                          for n in phi[o][d]:
                          
                              index = startVar.index('z_'+str(n)+'_'+str(T[o])+'_'+str(r+1)+'_'+str(c+1))
                              print( 'z_'+str(n)+'_'+str(T[o])+'_'+str(r+1)+'_'+str(c+1)+' = '+str(startVal[index])) 
                          
                          for a in range(o+1,d+1):                                      
                              index = startVar.index('w_'+str(o+1)+'_'+str(d+1)+'_'+str(a+1)+'_'+str(r+1)+'_'+str(c+1))
                              print( 'w_'+str(o+1)+'_'+str(d+1)+'_'+str(a+1)+'_'+str(r+1)+'_'+str(c+1)+' = '+str(startVal[index])) 
                          
                          
                          
                      expr.append(cplex.SparsePair(cols,coefs))   
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 16: assegura que todos os $N_{o}$ contêineres do pátio $o$ já foram embarcados no navio:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I7(model,omega,N,R,C,T):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):    
          for n in omega[o]:
              for t in range(T[o]-1):
                  if debug:
                      rest_names.append('restI7_'+str(n)+'_'+str(t+1)+'_'+str(o+1))
                  rhs.append(1.0)
                  cols = [] 
                  coefs = []
                  for r in range(R):
                      for c in range(C):
                          cols.append(zind[(n,T[o],r+1,c+1)])
                          coefs.append(1.0)                                     
                       
                  expr.append(cplex.SparsePair(cols,coefs))   
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 17:garante que, durante o processo de carregamento do navio, nenhum contêiner seja alocado em uma 
# posição flutuante ou que ocupe a posição de um contêiner que já estava no navio ou foi remanejado:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I8(model,omega,N,R,C,T,P):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):
        for r in range(R-1):
            for c in range(C):
                for t in range(T[o]):
                    if debug:
                        rest_names.append('restI8_'+str(r+1)+'_'+str(c+1)+'_'+str(t+1)+'_'+str(o+1))
                    rhs.append(0.0)
                    cols = [] 
                    coefs = []
                    for n in omega[o]: 
                        cols.append(zind[(n,t+1,r+1,c+1)])
                        coefs.append(-1.0)
                        cols.append(zind[(n,t+1,r+2,c+1)])
                        coefs.append(1.0)                      
                    
                    for oo in range(o):
                        for d in range(o+1,P):
                            for a in range(o+1,d+1):
                                cols.append(wind[(oo+1,d+1,a+1,r+1,c+1)])
                                coefs.append(-1.0)
                    
                    for d in range(o+1,P):
                        cols.append(qind[(o+1,d+1,r+1,c+1)])
                        coefs.append(-1.0)                            
                    
                    expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 18: contabiliza o número total de contêineres que foram remanejados no porto $o$:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I9(model,omega,N,R,C,P):
    
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):
        for d in range(o+1,P):
            if debug:
                rest_names.append('restI9_'+str(o+1)+'_'+str(d+1))
            rhs.append(0.0)
            cols = [] 
            coefs = []
            for oo in range(o):
                for r in range(R):
                    for c in range(C):
                        cols.append(wind[(oo+1,d+1,o+1,r+1,c+1)])
                        coefs.append(1.0)
            for r in range(R):
                for c in range(C):
                    cols.append(qind[(o+1,d+1,r+1,c+1)])
                    coefs.append(-1.0)
                            
                    
            expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model

#------------------------------------------------------------------------------------------------------------------#
# Constraint 19: mantém a estabilidade do navio:
#------------------------------------------------------------------------------------------------------------------#
def restricao_I10(model, omega, N, R, C, TT):
    rhs = []
    rest_names = []
    expr = []

    v = np.sum(TT, axis=0)
    theta = [0.0] * len(N)
    theta[0] = N[0] - v[0]

    for i in range(1, len(N)):
        theta[i] = theta[i - 1] + N[i] - v[i]

    for o in range(len(N)):
        temp = int(np.ceil(float(theta[o]) / float(C)))

        if temp < R:
            if debug:
                rest_names.append('restI10_' + str(o + 1))
            rhs.append(0.0)
            cols = []
            coefs = []

            for c in range(C):
                for r in range(temp, R):
                    cols.append('u_' + str(o + 1) + '_' + str(r + 1) + '_' + str(c + 1))
                    coefs.append(1.0)

            expr.append(cplex.SparsePair(cols, coefs))

    senses = ["E"] * len(expr)

    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs, names=rest_names)
    return model
    
#------------------------------------------------------------------------------------------------------------------#
# Constraint 20: restrição de conservação de fluxo, e indica que o número total de contêineres no porto $o$ deve ser
# igual ao número de contêineres que foram embarcados nos portos $p=1,2,...,o$ menos os contêineres que foram 
# desembarcados nos portos $p=1,2,...,o$ :
#------------------------------------------------------------------------------------------------------------------#
def restricao_N1(model,P,R,C,TT):   
    rhs = []
    rest_names = []
    expr = []
    for o in range(P-1):
        for d in range(o+1,P):
            if debug:
                rest_names.append('restN1_'+str(o+1)+'_'+str(d+1))
            rhs.append(float(TT[o,d]))
            cols = [] 
            coefs = []
            for a in range(o+1,d+1):        
                for r in range(R):
                    for c in range(C):
                        cols.append(wind[(o+1,d+1,a+1,r+1,c+1)])
                        coefs.append(1.0)
            
            for m in range(o):
                for r in range(R):
                    for c in range(C):            
                        cols.append(wind[(m+1,d+1,o+1,r+1,c+1)])
                        coefs.append(-1.0)                                
                    
            expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model  
    
#------------------------------------------------------------------------------------------------------------------#
# Constraint 21: garante que cada posição $(r, c)$ tenha no máximo um único contêiner:
#------------------------------------------------------------------------------------------------------------------#
def restricao_N2(model,R,C,P):   
    rhs = []
    rest_names = []
    expr = []
    for o in range(P-1):
        for r in range(R):
            for c in range(C):
                if debug:
                    rest_names.append('restN2_'+str(o+1)+'_'+str(r+1)+'_'+str(c+1))
                rhs.append(0.0)
                cols = [uind[(o+1,r+1,c+1)]]
                coefs = [-1.0]
                for m in range(o+1):
                    for d in range(o+1,P):
                        for a in range(o+1,d+1):
                            cols.append(wind[(m+1,d+1,a+1,r+1,c+1)])
                            coefs.append(1.0)                                                
                        
                expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["E"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model      
    
#------------------------------------------------------------------------------------------------------------------#
# Constraint 22: garante que existem contêineres embaixo do contêiner que ocupa a célula $(r, c)$:
#------------------------------------------------------------------------------------------------------------------#
def restricao_N3(model,N,R,C):   
    rhs = []
    rest_names = []
    expr = []
    for o in range(len(N)):
        for r in range(R-1):
            for c in range(C):
                if debug:
                    rest_names.append('restN3_'+str(o+1)+'_'+str(r+1)+'_'+str(c+1))
                rhs.append(0.0)
                cols = [uind[(o+1,r+1,c+1)]]
                coefs = [-1.0]
                cols.append(uind[(o+1,r+2,c+1)])
                coefs.append(1.0)
                                              
                        
                expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model      
    
#------------------------------------------------------------------------------------------------------------------#
# Constraint 23: é responsável por definir como um contêiner pode ser desembarcado no porto $d$ ao impor que se
# um contêiner que ocupa a posição $(r, c)$, então ele será desembarcado no porto $d$, se não houver um contêiner
# na posição $(r+1, c)$ acima dele:
#------------------------------------------------------------------------------------------------------------------#
def restricao_N4(model,P,R,C):   
    rhs = []
    rest_names = []
    expr = []
    
    for d in range(1,P):
        for r in range(R-1):
            for c in range(C):
                if debug:
                    rest_names.append('restN4_'+str(d+1)+'_'+str(r+1)+'_'+str(c+1))
                rhs.append(1.0)
                cols = [] 
                coefs = []
                for o in range(d):
                    for e in range(d,P):
                        cols.append(wind[(o+1,e+1,d+1,r+1,c+1)])
                        coefs.append(1.0)
                    for u in range(d+1,P):
                        for a in range(d+1,u+1):
                            cols.append(wind[(o+1,u+1,a+1,r+2,c+1)])
                            coefs.append(1.0)           
                        
                expr.append(cplex.SparsePair(cols,coefs))    
    
    senses=["L"]*len(expr)
                    
    model.linear_constraints.add(lin_expr=expr, senses=senses, rhs=rhs,names=rest_names)               
    return model 
 
#------------------------------------------------------------------------------------------------------------------#
# Obter Navio: Retorna como foi a configuração do navio em cada porto 
#------------------------------------------------------------------------------------------------------------------#

def obterNavio(model,R,C,Npatios):

    Navio = []
    for o in range(Npatios):
        Navio.append(np.zeros((R,C)))
        for r in range(R):
            for c in range(C):
                if model.solution.get_values('u_'+str(o+1)+'_'+str(r+1)+'_'+str(c+1)) == 1 :
                    for k in range(o+1):
                        for j in range(o+1,Npatios+1):
                            for v in range(o+1,j+1):
                                if model.solution.get_values('w_'+str(k+1)+'_'+str(j+1)+'_'+str(v+1)+'_'+str(r+1)+'_'+str(c+1)) == 1 :
                                    Navio[o][r,c]=j+1
                                 
    return Navio                                 
#------------------------------------------------------------------------------------------------------------------#
# Obter Patio: Retorna como foi a movimentacao dos patios em cada porto 
#------------------------------------------------------------------------------------------------------------------#
                              
def obterP(model,H,W,omega,T):

    Pt = []                  
 #   for t in range(T):    
 #       Pt.append(np.zeros((H,W)))
 #       for i in range(W):
 #           for j in range(H):
 #               for n in omega:
 #                  if model.solution.get_values('b_'+str(i+1)+'_'+str(j+1)+'_'+str(n)+'_'+str(t+1)) == 1 :                                                     
 #                      Pt[t][H-j-1,i] == n 
 #   return Pt


    for t in range(1,T+1):
        Pt.append(np.zeros((H,W)))
        for i in range(1,W+1):
            for j in range(1,H+1):
                for n in omega:               
                    if model.solution.get_values('b_'+str(i)+'_'+str(j)+'_'+str(n)+'_'+str(t)) == 1 :
                       Pt[t-1][H-j,i-1] = n 
    return Pt
                    
#------------------------------------------------------------------------------------------------------------------#
# FIM
#------------------------------------------------------------------------------------------------------------------#

    
if __name__ == "__main__":
     i=1
     for n in range(12):
         if n==3 or n==4 or n==7 or n==8 or n==10 or n==11 or n==12:
             name_instance = 'SolucaoInicial_Aleatoria_Instancia_' + str(i) 
             print(name_instance)
             _ =  SolucaoInicial(name_instance)
             i = i+1
         else:
             name_instance = 'SolucaoInicial_Instancia_' + str(i) 
             print(name_instance)
             _ =  SolucaoInicial(name_instance)
             i = i+1
          
   #  _ = SolucaoInicial('SolucaoInicial_Instancia_1.mat')