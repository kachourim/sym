import numpy as np
import sympy as sp
import matplotlib.pyplot as pl
import tkinter as tk
from numpy.linalg import inv
from sympy.core.function import _coeff_isneg as v_neg


# ============================================
# Define symmetry operations
def rx(a):
    return np.array([[1, 0, 0], [0 , sp.cos(a), -sp.sin(a)], [0, sp.sin(a), sp.cos(a)]])

def ry(a):
    return np.array([[sp.cos(a), 0, sp.sin(a)], [0, 1, 0], [-sp.sin(a), 0, sp.cos(a)]])
    
def rz(a):
    return np.array([[sp.cos(a), -sp.sin(a), 0], [sp.sin(a), sp.cos(a), 0], [0, 0, 1]])

px   = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
py   = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]])
pz   = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]])
iv   = -np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
c2x  = rx(sp.pi)
c2y  = ry(sp.pi)
c2z  = rz(sp.pi)
c4z  = rz(sp.pi/2)
pxy  = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
pxy2 = np.array([[0, -1, 0], [-1, 0, 0], [0, 0, 1]]) 
pxz  = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
pxz2 = np.array([[0, 0, -1], [0, 1, 0], [-1, 0, 0]])
pyz  = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]])
pyz2 = np.array([[1, 0, 0], [0, 0, -1], [0, -1, 0]])
c2xy = -pxy2
c2xz = -pxz2
c2yz = -pyz2

# ============================================
# Define material paramter symmetry and reciprocity conditions

def TSym2(T):
    eq = T - T.transpose()
    S = sp.solve(eq.flatten())
    
    return np.array(sp.Matrix(T).subs(S))
    
def TSym3(T):
    for n in range(3):
        eq = T[n] - T[n].transpose()
        S = sp.solve(eq.flatten())
        T[n] = np.array(sp.Matrix(T[n]).subs(S))
    
    return T
    
def TSym4(T):
    for i in range(3):
        for j in range(3): 
            eq = T[:,:,i,j] - T[:,:,i,j].transpose()
            S = sp.solve(eq.flatten())
            T[:,:,i,j] = np.array(sp.Matrix(T[:,:,i,j]).subs(S))
            
            eq = T[i,j,:,:] - T[i,j,:,:].transpose()
            S = sp.solve(eq.flatten())
            T[i,j,:,:] = np.array(sp.Matrix(T[i,j,:,:]).subs(S))
            
    return T
    
def TRep4(T):               
    for i in range(3):
        for j in range(3):        
            for k in range(3):
                for l in range(3):
                    T[i,j,k,l] = T[k,l,i,j]
                    
    return T

# ============================================
#  Main part

def T2solve(q,T,a):

    if q == 'i':
        q = 'iv'

    Q = eval(q)

    if a:
        b = 1
    else:
        b = sp.Matrix(Q).det()
    
    eq = T - Q.dot(T.dot(b*Q.transpose()))
    S = sp.solve(eq.flatten())
    
    return np.array(sp.Matrix(T).subs(S))
    
def T3solve(q,T,a):

    if q == 'i':
        q = 'iv'

    Q = eval(q)
   
    if a:
        b = 1
    else:
        b = sp.Matrix(Q).det()
    
    M = np.zeros((3,3,3), dtype=object)
    
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        for n in range(3):
        
                            M[i,j,k] = M[i,j,k] + b*Q[i,l]*Q[j,m]*Q[k,n]*T[l,m,n]
                            
    eq = T - M
    S = sp.solve(eq.flatten())
    
    M = np.zeros((3,3,3), dtype=object)
    
    for i in range(3):
       M[i] = np.array(sp.Matrix(T[i]).subs(S))

    return M

def T4solve(q,T,a):

    if q == 'i':
        q = 'iv'
        
    Q = eval(q)

    if a:
        b = 1
    else:
        b = sp.Matrix(Q).det()
    
    M = np.zeros((3,3,3,3), dtype=object)
    
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        for n in range(3):
                            for o in range(3):
                                for p in range(3):
        
                                    M[i,j,k,l] = M[i,j,k,l] + b*Q[i,m]*Q[j,n]*Q[k,o]*Q[l,p]*T[m,n,o,p]
                            
    eq = T - M
    S = sp.solve(eq.flatten())
    
    M = np.zeros((3,3,3,3), dtype=object)
    
    for i in range(3):
        for j in range(3):
            M[i,j,:,:] = np.array(sp.Matrix(T[i,j,:,:]).subs(S))

    return M

def T3Reshape(M):
    return np.array([[M[0, 0, 1], M[0, 0, 2], M[0, 1, 2], M[0, 0, 0], M[0, 1, 1], M[0, 2, 2]],
                     [M[1, 0, 1], M[1, 0, 2], M[1, 1, 2], M[1, 0, 0], M[1, 1, 1], M[1, 2, 2]], 
                     [M[2, 0, 1], M[2, 0, 2], M[2, 1, 2], M[2, 0, 0], M[2, 1, 1], M[2, 2, 2]]])

def T4Reshape(M):
    return np.array([[M[0, 1, 0, 1], M[0, 1, 0, 2], M[0, 1, 1, 2], M[0, 1, 0, 0], M[0, 1, 1, 1], M[0, 1, 2, 2]], 
                     [M[0, 2, 0, 1], M[0, 2, 0, 2], M[0, 2, 1, 2], M[0, 2, 0, 0], M[0, 2, 1, 1], M[0, 2, 2, 2]], 
                     [M[1, 2, 0, 1], M[1, 2, 0, 2], M[1, 2, 1, 2], M[1, 2, 0, 0], M[1, 2, 1, 1], M[1, 2, 2, 2]], 
                     [M[0, 0, 0, 1], M[0, 0, 0, 2], M[0, 0, 1, 2], M[0, 0, 0, 0], M[0, 0, 1, 1], M[0, 0, 2, 2]], 
                     [M[1, 1, 0, 1], M[1, 1, 0, 2], M[1, 1, 1, 2], M[1, 1, 0, 0], M[1, 1, 1, 1], M[1, 1, 2, 2]], 
                     [M[2, 2, 0, 1], M[2, 2, 0, 2], M[2, 2, 1, 2], M[2, 2, 0, 0], M[2, 2, 1, 1], M[2, 2, 2, 2]]])

def AssignVal(T):

    global vm

    T = sp.Matrix(T)
    SymList = list(T.free_symbols)
    
    v = np.max(vm)
    
    c = 0
    for i in range(v+1,v+1+len(SymList)):
        T = T.subs(SymList[c],i)
        c = c + 1
  
    T = np.array(T,dtype=int)
    vm.append(np.max(T))
    
    return T

def PlotM(M):

    global SymString

    Mv = 1-np.log(np.abs(M[M !=0 ])+1)
    Mmin = np.min(Mv) - 0.5

    fig, ax = pl.subplots(figsize=(10, 8))
    
    ax.imshow(1-np.log(np.abs(M)+1), cmap='hot', vmin=Mmin, vmax=0.5)
    
    # create grids
    ax.set_xticks(np.arange(-.5, 17, 1), minor=True)
    ax.set_yticks(np.arange(-.5, 17, 1), minor=True)   
    ax.grid(which='minor', color='k', linestyle='-', linewidth=1)
    
    pl.hlines(y=np.array([2.5, 5.5, 11.5]), xmin=np.full(3, 0)-0.5, xmax=np.full(3, 17.5), color="k", linewidth=4)
    pl.vlines(x=np.array([2.5, 5.5, 11.5]), ymin=np.full(3, 0)-0.5, ymax=np.full(3, 17.5), color="k", linewidth=4)
    
    # create labels
    ax.set_xticks(np.arange(18))
    ax.set_yticks(np.arange(18))
    
    ylabel = ['$p_x$','$p_y$','$p_z$','$m_x$','$m_y$','$m_z$','$Q_{xy}$','$Q_{xz}$','$Q_{yz}$','$Q_{xx}$','$Q_{yy}$','$Q_{zz}$','$S_{xy}$','$S_{xz}$','$S_{yz}$','$S_{xx}$','$S_{yy}$','$S_{zz}$']
    ax.set_yticklabels(ylabel, fontsize=16)
    
    xlabel = ['$E_x$','$E_y$','$E_z$','$H_x$','$H_y$','$H_z$','$\diamond_{xy}^e$','$\diamond_{xz}^e$','$\diamond_{yz}^e$','$\diamond_{xx}^e$','$\diamond_{yy}^e$','$\diamond_{zz}^e$','$\diamond_{xy}^h$','$\diamond_{xz}^h$','$\diamond_{yz}^h$','$\diamond_{xx}^h$','$\diamond_{yy}^h$','$\diamond_{zz}^h$']
    ax.set_xticklabels(xlabel, fontsize=16)
    
    # create title
    pl.title("Material tensors for: {}".format(SymString), fontsize=16, fontweight="bold")
    
    # add annotations
    for i in range(len(M)):
        for j in range(len(M)):
            if M[i,j] !=0:
                text = ax.text(j, i, M[i, j], ha="center", va="center", color="k")
    
    pl.show(block = False)

def PlotT(T,T_val,n):

    global SymString
    global phi

    Mv = 1-np.log(np.abs(T_val[T_val !=0 ])+1)
    Mmin = np.min(Mv) - 0.5

    fig, ax = pl.subplots(figsize=(6, 5))
    
    ax.imshow(1-np.log(np.abs(T_val)+1), cmap='hot', vmin=Mmin, vmax=0.5)
    
    # create grids
    ax.set_xticks(np.arange(-.5, 4, 1), minor=True)
    ax.set_yticks(np.arange(-.5, 4, 1), minor=True)   
    ax.grid(which='minor', color='k', linestyle='-', linewidth=3)
    
    # remove ticks
    ax.set_xticks([])
    ax.set_yticks([]) 
    
    # create title
    if n:
        pl.title("Scattering matrix for: {} and $\phi$={}°\n at normal incidence".format(SymString,np.round(phi*180/np.pi,1)), fontsize=16, fontweight="bold")
    else:
        pl.title("Scattering matrix for: {} and $\phi$={}°\n at oblique incidence".format(SymString,np.round(phi*180/np.pi,1)), fontsize=16, fontweight="bold")

    # add annotations
    for i in range(len(T)):
        for j in range(len(T)):
            if T_val[i,j] !=0:
                text = ax.text(j, i, T[i, j], ha="center", va="center", color="k", fontweight="bold", fontsize=14)
    
    pl.show(block = False)


def FindSym(L):

    global vm

    # 2nd rank tensors
    ep = np.array(sp.symbols('a1:10')).reshape(3,3)
    mu = np.array(sp.symbols('b1:10')).reshape(3,3)
    xi = np.array(sp.symbols('c1:10')).reshape(3,3)
    ep = TSym2(ep)
    mu = TSym2(mu)

    # 3rd rank tensors
    Xpee = np.array(sp.symbols('d1:28')).reshape(3,3,3)
    Xpem = np.array(sp.symbols('e1:28')).reshape(3,3,3)
    Xpme = np.array(sp.symbols('f1:28')).reshape(3,3,3)
    Xpmm = np.array(sp.symbols('g1:28')).reshape(3,3,3)
    Xpee = TSym3(Xpee)
    Xpem = TSym3(Xpem)
    Xpme = TSym3(Xpme)
    Xpmm = TSym3(Xpmm)

    # 4th rank tensors
    Qpee = np.array(sp.symbols('h1:82')).reshape(3,3,3,3)
    Qpem = np.array(sp.symbols('i1:82')).reshape(3,3,3,3)
    Spmm = np.array(sp.symbols('j1:82')).reshape(3,3,3,3)
    Qpee = TSym4(Qpee)
    Qpee = TRep4(Qpee)
    Qpem = TSym4(Qpem)
    Spmm = TSym4(Spmm)
    Spmm = TRep4(Spmm)

    # Find invariant tensors that correspond to symmetries
    for q in L:
        
        ep = T2solve(q, ep, 1)
        mu = T2solve(q, mu, 1)
        xi = T2solve(q, xi, 0)
        
        Xpee = T3solve(q, Xpee, 1)
        Xpem = T3solve(q, Xpem, 0)
        Xpme = T3solve(q, Xpme, 0)
        Xpmm = T3solve(q, Xpmm, 1)
        
        Qpee = T4solve(q, Qpee, 1)
        Qpem = T4solve(q, Qpem, 0)
        Spmm = T4solve(q, Spmm, 1)
        
    # Reshape arrays
    Xpee = T3Reshape(Xpee)
    Xpem = T3Reshape(Xpem)
    Xpme = T3Reshape(Xpme)
    Xpmm = T3Reshape(Xpmm)
    Qpee = T4Reshape(Qpee)
    Qpem = T4Reshape(Qpem)
    Spmm = T4Reshape(Spmm)

    # Convert symbols into values
    vm = [0]
    ep   = AssignVal(ep)
    xi   = AssignVal(xi)
    mu   = AssignVal(mu)
    Xpee = AssignVal(Xpee)
    Xpme = AssignVal(Xpme)
    Qpee = AssignVal(Qpee)
    Xpem = AssignVal(Xpem)    
    Xpmm = AssignVal(Xpmm) 
    Qpem = AssignVal(Qpem)
    Spmm = AssignVal(Spmm)

    # Combine data into final matrix
    m1 = np.concatenate((ep,xi,Xpee,Xpem),axis=1)
    m2 = np.concatenate((-xi.transpose(),mu,Xpme,Xpmm),axis=1)
    m3 = np.concatenate((Xpee.transpose(),-Xpme.transpose(),Qpee,Qpem),axis=1)
    m4 = np.concatenate((-Xpem.transpose(),Xpmm.transpose(),-Qpem.transpose(),Spmm),axis=1)

    M = np.concatenate((m1,m2,m3,m4))

    # Plot
    PlotM(M)
    
def FindT(L):

    global vm
    global phi
    
    # initialize T matrix
    R11, R12, T13, T14 = sp.symbols('R11, R12, T13, T14') 
    R21, R22, T23, T24 = sp.symbols('R21, R22, T23, T24')
    T31, T32, R33, R34 = sp.symbols('T31, T32, R33, R34')
    T41, T42, R43, R44 = sp.symbols('T41, T42, R43, R44')
    
    T = np.array([[R11, R12, T13, T14], [R21, R22, T23, T24],
                  [T31, T32, R33, R34], [T41, T42, R43, R44]])
                  
    It = np.array([[1, 0],[0, 1]])
    
    POI = np.array([0, 1, 0]) # original orientation of the plane of incidence (xz-plane)
  
    # Find invariant scattering matrix that corresponds to symmetries
    for q in L:
    
        if q == 'i':
            q = 'iv'
        
        if q == 'c4z':
            q = 'c2z'
            if 'px' in L and 'py' in L:
                L.append("pxy")
                L.append("pxy2")
        
        Q = eval(q)
    
        Q = np.array(rz(-phi).dot(Q.dot(rz(phi))),dtype=float)
 
        m  = np.array([[Q[1,1] , 0],[0, Q[0,0]]])

        c1 = 0
        if Q[0,2] == 0 and Q[1,2] == 0:
            c1 = 1
        
        c2 = 0
        
        BPOI = rz(-phi).dot(rz(phi).dot(POI))
        BPOI[np.abs(BPOI) < 1e-5] = 0

        APOI = Q.dot(POI)
        APOI[np.abs(APOI) < 1e-5] = 0
        
        if (APOI == BPOI).all() or (APOI == -BPOI).all():
            c2 = 1
        
        c = c1*c2
        
        ap = c/2*(1 + Q[2,2])*m + (1 - c)*It
        am = c/2*(1 - Q[2,2])*m
        b  = c/2*(1 + Q[0,0])   +  1 - c
        
        m1 = np.concatenate((ap,am),axis=1)
        m2 = np.concatenate((am,ap),axis=1)
        
        M = np.concatenate((m1,m2))

        M[np.abs(M) < 1e-5] = 0
        
        eq = M.dot(b*T.dot(inv(M)) + (1-b)*T.transpose().dot(inv(M))) - T
        S = sp.solve(eq.flatten())
        
        if S:
            # reverse solutions
            Si = {}       
            for k, v in S.items():
                if v:         
                    if v_neg(v):
                        Si[-v] = -k
                    else:
                        Si[+v] = +k                   
                else:
                    Si[k] = v
            
            T = np.array(sp.Matrix(T).subs(Si)) 
            
        
    vm = [0]
    T_val = AssignVal(T)  

    PlotT(T,T_val,0)   


    
def FindTNI(L):

    global vm
    global phi
    
    # initialize T sp.Matrix
    R11, R12, T13, T14 = sp.symbols('R11, R12, T13, T14') 
    R21, R22, T23, T24 = sp.symbols('R21, R22, T23, T24')
    T31, T32, R33, R34 = sp.symbols('T31, T32, R33, R34')
    T41, T42, R43, R44 = sp.symbols('T41, T42, R43, R44')
    
    T = np.array([[R11, R12, T13, T14], [R21, R22, T23, T24],
                  [T31, T32, R33, R34], [T41, T42, R43, R44]])
                  
    It = np.array([[1, 0],[0, 1]])
    
    POI = np.array([0, 1, 0]) # original orientation of the plane of incidence (xz-plane)

    # at normal incidence, phi = 0
    phi = 0
  
    # Find invariant scattering sp.Matrix that corresponds to symmetries
    for q in L:
    
        if q == 'i':
            q = 'iv'
        
        if q == 'c4z':
            if 'px' in L and 'py' in L:
                L.append("pxy")
                L.append("pxy2")
                
            M = np.array([[0, -1, 0, 0],[1, 0, 0, 0],[0, 0, 0, -1],[0, 0, 1, 0]])
            
        else:
        
            Q = eval(q)
        
            Q = np.array(rz(-phi).dot(Q.dot(rz(phi))),dtype=float)
     
            m  = np.array([[Q[1,1] , 0],[0, Q[0,0]]])

            c1 = 0
            if Q[0,2] == 0 and Q[1,2] == 0:
                c1 = 1
            
            c2 = 0
            
            BPOI = rz(-phi).dot(rz(phi).dot(POI))
            BPOI[np.abs(BPOI) < 1e-5] = 0

            APOI = Q.dot(POI)
            APOI[np.abs(APOI) < 1e-5] = 0
            
            if (APOI == BPOI).all() or (APOI == -BPOI).all():
                c2 = 1
            
            c = c1*c2
            
            ap = c/2*(1 + Q[2,2])*m + (1 - c)*It
            am = c/2*(1 - Q[2,2])*m
            
            m1 = np.concatenate((ap,am),axis=1)
            m2 = np.concatenate((am,ap),axis=1)
            
            M = np.concatenate((m1,m2))

            M[np.abs(M) < 1e-5] = 0
        
        eq = M.dot(T.dot(inv(M))) - T
        S = sp.solve(eq.flatten())
        
        if S:
            # reverse solutions
            Si = {}       
            for k, v in S.items():
                if v:         
                    if v_neg(v):
                        Si[-v] = -k
                    else:
                        Si[+v] = +k                   
                else:
                    Si[k] = v
            
            T = np.array(sp.Matrix(T).subs(Si)) 
            
        
    vm = [0]
    T_val = AssignVal(T)  

    PlotT(T,T_val,1)   

        
def SymsT(SymS):

    global SymString

    # get symmetries from entry box
    SymString = ST.get()
    SymTerms = SymString.split()
    
    FindSym(SymTerms)


def SymsS(SymS):

    global SymString
    global phi

    # get symmetries from entry box
    SymString = ST.get()
    SymTerms = SymString.split()
    
    #get phi angle
    phi = eval(ST2.get())*np.pi/180
 
    FindT(SymTerms)
    FindTNI(SymTerms) # at normal incidence

# ===============================
# Create main window
# ===============================
ws = tk.Tk()
ws.title('Spatial symmetries, multipolar tensors and scattering matrix')
ws.geometry('600x100')

# close window by pressing on Esc
ws.bind("<Escape>", lambda e: ws.destroy())

# ===============================
# Create input frames
# ===============================
Input_frame = tk.Frame(ws)
Input_frame.pack(padx=10,pady=10, fill=tk.X)

lbl = tk.Label(Input_frame, text="List of symmetries: ", font=("Courier", 14))
lbl.grid(column=0,row=0,ipadx=5)

Input_frame2 = tk.Frame(ws)
Input_frame2.pack(padx=10,pady=10, fill=tk.X)

lbl2 = tk.Label(Input_frame2, text="Incidence plane orientation (°): ", font=("Courier", 14))
lbl2.grid(column=0,row=0,ipadx=5)

# actions when pressing enter
ST = tk.Entry(Input_frame, font=("Courier", 14))
ST.bind('<Return>',SymsT)
ST.bind('<KP_Enter>',SymsT)
ST.grid(column=1,row=0,sticky=tk.W+tk.E)

ST2 = tk.Entry(Input_frame2, font=("Courier", 14))
ST2.bind('<Return>',SymsS)
ST2.bind('<KP_Enter>',SymsS)
ST2.grid(column=1,row=0,sticky=tk.W+tk.E)
ST2.insert(0,0)

Input_frame.grid_columnconfigure(1, weight=1)
Input_frame2.grid_columnconfigure(1, weight=1)

ws.mainloop()
