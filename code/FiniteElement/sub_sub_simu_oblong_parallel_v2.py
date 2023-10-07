#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Use the python interface to FEniCS for finite element simulations of
the Toner-Tu equations in oblong geometry.

Authors:
    Camille Jorge <camille.jorge@ens-lyon.fr>
    Yoann Poupart <yoann.poupart@ens-lyon.fr>

Licence:
    Copyright (C) 2023 ENS de Lyon
    Contributors: 
        Camille Jorge <camille.jorge@ens-lyon.fr>
        Yoann Poupart <yoann.poupart@ens-lyon.fr>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
     
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
     
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from fenics import *
from dolfin import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt
import random


def fig_defaut(u, p, ifig=-1, plotfigure=False, nom='no_name'):
    if plotfigure:
        plt.figure(figsize=(13, 13))

        pl = plot(project(u / dot(u, u) ** 0.5))
        pl.set_clim(0.0, 10.0)
        plot(p, alpha=0.3)

        axes = plt.gca()
        axes.set_xlabel("x (mm)", fontsize=14)
        axes.set_ylabel("y (mm)", fontsize=14)
        plt.savefig(f"./Image/{nom}_v_{ifig:05}.png", bbox_inches='tight', format='png')
        plt.close()


def get_value(x, y, XX, YY, L):
    for i in range(len(L)):
        if abs(XX[i] - x) < 1e-3 and abs(YY[i] - y) < 1e-3:
            return L[i]
    print(x, y)
    raise ValueError


def get_value_int(x, y, XX, YY, L):
    for i in range(len(L)):
        if abs(XX[i] - 1 - x) < 1e-3 and abs(YY[i] - 1 - y) < 1e-3:
            return L[i]
    print(x, y)
    raise ValueError


class RandomVector(UserExpression):
    """ 
    Class based on UserExpression to generate random numbers for random vectors
    using Python. 
    """
    # def __init__(self,r=2, degree=2, **kwargs):
    def __init__(self,r=1.5, degree=2,n=1,**kwargs):
        # random.seed(2 + MPI.rank(MPI.comm_world) )
        random.seed(2 + n )
        super().__init__(**kwargs)
        self.r = r
    def eval(self, values, x):
        # values[0] = ((random.rand() - 0.5) * self.r, (random.rand() - 0.5) * self.r)
        values[0] = (random.random() - 0.5) * self.r
        values[1] = (random.random() - 0.5) * self.r
    def value_shape(self):
        return (2,)
        

# Geometry of the mesh 
def make_outer(N,WIDTH,L):
    
    import numpy as np
    
    x = [WIDTH/2*np.cos(2*np.pi/N*n-np.pi/2) + (2*np.pi/N*n-np.pi/2 < np.pi/2)*L for n in range(N+1)]
    y = [WIDTH/2*np.sin(2*np.pi/N*n-np.pi/2)  for n in range(N+1)]
    
    vertices = [Point(x[i],y[i]) for i in range(N+1)]
    
    return vertices


def run_simu(N,SIZEX,SIZEY,DT,LAMBDA,LAMBDA_2,LAMBDA_3,FAC_D_RHO,FAC_D_V,SIGMA,SIGMA_2,ALPHA,BETA,TERM_1,TERM_2,TERM_3,TERM_SURF,NUM_STEP,NB_SOL,directory,i):
    # définition de la simulation
    # NX = int(SIZEX / DX)# Nombre d'élément par ligne
    # NY = int(SIZEY / DY)# Nombre d'élément par colonne
    DEG_V = 2  # degré des éléments V
    DEG_P = 2  # degré des éléments rho
    

    # définition de la sauvegarde
    NB_IM = 5  # 100
    NOM = 'TT_'
    
    # définition de la physique
    RHO_0 = Constant(0.3)

    # Define the mesh
    domain = Polygon(make_outer(N,SIZEY,SIZEX))
    mesh = generate_mesh(domain,20)


    V1 = VectorElement("Lagrange", mesh.ufl_cell(), DEG_V)
    V2 = FiniteElement("Lagrange", mesh.ufl_cell(), DEG_P)
    V = VectorFunctionSpace(mesh, 'P', DEG_V)
    Vs = FunctionSpace(mesh, 'P', DEG_P)

    # Make a mixed space
    TH = V1 * V2
    W = FunctionSpace(mesh, TH)
    
    X = Expression(('x[0]', 'x[1]'), element=V.ufl_element())
    xv = interpolate(X, V)
    
    xp = interpolate(X, VectorFunctionSpace(mesh, 'P', DEG_P))
    np.savetxt(directory+"/coordo_v.txt", np.reshape(np.array(xv.vector()), (-1, 2)))
    np.savetxt(directory+"/coordo_p.txt", np.reshape(np.array(xp.vector()), (-1, 2)))

    # Initial conditions
    print("Defining initial conditions...")
    
    
    r = 2
    
    # u_in = Expression(('r*(rand()-0.5)', 'r*(rand()-0.5)'), degree=DEG_V, r=r)
    u_in = RandomVector(r=r, degree=DEG_V,n=i)
    print(f' Value shape {u_in.value_shape}')
    u_0 = interpolate(u_in, W.sub(0).collapse())

    p_in = Expression('RHO_0', degree=0, RHO_0=RHO_0)
    p_0 = interpolate(p_in, W.sub(1).collapse())
    
    # Declare functions
    w = Function(W)
    w_1 = Function(W)

    assign(w, [u_0, p_0])
    assign(w_1, [u_0, p_0])

    w_t = TestFunction(W)

    (u, p) = split(w)
    (u_1, p_1) = split(w_1)
    (v, q) = split(w_t)

    # Declare constants
    n = FacetNormal(mesh) #vecteur normal

    lamb = Constant(LAMBDA)
    lamb_2 = Constant(LAMBDA_2)
    lamb_3 = Constant(LAMBDA_3)  # SD
    alph = Constant(ALPHA)  # s-1
    rho_c = Constant(3e-3)  # SD
    a2 = Constant(alph * (RHO_0 - rho_c))  # s-1
    a4 = Constant(BETA)  # mm-2 s
    V_0 = (a2 / a4) ** 0.5
    D = Constant(1e-2 * FAC_D_V)  # mm2 s-1
    D2 = Constant(4e-6 * FAC_D_RHO)  # mm2 s-1
    sigma = Constant(SIGMA)  # mm2 s-2
    sigma_2 = Constant(SIGMA_2)
    Dt = Constant(DT)
    sigma_penalty = Constant(2)  # penalty parameter
    h = MaxCellEdgeLength(mesh)

    # Declare the form
    print("Defining the variational form...")
    
            
    F = dot((u - u_1) / Dt, v) * dx \
        + lamb * dot(dot(u, nabla_grad(u)), v) * dx \
        + D * (dot(grad(u[0]), grad(v[0])) + dot(grad(u[1]), grad(v[1]))) * dx \
        - (alph * (p - rho_c) - a4 * dot(u, u)) * dot(u, v) * dx \
        + sigma * dot(grad(p), v) * dx \
        + (((p - p_1) / Dt + div(p * u)) * q)*dx + 1./h**sigma_penalty*inner(dot(u, n), dot(v, n))*ds 
        #+ D2 * dot(grad(p), grad(q))) * dx \
        #- D2 * dot(grad(p), n) * q * ds


    if TERM_1:
        F += lamb_2 * div(u) * dot(u, v) * dx
    if TERM_2:
        F += lamb_3 * dot(grad(dot(u, u)), v) * dx
    if TERM_3:
        F += sigma_2 * dot(u, v) * dot(u, grad(p)) * dx
    if TERM_SURF:
    	F += -D*(dot(grad(u[0]),n)*v[0]+dot(grad(u[1]),n)*v[1])*ds \
    	     -D2*dot(grad(p),n)*q*ds

    print("Begining simulations...")
    for i in range(0, NUM_STEP):
        plotfigure = False
        
        solve(F == 0, w)
        assign(w_1, w)

        p_ana = project(w[2], Vs)
        u_ana = project(as_vector([w[0], w[1]]), V)

        if i % NB_SOL == 0:
            np.savetxt(directory+"/soluce_v_"+ str(i) + ".txt", np.reshape(np.array(u_ana.vector()), (-1, 2)))
            np.savetxt(directory+"/soluce_p_"+ str(i) + ".txt", p_ana.vector()[:])
