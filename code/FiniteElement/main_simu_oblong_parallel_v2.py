#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Launch finite element simulations of Toner-Tu equations in oblong geometry.

Authors:
    Camille Jorge <camille.jorge@ens-lyon.fr>
    Yoann Poupart <yoann.poupart@ens-lyon.fr>

Usage:
    python3 main_simu_oblong_parallel_v2.py

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

import sub_simu_oblong_parallel_v2
from joblib import Parallel, delayed
from multiprocessing import Pool
import numpy as np

if __name__ == '__main__':
    # définition de la physique
    DT = 1e-3  # pas de temps
    N=50#
    SIZEY = 0.2
    
    
    FAC_D_RHO = 1
    FAC_D_V = 0.2 #defaut 1, coefficient pour Quincke
    ALPHA = 100  # defaut 100
    BETA = 10  # defaut 10
    LAMBDA = 0.7
    SIGMA = 5  # defaut 5
    LAMBDA_2 = 0.17 # default 0.7, coefficient pour Quincke
    LAMBDA_3 = -0.17 # default 0.7, coefficient pour Quincke
    SIGMA_2 = 0.1# defaut 5
    
    TERM_SURF = False
    TERM_1 = True
    TERM_2 = True
    TERM_3 = True

    # définition de la sauvegarde
    NB_SOL = 500  # 10
    
    SIZEX_LISTE = [0.2, 0.18, 0.16, 0.14, 0.12, 0.1, 0.08, 0.06, 0.04, 0.02, 0] 
    NUM_STEP = 3001
    it = 100
    
    # sub_simu_rectangle_parallel_args = [(i,) for i in range(it)]
    for SIZEX in SIZEX_LISTE:
        Parallel(n_jobs=10)(
            delayed(sub_simu_oblong_parallel_v2.para)(
                SIZEX, SIZEY, N, DT, LAMBDA, LAMBDA_2, LAMBDA_3, FAC_D_RHO,
                FAC_D_V, SIGMA, SIGMA_2, ALPHA, BETA, TERM_1, TERM_2, TERM_3,
                TERM_SURF, NUM_STEP, NB_SOL, i
            )
            for i in range(it)
        )
            
