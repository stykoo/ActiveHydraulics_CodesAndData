"""
Sub module for  finite element simulations of Toner-Tu equations
in oblong geometry.

Author:
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

import os
import sub_sub_simu_oblong_parallel_v2

def para(SIZEX, SIZEY, N, DT, LAMBDA, LAMBDA_2, LAMBDA_3, FAC_D_RHO, FAC_D_V,
         SIGMA, SIGMA_2, ALPHA, BETA, TERM_1, TERM_2, TERM_3, TERM_SURF,
         NUM_STEP, NB_SOL, i):
    # crée le dossier où seront stockés les données
    directory = \
        "solution" + "_D="+ str(FAC_D_V*1e-2) + "_sigma=" + str(SIGMA) \
        + "_lambda=" + str(LAMBDA) + "_L=" + str(SIZEX) + "_W=" + str(SIZEY) \
        + "_DT=" + str(DT) + "_lambda2=" + str(LAMBDA_2) +"_lambda3=" \
        + str(LAMBDA_3)+ "_sigma2=" + str(SIGMA_2) + "_it=" + str(i+1)
    
    
    if os.path.exists(directory)==False:
        os.mkdir(directory)
        #lance la simu
        sub_sub_simu_oblong_parallel_v2.run_simu(
                N, SIZEX, SIZEY, DT, LAMBDA, LAMBDA_2, LAMBDA_3, FAC_D_RHO,
                FAC_D_V, SIGMA, SIGMA_2, ALPHA, BETA, TERM_1, TERM_2, TERM_3,
                TERM_SURF, NUM_STEP, NB_SOL, directory, i)
                                                
