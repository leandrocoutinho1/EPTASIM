import eptasim
from eptasim.rocketAtributes.rocketclass import rocket
from eptasim.plots.rocketplotfunc import rocketplot
from eptasim.ga.gaAtributes import ga_config, ga_run
from eptasim.flightsim.simAtributes import sim_config, Rocket_sim
import numpy
import matplotlib.pyplot as plt

########################################### Rocket Definition ##################################
## Input de dados para construção do Foguete ##

Ankaa = rocket(dry_mass=4.961,  # Massa vazia do Foguete
                initial_cg=607,  # Cg inicial
                diammeter=89,  # Diâmetro do Foguete
                body_length=1050,  # Comprimento do Foguete
                boattail_present=True)  # Presença ou não de boattail

## Curva de Empuxo ##
Ankaa.engine(thrust_curve="Thrust curve\meteor-RASP_Projeto ANKAA 9 def.eng", prop_mass=0.897)

Ankaa.nosecone(nosecone_type='elliptical',            ## Formato do nosecone ##
                LD_ratio_nose=1.12359)                ## Razão LD do nosecone (comprimento / diâmetro) ##

Ankaa.boattail(LD_ratio_boattail=1.3727,              ## Razão LD do boattail (comprimento / diâmetro) ##
                boattail_Aft_Diammeter=53.6334)       ## Comprimento do boattail ##

## Definição dos parâmetros que serão utilizados para otimização das aletas ##

Ankaa_config = ga_config(Ankaa)
Ankaa_config.var_boundaries(root_chord=numpy.array([90, 170]),    # Corda da raiz
                             tip_chord=numpy.array([50, 110]),     # Corda da ponta
                             span=numpy.array([75, 130]),         # Envergadura
                             sweep_angle=numpy.array([20, 50]),    # Ângulo de Enflechamento
                             static_margin_min=1.7,                # Margem Estática Mínima
                             static_margin_max=1.9,                # Margem Estática Máxima
                             rail_exit_speed=36,                   # Velocidade de Saída do Trilho
                             crosswind=0)                          # Vento

Ankaa_config.algorithm_param(max_num_iteration=50,  # Número de Iterações
                              population_size=500)  # Número de Indivíduos por Iteração

fins = ga_run(Ankaa_config)

####################################### Fins #######################################################

Ankaa.fins(fins_number=4,
            fins_thickness=5,
            fins_dimensions=fins)




####################################### Results #######################################################
# First   [104.19182795  30.41668783  60.64777282  59.47422396]
# Objective function:
# 0.4343589

# Second   [65.53322904 33.16110239 60.27849458 39.57023285]
# Objective function:
# 0.46979540000000003

# Third    [106.30650161  25.07835107  60.83600627  59.54896439]
# Objective function:
# 0.4338089

# Fourth   [90.31422119 36.05386849 62.16880296 59.55674173]
# Objective function:
# 0.4349896000000001

# Fifth   [107.56591859  35.10567245  60.17071182  59.92493509]
# Objective function:
# 0.4342865

# Sixth   [100.2096026   35.1095293   61.00012768  59.93761814]
# Objective function:
# 0.43417780000000006

# Seventh   [93.15795645 41.75630846 61.18205772 59.50787941]
# Objective function:
# 0.43533479999999997

# Eighth  [79.67715271 57.56616961 61.85939138 59.74159999]
# Objective function:
# 0.4365631

# Ninth    [113.53935885  20.76352364  58.54357396  59.7019139 ]
# Objective function:
# 0.4309731

# Tenth    [106.12254492  40.56840565  58.09663524  59.62529959]
# Objective function:
# 0.432885

# Eleventh    [85.63999122 58.81795498 55.26232903 39.46894078]
# Objective function:
# 0.46768319999999997

# Twelfth    [75.00132758 40.0639049  60.04048661 39.99613397]
# Objective function:
# 0.47115550000000006

#Thirteenth   [92.47913329 35.45309509 59.83137392 59.23324427]
# Objective function:
# 0.4328614

# Fourteenth    [93.1000256  40.70395759 61.38452535 59.71848659]
# Objective function:
# 0.43499479999999996
#######################################################################################################

rocketplot(Ankaa, aoa=0)

#print(Ankaa.static_margin(0, 10))






