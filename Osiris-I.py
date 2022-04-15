import eptasim
from eptasim.rocketAtributes.rocketclass import rocket
from eptasim.plots.rocketplotfunc import rocketplot
from eptasim.ga.gaAtributes import ga_config, ga_run
from eptasim.flightsim.simAtributes import sim_config, Rocket_sim
import numpy
import matplotlib.pyplot as plt


########################################### Rocket Definition ##################################
## Input de dados para construção do Foguete ##

Osiris = rocket(dry_mass=5 ,               # Massa vazia do Foguete
                initial_cg= 716,           # Cg inicial
                diammeter= 80,             # Diâmetro do Foguete
                body_length=919,           # Comprimento do Foguete
                boattail_present= True)    # Presença ou não de boattail

## Curva de Empuxo ##
Osiris.engine(thrust_curve = "Thrust curve\Teste 1km - Classe J .eng", prop_mass=0.76)

Osiris.nosecone(nosecone_type = 'elliptical',  ## Formato do nosecone ##
                LD_ratio_nose=3)               ## Razão LD do nosecone (comprimento / diâmetro) ##

Osiris.boattail(LD_ratio_boattail=0.4375,      ## Razão LD do boattail (comprimento / diâmetro) ##
                boattail_Aft_Diammeter= 58.3)  ## Comprimento do boattail ##
        


######################################## Optimization #############################################
## Definição dos parâmetros que serão utilizados para otimização das aletas ##

Osiris_config = ga_config(Osiris)
Osiris_config.var_boundaries(root_chord=numpy.array([80,150]),     # Corda da raiz
                            tip_chord= numpy.array([40,90]),       # Corda da ponta
                            span= numpy.array([75,90]),            # Envergadura
                            sweep_angle=numpy.array([30,40]),      # Ângulo de Enflechamento
                            static_margin_min=1.7,                 # Margem Estática Mínima
                            static_margin_max= 1.9,                # Margem Estática Máxima
                            rail_exit_speed= 36,                   # Velocidade de Saída do Trilho
                            crosswind= 0)                          # Vento

Osiris_config.algorithm_param(max_num_iteration=50,                # Número de Iterações
                           population_size= 500)                   # Número de Indivíduos por Iteração

fins = ga_run(Osiris_config)




####################################### Fins #######################################################

Osiris.fins(fins_number= 4,
            fins_thickness= 5,
            fins_dimensions= fins)
# First   [80.21577366 59.90776674 77.14089327 39.9997227 ]

# Objective function:
# 0.5708552

rocketplot(Osiris, aoa= 0)

print(Osiris.static_margin(0, 10))






