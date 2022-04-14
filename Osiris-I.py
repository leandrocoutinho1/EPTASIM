import eptasim
from eptasim.rocketAtributes.rocketclass import rocket
from eptasim.plots.rocketplotfunc import rocketplot
from eptasim.ga.gaAtributes import ga_config, ga_run
from eptasim.flightsim.simAtributes import sim_config, Rocket_sim
import numpy
import matplotlib.pyplot as plt


########################################### Rocket Definition ##################################
Osiris = rocket(dry_mass=5 , 
                initial_cg= 716, 
                diammeter= 80, 
                body_length=919,    
                boattail_present= True)

Osiris.engine(thrust_curve = "Thrust curve\Teste 1km - Classe J .eng", prop_mass=0.76)

Osiris.nosecone(nosecone_type = 'elliptical',
                LD_ratio_nose=3)

Osiris.boattail(LD_ratio_boattail=0.4375, 
                boattail_Aft_Diammeter= 58.3)
        


######################################## Optimization #############################################

Osiris_config = ga_config(Osiris)
Osiris_config.var_boundaries(root_chord=numpy.array([80,150]),
                            tip_chord= numpy.array([40,90]),
                            span= numpy.array([75,90]),
                            sweep_angle=numpy.array([30,40]),
                            static_margin_min=1.7,
                            static_margin_max= 1.9,
                            rail_exit_speed= 36,
                            crosswind= 0)

Osiris_config.algorithm_param(max_num_iteration=50,
                           population_size= 500)

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






