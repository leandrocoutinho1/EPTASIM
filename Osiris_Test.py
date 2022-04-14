import os
import eptasim
from eptasim.rocketAtributes.rocketclass import rocket
from eptasim.plots.rocketplotfunc import rocketplot
from eptasim.ga.gaAtributes import ga_config, ga_run
from eptasim.flightsim.simAtributes import sim_config, Rocket_sim
import numpy
import matplotlib.pyplot as plt

###########################################  Rocket Definition ##################################
Osiris = rocket(dry_mass=4.068,
                initial_cg=702,
                diammeter=80,
                body_length=925,
                boattail_present=True)

Osiris.engine(thrust_curve='.\Teste 1km - Classe J .eng',
              prop_mass=0.76)

Osiris.nosecone(nosecone_type='elliptical',
                LD_ratio_nose=3)

Osiris.boattail(LD_ratio_boattail=0.5125,
                boattail_Aft_Diammeter=45)

####################################### Fins #######################################################

Osiris.fins(fins_number=4,
            fins_thickness=5,
            fins_dimensions=numpy.array([82, 53.46, 77.58, 44.96]))

rocketplot(Osiris, aoa=0)
