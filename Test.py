import eptasim
from eptasim.rocketAtributes.rocketclass import rocket
from eptasim.plots.rocketplotfunc import rocketplot
from eptasim.ga.gaAtributes import ga_config, ga_run
from eptasim.flightsim.simAtributes import sim_config, Rocket_sim
import numpy
import matplotlib.pyplot as plt


Vesper = rocket(dry_mass=2, initial_cg=700,
                diammeter=80, body_length=1000)     #
Vesper.nosecone(nosecone_type='elliptical', LD_ratio_nose=4)


Vesper.engine('D:\Google Drive\Programação\Python\EPTASIM\I348 - Teste 10 - 1km -13-08-21.eng', prop_mass=0.25)

Vesper.boattail(LD_ratio_boattail=1.2, boattail_Aft_Diammeter=50)


# fins = ga_run(config_1)

Vesper.fins(fins_number=4, fins_thickness=3,
            fins_dimensions=numpy.array([100, 50, 100, 45]))  # Incicializa a aleta


config_1 = ga_config(Vesper)
config_1.var_boundaries(root_chord=numpy.array([90, 130]), tip_chord=numpy.array(
    [50, 100]), span=numpy.array([40, 100]), sweep_angle=numpy.array([40, 60]),
    static_margin_min=1.5, static_margin_max=2, rail_exit_speed=30, crosswind=3)
config_1.algorithm_param(max_num_iteration=5000, population_size=1000)

fins = ga_run(config_1)

Vesper.fins(fins_number=4, fins_thickness=3, fins_dimensions=fins)
rocketplot(Vesper)

Simulacao_1 = sim_config(Vesper, t0=0, tf=30)

Simulacao_1.launch_options()

x, v, t = Rocket_sim(Simulacao_1)

f1 = plt.plot(t, x[2, :])
plt.show()
