import numpy
from geneticalgorithm import geneticalgorithm as ga
from eptasim.flightsim.simAtributes import sim_config, Rocket_sim
import math


class ga_config:
    def __init__(self, rocket):
        self.rocket = rocket

    def var_boundaries(self, root_chord, tip_chord, span, sweep_angle, static_margin_min=1.5, static_margin_max=2, rail_exit_speed=30, crosswind=4) -> numpy.array:
        self.varbound = numpy.array([root_chord, tip_chord, span, sweep_angle])
        self.static_margin_min = static_margin_min
        self.static_margin_max = static_margin_max
        self.rail_exit_speed = rail_exit_speed
        self.crosswind = crosswind

    def algorithm_param(self, max_num_iteration=1000, population_size=500, mutation_probability=0.08,
                        elit_ratio=0.01, crossover_probability=0.5, parents_portion=0.3,
                        crossover_type='uniform', max_iteration_without_improv=None, function_timeout=100):

        self.algorithm_param = {'max_num_iteration': max_num_iteration,
                                'population_size': population_size,
                                'mutation_probability': mutation_probability,
                                'elit_ratio': elit_ratio,
                                'crossover_probability': crossover_probability,
                                'parents_portion': parents_portion,
                                'crossover_type': crossover_type,
                                'max_iteration_without_improv': max_iteration_without_improv,
                                'function_timeout': function_timeout}                                               # This array will be given to the GA function


def ga_run(config_object):
    rocket = config_object.rocket

    def cost(X):
        rocket.fins(fins_number=4, fins_thickness=3, fins_dimensions=X)
        sm = rocket.static_margin(mach=0, aoa=0)

        if sm < config_object.static_margin_min or sm > 1.2*config_object.static_margin_max or X[0] <= 1.2*X[1]:
            return math.inf
        elif sm >= config_object.static_margin_min:
            top_mach = 173

            if config_object.crosswind >0:
                rail_aoa = math.degrees(math.atan(
                    config_object.rail_exit_speed/config_object.crosswind))
            else:
                rail_aoa = 0

            SM_rail = rocket.static_margin(0, aoa=rail_aoa)

            if SM_rail >= config_object.static_margin_min:
                    

                drag_1 = rocket.drag(0.2*top_mach)
                drag_2 = rocket.drag(0.4*top_mach)
                drag_3 = rocket.drag(0.6*top_mach)
                drag_4 = rocket.drag(0.8*top_mach)
                drag_5 = rocket.drag(top_mach)
                return (drag_1[0]*0.08) + (drag_2[0]*0.16) + (0.33*drag_3[0]) + (0.27*drag_4[0]) + (0.16*drag_5[0])

            else:
                return math.inf

    model = ga(function=cost,
               dimension=4,
               variable_type='real',
               variable_boundaries=config_object.varbound,
               algorithm_parameters=config_object.algorithm_param)                                              # Model config

    # Run
    model.run()

    # Best results
    return model.output_dict['variable']
