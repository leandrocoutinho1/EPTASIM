import numpy
import numpy.linalg as la
import math


class sim_config:                                                       # Config object
    def __init__(self, rocket, t0=0, tf=20, dt=0.01):
        self.rocket = rocket
        self.t0 = t0
        self.tf = tf
        self.dt = dt

    def launch_options(self, Launch_angle=0, Launch_Heading=0):
        self.Launch_angle = Launch_angle
        self.Launch_Heading = Launch_Heading

    def wind(self, Wind_speed=0, Wind_direction=0):
        self.Wind_speed = Wind_speed
        self.Wind_direction = Wind_direction


def Rocket_sim(sim_config):

    # Rocket attitude vector paralel to the wind
    def Rocket_attitude(rocket, v, t):

        # Speed
        v_norm = la.norm(v)
        # Attitude vector (Unit Vector)
        Attitude_vector = (v/v_norm)

        return Attitude_vector

    # Rocket Drag Force
    def Drag_force(rocket, v, t, alt):
        rocket_area = math.pi * ((rocket.diammeter/2000)**2)  # Area in m^2
        attitude = Rocket_attitude(rocket, v, t)
        Speed = la.norm(v)
        cd_rho = rocket.drag(alt, Speed)[0]
        Drag_vector = -0.5 * cd_rho * math.pow(Speed,
                                               2) * rocket_area * attitude

        return Drag_vector

    def Thrust(rocket, v, t):
        Thrust_vector = Rocket_attitude(rocket, v, t) * rocket.Thrust(t)
        return Thrust_vector

    def acceleration(rocket, t, v, alt):
        if alt > 0:                                                                 # Rocket on the ground ?
            g = numpy.array([0, 0, -9.81])
        # Gravity turns on only after rocket leaves the ground
        else:
            g = numpy.array([0, 0, 0])

        mass = rocket.mass(t)

        Acceleration_vector = g + ((Thrust(rocket, v, t)/mass) +
                                   (Drag_force(rocket, v, t, alt)/mass))

        return Acceleration_vector

    # RungeKutta
    def RK4(rocket, f, v0, t0, tf, dt):
        # Time Array (From t0 to tf)
        t = numpy.arange(t0, tf, dt)
        # Time intervals number
        nt = t.size
        # number of variables in V0 vector (Usually 3)
        n_row, n_col = v0.shape
        v = numpy.zeros((n_row, nt))
        # Space Vector
        x = numpy.zeros((n_row, nt))

        # Initial Velocity Vector
        v[:, 0] = v0[:, 0]

        # Calculations
        for k in range(0, nt-1):
            alt = x[2, k]
            k1 = dt * f(rocket, t[k], v[:, k], alt)
            k2 = dt * f(rocket, t[k] + dt/2, v[:, k] + k1/2, alt)
            k3 = dt * f(rocket, t[k] + dt/2, v[:, k] + k2/2, alt)
            k4 = dt * f(rocket, t[k] + dt, v[:, k] + k3, alt)

            dv = (k1 + (2*k2) + (2*k3) + (k4)) / \
                6                                    # diferential

            # State vector in a given t
            v[:, k+1] = v[:, k] + dv

            dx = v[:, k]*dt

            # Space Vector
            x[:, k+1] = x[:, k] + dx
            if x[2, k+1] <= 0:  # Hits gorund
                break
        return x, v, t

    #################### SIM CONFIG #########################

    try:
        # Launch angle (0° = vertical)
        Launch_angle = sim_config.Launch_angle
        # Launch Heading (0° = North)
        Launch_Heading = sim_config.Launch_Heading
        # Initial Launch Speed > 0 otherwise = Bug
        Launch_speed = 0.0001
        t0 = sim_config.t0
        tf = sim_config.tf
        dt = sim_config.dt
        rocket = sim_config.rocket
    except:
        print("### No launch angle or heading data ####")
        return

    ###################### DATA ARAYS #########################

    vx = (Launch_speed * math.sin(math.radians(Launch_angle))
          * math.cos(math.radians(Launch_Heading)))                                 # To polar coordinates
    vy = (Launch_speed * math.sin(math.radians(Launch_angle))
          * -math.sin(math.radians(Launch_Heading)))                                # To Polar Coordinates
    vz = (Launch_speed * math.cos(math.radians(Launch_angle))
          )                      # To polar coordinates

    v0 = numpy.array([[vx],
                      [vy],
                      [vz]])                                                            # Initial Velocity

    x, v, t = RK4(rocket, acceleration, v0, t0, tf, dt)

    return x, v, t
