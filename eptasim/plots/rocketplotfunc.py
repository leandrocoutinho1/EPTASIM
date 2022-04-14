import matplotlib.pyplot as plt
import numpy as np
import math


def rocketplot(rocket, aoa=0, mach=0, t=0):
    try:

        ######################### Variables ##########################################
        Bodytube_Lenght = rocket.body_length
        Bodytube_Diammeter = rocket.diammeter
        Nosecone_type = rocket.nosecone_type
        LD_Ratio = rocket.LD_ratio_nose
        Fins_Number = rocket.fins_number
        Boattail = rocket.boattail_present
        Nosecone_Length = rocket.nosecone_length

        # Rocket Actual CG
        Actual_CG = rocket.cg(t)

        f1 = plt.figure(1, figsize=(20, 10))
        plt.axes().set_aspect('equal')

        ######################## Nosecone ############################################

        if Nosecone_type == 'cone':  # Cone
            X_NoseCone = np.array([0, Nosecone_Length, Nosecone_Length, 0])
            Y_NoseCone = np.array(
                [0, 0.5*Bodytube_Diammeter, -0.5*Bodytube_Diammeter, 0])

        elif Nosecone_type == 'ogive':  # Ogive
            pass

        elif Nosecone_type == 'parabolic':  # Parabolic
            pass

        elif Nosecone_type == 'elliptical':  # Ellipsoid
            T = np.linspace(90, 270)
            T = T * math.pi/180
            X_NoseCone = (Nosecone_Length * np.cos(T)) + Nosecone_Length
            Y_NoseCone = (0.5 * Bodytube_Diammeter) * np.sin(T)

        plt.plot(X_NoseCone, Y_NoseCone)

        ######################## Bodytube ############################################

        X_Bodytube = np.array([Nosecone_Length, Nosecone_Length
                               + Bodytube_Lenght, Nosecone_Length
                               + Bodytube_Lenght, Nosecone_Length, Nosecone_Length])                                       # X values Bodytube
        Y_Bodytube = np.array([0.5*Bodytube_Diammeter, 0.5*Bodytube_Diammeter, -0.5*Bodytube_Diammeter, -
                               0.5*Bodytube_Diammeter, 0.5*Bodytube_Diammeter])                                 # Y Values Bodytube

        # Plot Bodytube
        plt.plot(X_Bodytube, Y_Bodytube)

        ########################### Fins #############################################
        root_chord = rocket.fins_dimensions[0]
        tip_chord = rocket.fins_dimensions[1]
        semispan = rocket.fins_dimensions[2]
        sweep_angle = rocket.fins_dimensions[3] * math.pi/180  # to radians
        sweep_distance = semispan * math.tan(sweep_angle)

        ### Upper fins ####

        X_Fins = np.array([Nosecone_Length + Bodytube_Lenght-root_chord,
                           Nosecone_Length + Bodytube_Lenght-root_chord + sweep_distance,
                           Nosecone_Length + Bodytube_Lenght-root_chord + sweep_distance + tip_chord,
                           Nosecone_Length + Bodytube_Lenght])                          # X coordinates Fins
        Y_Upper_Fins = np.array([0.5*Bodytube_Diammeter,
                                0.5*Bodytube_Diammeter+semispan,
                                0.5*Bodytube_Diammeter+semispan,
                                0.5*Bodytube_Diammeter])                               # Y coordinates Upper Fins
        Y_Lower_Fins = np.array([-0.5*Bodytube_Diammeter,
                                -0.5*Bodytube_Diammeter-semispan,
                                -0.5*Bodytube_Diammeter-semispan,
                                -0.5*Bodytube_Diammeter])                             # Y coordinates Lower fins

        plt.plot(X_Fins, Y_Upper_Fins)
        plt.plot(X_Fins, Y_Lower_Fins)

        ####################### Boattail #############################################

        if Boattail:
            D1_Boattail = Bodytube_Diammeter
            # Boattail Length
            Boattail_Length = rocket.boattail_length
            # Boattail Lower CG
            D2_Boattail = rocket.boattail_Aft_Diammeter
            X_Boattail = np.array([Nosecone_Length + Bodytube_Lenght,
                                   Nosecone_Length + Bodytube_Lenght + Boattail_Length,
                                   Nosecone_Length + Bodytube_Lenght + Boattail_Length,
                                   Nosecone_Length + Bodytube_Lenght])                    # X coordinates Boattail

            # Y coordinates Boattail
            Y_Boattail = np.array(
                [0.5*D1_Boattail, 0.5*D2_Boattail, -0.5*D2_Boattail, -0.5*D1_Boattail])

            plt.plot(X_Boattail, Y_Boattail)

        else:
            Boattail_Length = 0

        ############################ Plot ############################################

        #CG#################
        plt.plot(Actual_CG, 0, 'r+', label="CG: %.2f" % Actual_CG)
        #CP################
        X_CP = (rocket.static_margin(mach, aoa)
                * Bodytube_Diammeter) + Actual_CG
        SM = rocket.static_margin(mach, aoa)
        plt.annotate('Static Margin: %.2f' % SM, xy=(0, 0), xytext=(0, -90))

        plt.plot(X_CP, 0, 'bo', label="CP: %.2f" % X_CP)

        ##################

        plt.ylim(-Bodytube_Diammeter - semispan, Bodytube_Diammeter + semispan)
        plt.legend()
        plt.show()
    except:
        print('### Rocket not completly defined ###')
