from scipy.interpolate import interp1d
import numpy
import math


########## Variaveis Troposfericas ##########
# gravidade local (constante mudar caso seja necessario no futuro) [m/s^2]
global atmo_g
global atmo_M
# taxa de gradiente adiabatico (negativado para contas)
global atmo_L_0
global atmo_R
global atmo_rho
global atmo_T
global atmo_P
global atmo_data

atmo_g = 9.80665
atmo_M = 0.0289644
atmo_L_0 = -0.0065
atmo_R = 8.31447  # constante do gás ideal
atmo_rho = 1.22497705587732  # densidade CNTP do ar seco
atmo_T = 288.15  # Temperatura CNTP  [K]
atmo_P = 101325  # Pressão CNTP [Pa]
atmo_data = numpy.array([atmo_g,
                         atmo_M,
                         atmo_L_0,
                         atmo_R,
                         atmo_rho,
                         atmo_T,
                         atmo_P])


# Object rocket as argument
def fins_volume(rocket):

    try:                                                                           # Checking if fins are defined
        fins_data = rocket.fins_dimensions
        fins_thickness = rocket.fins_thickness
        fins_number = rocket.fins_number
        root_chord = fins_data[0]
        tip_chord = fins_data[1]
        semispan = fins_data[2]
        # to radians
        sweep_angle = fins_data[3] * math.pi/180
        # Distace between root leading edge to tip LE
        sweep_distance = semispan * math.tan(sweep_angle)

        # Upper and lower surf areas (Trapezoidal fins)
        Surf_areas = (root_chord + tip_chord)*0.5*semispan

        # Considering squase shaped profile
        Total_Fins_Volume = fins_thickness * (Surf_areas) * fins_number

        return Total_Fins_Volume
    except:
        print('### No Fins defined ###')


def fins_mass(rocket):                                                  # Rocket object

    try:                                                                # Checks if the fins are defined
        # Rocket number of fins
        Fins_Number = rocket.fins_number

        Total_Fins_Volume = fins_volume(
            rocket)                         # Fins volume

        Fins_Density = 0.00185  # 0.00122 #g/mm3                         # Get from firebase

        # Fins Mass

        Fins_Mass = Total_Fins_Volume * Fins_Density/1000               # Fins Mass in Kg
        return Fins_Mass
    except:
        print('### No fins defined ###')


def area(rocket):

    # Variables
    Bodytube_Lenght = rocket.body_length
    Bodytube_Diammeter = rocket.diammeter
    Nosecone_type = rocket.nosecone_type
    LD_Ratio = rocket.LD_ratio_nose
    Fins_Number = rocket.fins_number
    Boattail = rocket.boattail_present
    Fins_Thicknees = rocket.fins_thickness

    # odytube

    Wet_Area_Body = math.pi * Bodytube_Diammeter * \
        Bodytube_Lenght                 # Bodytube wet area (cylinder)

    # Fins

    root_chord = rocket.fins_dimensions[0]
    tip_chord = rocket.fins_dimensions[1]
    semispan = rocket.fins_dimensions[2]
    # to radians
    sweep_angle = rocket.fins_dimensions[3] * math.pi/180
    # Distace between root leading edge to tip LE
    sweep_distance = semispan * math.tan(sweep_angle)

    # Upper and lower
    Surfaces = (root_chord+tip_chord) * semispan

    # Leading edge surface
    Leading_edge = semispan/math.cos(sweep_angle) * Fins_Thicknees

    # Tip Surface
    Tip = Fins_Thicknees * tip_chord

    Trailing_edge = (math.sqrt(((sweep_distance + tip_chord - root_chord)**2)
                               + (semispan**2)) * Fins_Thicknees)                            # Trailing Edge Surface

    Wet_Area_Fins = (Surfaces + Leading_edge + Tip +
                     Trailing_edge) * Fins_Number   # Total Wet area

    # osecone

    # Calculating nosecone lenght based on LD ratio
    Nosecone_Length = rocket.nosecone_length
    # Calculating based on the shape
    if Nosecone_type == 'cone':                                                          # Cone
        Wet_Area_Nosecone = (math.sqrt((Nosecone_Length**2) +
                                       ((0.5*Bodytube_Diammeter)**2)) * math.pi *
                             (0.5*Bodytube_Diammeter))

    elif Nosecone_type == 'ogive':                                                        # ogive
        pass
    elif Nosecone_type == 'parabolic':                                                        # Parabolic
        pass
    elif Nosecone_type == 'elliptical':                                                        # Ellipsoid
        e = math.sqrt(1 - (((Bodytube_Diammeter/2)**2)/(Nosecone_Length**2)))

        Wet_Area_Nosecone = (math.pi * ((Bodytube_Diammeter/2)**2) *                 # Area molhada nosecone
                             (1 + ((Nosecone_Length/((Bodytube_Diammeter/2)*e))
                                   * math.asin(e))))

    # If there's a boattail
    if Boattail is not None:
        LD_Ratio_Boattail = rocket.LD_ratio_boattail
        # Upper Diammeter (Attached to the bodytube)
        D1_Boattail = Bodytube_Diammeter
        # Boattail Lenght based on LD ratio
        Boattail_Length = LD_Ratio_Boattail * Bodytube_Diammeter
        # Lower Diammeter
        D2_Boattail = rocket.boattail_Aft_Diammeter

        geratix = (math.sqrt((((D1_Boattail - D2_Boattail)*0.5)**2)                   # Cone generatrix
                             + (Boattail_Length**2)))

        # Boattail Wet area (only outer surface)
        Wet_Area_Boattail = math.pi * geratix * \
            ((D1_Boattail + D2_Boattail)*0.5)
    else:
        # if Boattail = False Area == 0
        Wet_Area_Boattail = 0

    Wet_Area_Rocket = (Wet_Area_Nosecone + Wet_Area_Body                            # Total wet area of the rocket
                       + Wet_Area_Fins + Wet_Area_Boattail)
    Wet_data = numpy.array([Wet_Area_Nosecone,
                            Wet_Area_Body, Wet_Area_Fins,
                            Wet_Area_Boattail,
                            Wet_Area_Rocket])                                          # Separate Wet areas

    return Wet_data


def rocket_drag(rocket, atmo_data, altitude, v):

    # Maioria dos cálculos foram baseados pela documentação do OpenRocket na seção de arrasto http://servidor.demec.ufpr.br/CFD/bibliografia/aerodinamica/Niskanen_2013.pdf !!!!!!!!!

    # Pi
    pi = math.pi

    # Troposfera (wikipedia ISA em pt [tinha um erro nas formulas, o expoente faltava M]) + https://www.digitaldutch.com/atmoscalc/US_Standard_Atmosphere_1976.pdf +https://pt.wikipedia.org/wiki/Atmosfera_padrão_internacional

    # Temperatura em uma altura x [K]
    Tx = atmo_data[-2]+(atmo_data[2])*altitude
    # Divisão da Tx
    Tdiv = Tx/atmo_data[-2]
    Px = (atmo_data[-1]*(math.pow((Tdiv), (-atmo_data[0]
                                           * atmo_data[1]/(atmo_data[3]*atmo_data[2])))))                              # Pressão em uma altura x [N/m^2]
    rhox = (atmo_data[-3]*math.pow((Tdiv), -(atmo_data[0]
                                             * atmo_data[1]/(atmo_data[3]*atmo_data[2])+1)))                           # Densidade em uma altura x [kg/m^3]

    Wet_data = area(rocket)

    # Geometrias
    # Espessura da aleta [m]  (Maior espessura pensando em uma aleta de secao quadrada)
    fin_t = rocket.fins_thickness/1e3
    # Diametro do foguete [m]
    body_diam = rocket.diammeter/1e3
    # Tamanho do bodytubo [m]
    body_tubeL = rocket.body_length/1e3
    # Area molhada da coifa [m^2]
    coifa_area = Wet_data[0]/1e6
    # Area molhada do bodytube [m^2]
    body_area = Wet_data[1]/1e6
    # Raiz dos fins [m]
    fin_root = rocket.fins_dimensions[0]/1e3
    # Corda da ponta dos fins [m]
    fin_tipc = rocket.fins_dimensions[1]/1e3
    # Angulo de enflechamento [graus]
    fin_sweep = rocket.fins_dimensions[3]
    # Taper ratio [ ]
    taperatio = fin_tipc/fin_root
    # Mean aero chord [ ]
    fin_mac = fin_root*2/3*((1+taperatio+math.pow(taperatio, 2))/(1+taperatio))
    # Altura dos fins [m]
    fin_h = rocket.fins_dimensions[2]/1e3
    # Area total dos fins (wet) [m^2]
    fin_area = Wet_data[2]/1e6
    # Secção da aleta; 1=quadrada []
    fin_corte_fin = 1
    # Area do Boattail [m^2]
    boat_area = Wet_data[3]/1e6
    # Comprimento Boattail [m]
    if rocket.boattail_present is not None:
        boat_L = rocket.LD_ratio_boattail * rocket.diammeter/1e3
        # Menor raio do Boattail [m]
        boat_aft = rocket.boattail_Aft_Diammeter/1e3
    else:
        boat_L = 0
        boat_aft = 0

    # Tamanho total do foguete [m]
    body_L_tot = (rocket.body_length + rocket.nosecone_length + boat_L)/1e3

    # sutherland law
    # Sutherland temperature
    suth = 110.56
    # Viscosidade dinamica
    visc = 1.7894e-5*math.pow((Tx/273.11), (3/2))*((273.11+suth)/(Tx+suth))
    # Viscosidade cinematica
    kin = visc/rhox

    # Variavei/calculo delas
    # Velocidade Local do Som
    sound = 331.3+0.606*(Tx-273.15)
    # Traduz velocidad para mach
    mach = v/sound
    # Rugosidade superficial [m]       !!!!!!!!! MUDAR PRO DE VERDADE !!!!!!!!!!
    Rs = 20e-6
    # Reynolds                       (correção para compresibilidade existe?)/?????????????????????
    Re = body_L_tot*v/kin
    # Reynolds critico
    Rcrit = 51*math.pow((Rs/body_L_tot), -1.039)

    # Cálculo do Coeficiente de Fricção

    # (Eq. 3.81) Primeira parte openrocket "barrowmancalculator.java"
    if Re <= 1e4:
        Cf = 1.48e-2

    elif Re <= Rcrit:
        # (Eq. 3.81) Segunda parte openrocket
        Cf = 1/math.pow(1.5*math.log(Re)-5.6, 2)

    else:
        # (Eq. 3.81) Terceira parte openrocket
        Cf = 0.032*math.pow(Rs/body_L_tot, 0.2)

    # Correção de Mach
    # (Eq. 3.82) }
    c1 = (1-(0.1)*math.pow(mach, 2))
    # (Eq. 3.83) }Openrocket
    c2 = 1/math.pow(1+0.15*math.pow(mach, 2), 0.58)
    # (Eq. 3.84) }
    c3 = 1/(1+0.18*math.pow(mach, 2))

    # Este "if" inteiro foi retirado direto das linhas do programa Openrocket (...)
    if mach < 0.9:
        # ele interpola os valores entre mach 0.9 e 1.1 para tentar pegar a transição (...)
        Cfc = Cf * c1
    #                                                                               # imagino que seja meio que "achismo".
    elif mach < 1.1 and Re <= Rcrit:
        Cfc = Cf * c1 * (mach - 0.9) / 0.2 + c2 * (1.1 - mach) / 0.2

    elif mach < 1.1:
        Cfc = Cf * (c1 * (mach - 0.9) / 0.2 + c3 * (1.1 - mach) / 0.2)

    elif mach > 1.1 and Re <= Rcrit:
        Cfc = Cf * c2

    else:
        Cfc = Cf * c3

    Cfc_mult = (((1+1/(2*(body_L_tot/body_diam)))*(body_area+coifa_area+boat_area) +  # }
                (1+2*fin_t/fin_mac) * (fin_area))/pi/math.pow(body_diam/2, 2))        # }(Eq. 3.85) Openrocket
    # }
    cd_fric = Cfc*Cfc_mult

    # Cálculo do Coeficiente de Body Fore Pressure Drag
    # (SOMENTE DE NOSECONE ELIPTICO POR ENQUANTO)

    if mach < 1:
        # (Eq. B.1) Stagnation pressure
        stag = 1+math.pow(mach, 2)/4+math.pow(mach, 4)/40
    else:
        stag = (1.84-0.76/math.pow(mach, 2) +
                0.166/math.pow(mach, 4)+0.035/math.pow(mach, 6))                         # (Eq. B.1) Stagnation pressure

    # Openrocket Docu
    cd_stag = 0.85*stag

    # Aproximação criada pela equipe (PURO ACHISMO) para calculo de arrasto de forma do nosecone
    if mach < 1.2:
        a = 0.0008357
        b = 4.063

        cd_press = a*math.exp(mach*b)

    else:
        fit_elips = interp1d(numpy.array([1.2, 1.25, 1.3, 1.4, 1.6, 2.0, 2.4]),
                             numpy.array([0.110, 0.128, 0.140, 0.148, 0.152, 0.159, 0.162]))  # Aproximação recriada pelo estudo: https://ntrs.nasa.gov/api/citations/19630004995/downloads/19630004995.pdf pag.16

        cd_press = cd_stag*math.pow((fit_elips(mach)/cd_stag),
                                    (math.log((body_L_tot/body_diam+1), 4)))           # somente supersonico

    # Cálculo do arrasto de pressão das Aletas (Seção 3.4.4)

    if fin_corte_fin == 1:                                                        # Secção Quadrada
        # Arrasto 90graus em relacao ao corpo (Eq. 3.90)
        cd_perp = cd_stag

        if mach < 1:
            #
            cd_base = 0.12+0.13*math.pow(mach, 2)
        # (Eq. 3.94) + (Eq. 3.92)
        else:
            cd_base = 0.25/mach                                                       #

    # Secção Arredondada/Aerofolio pelo jeito
    else:
        # Arrasto 90graus em relacao ao corpo
        if mach < 0.9:
            cd_perp = math.pow((1-math.pow(mach, 2)), (-0.417)
                               )-1                       #
        # (Eq. 3.89)
        elif mach < 1:
            #
            cd_perp = 1-1.785*(mach-0.9)
        else:                                                                       #
            cd_perp = 1.214-0.502 / \
                math.pow(mach, 2)+0.1095/math.pow(mach, 4)            #

        # "For fins with rounded edges the value is taken as half of the total base drag"
        if mach < 1:
            cd_base = (0.12+0.13*math.pow(mach, 2)) / \
                2                                  #
        else:                                                                       #
            cd_base = 0.25/mach/2                                                     #

    # (Figura 3.2(a))
    gamma_l = fin_sweep

    # (Eq. 3.91)
    cd_le = cd_perp*math.pow(math.cos(math.radians(gamma_l)), 2)

    # (Eq. 3.93) "The reference area is the fin frontal area N · ts"
    cd_press_fin = (cd_le+cd_base)*fin_t*fin_h/(pi*math.pow(body_diam, 2)/4)*4

    # Cálculo arrasto base body  (Seção 3.4.5)

    # (Eq. 3.94) Igual a usada na seção anterior
    if mach < 1 and fin_corte_fin != 1:
        cd_base = (0.12+0.13*math.pow(mach, 2))
    elif fin_corte_fin != 1:
        cd_base = 0.25/mach

    # Por enquanto n lembro !!!!!!! (Achar cd base com boattail pra usar la em baixo? mudar pra baixo seria bom?)
    if boat_L != 0:
        cd_base = cd_base * math.pow(boat_aft, 2)/math.pow(body_diam, 2)

    # Cáclulo de arrasto pressão boattail

    if boat_L != 0:
        # Caso exista Boattail, o gamma (Eq. 3.88)
        gamma = boat_L/(body_diam-boat_aft)

        if gamma < 1:
            gamma = 1

        elif gamma < 3:
            gamma = (3-gamma)/2

        else:
            gamma = 0

        # Caso       exista Boattail, o gamma (Eq. 3.88)
        boat_press = cd_base*gamma*math.pow(boat_aft/body_diam, 2)
    else:
        # Caso não exista Boattail, o arrasto é zero
        boat_press = 0

    # Soma de Todos os componentes de Arrasto
    cd_tot = cd_fric+cd_press+cd_press_fin+cd_base+boat_press
    # Multiplica-se pelo "rhox" para simplificar o retorno da função
    cd_tot_multrho = cd_tot*rhox

    cd = numpy.array([round(cd_tot_multrho, 5), round(cd_fric, 5), round(
        cd_press, 5), round(cd_press_fin, 5), round(cd_base, 5), round(cd_tot, 5)])

    # Retorno do "Cd total multiplicado pelo Rho"
    return cd


def StaticMarginSubsonic(rocket, aoa):
    # If Mach <= 0.3 Thsi funcion gives the SM
    # For more info please check 'Technical Report 33' And 'WInd Instability - What Barrowman Left out'
    # Variables
    # Bodytube Lenght
    Bodytube_Lenght = rocket.body_length
    # Bodytube Diammeter
    Bodytube_Diammeter = rocket.diammeter
    # Cone / ogive / parabolic / ellipsoid
    Nosecone_type = rocket.nosecone_type
    # Nosecone LD Ratio
    LD_Ratio_Nosecone = rocket.LD_ratio_nose
    # Fins Number
    Fins_Number = rocket.fins_number
    # Is there a boattail?
    Boattail = rocket.boattail_present
    # Derised Angle of attack
    Angle_of_attack = aoa
    # Calculates the fins mass
    # Fins Mass
    fin_mass = fins_mass(rocket)

    Custom_part = rocket.custom_part

    # Pi
    pi = math.pi
    Angle_of_attack = Angle_of_attack * pi / \
        180                                        # to radians
    # Barrowman Extended constant
    K_Barrowman_Ext = 1.1
    # Refference Diammeter
    Ref_Diammeter = Bodytube_Diammeter

    # Correct values
    Correction_Constant = 0.993

    ################ Nosecone ##########################

    # Normal Force Coeff at 0 AoA
    CnAlpha_Nose_Barrwoman = 2

    Nosecone_Length = LD_Ratio_Nosecone * \
        Bodytube_Diammeter                          # Nosecone Length

    if Nosecone_type == 'cone':   # Cone
        x_cp_nose_barrowman = (2/3)*Nosecone_Length
        # Extended Barrowman Method if AoA != 0
        if Angle_of_attack != 0:
            Area_Nose = 0.5 * Ref_Diammeter * Nosecone_Length
            x_plane_nose = (2/3) * Nosecone_Length

    elif Nosecone_type == 'ogive':  # Ogive
        x_cp_nose_barrowman = 0.466 * Nosecone_Length
        # Extended Barrowman Method if AoA != 0
        if Angle_of_attack != 0:
            Area_Nose = (2/3) * Ref_Diammeter * Nosecone_Length
            x_plane_nose = (5/8) * Nosecone_Length

    elif Nosecone_type == 'parabolic':  # Parabolic
        x_cp_nose_barrowman = 0.5 * Nosecone_Length
        # Extended Barrowman Method if AoA != 0
        if Angle_of_attack != 0:
            Area_Nose = (2/3) * Ref_Diammeter * Nosecone_Length
            x_plane_nose = (3/5) * Nosecone_Length

    elif Nosecone_type == 'elliptical':  # Elliptical
        x_cp_nose_barrowman = (1/3) * Nosecone_Length
        if Angle_of_attack != 0:
            Area_Nose = pi * (0.5*Ref_Diammeter) * Nosecone_Length * 0.5
            x_plane_nose = 0.625 * Nosecone_Length
    # Extended Barrowman Method if AoA != 0
    if Angle_of_attack != 0:
        CnAlpha_sq_nose = (4 * (K_Barrowman_Ext/pi) *
                           (Area_Nose/(Ref_Diammeter**2)) * Angle_of_attack)

     ############### Bodytube ##########################
    # Extended Barrowman Method if AoA != 0
    if Angle_of_attack != 0:
        x_plane_body = Nosecone_Length + (Bodytube_Lenght/2)

    # Planiform Area
        Area_Body = Bodytube_Lenght * Bodytube_Diammeter

        CnAlpha_sq_Body = (4 * (K_Barrowman_Ext/pi)
                           * (Area_Body/(Ref_Diammeter**2))
                           * Angle_of_attack)

     ################ Fins #############################

    root_chord = rocket.fins_dimensions[0]
    tip_chord = rocket.fins_dimensions[1]
    semispan = rocket.fins_dimensions[2]
    # to radians
    sweep_angle = rocket.fins_dimensions[3] * pi/180

    sweep_distance = semispan * math.tan(sweep_angle)

    # For calculation
    cateto = sweep_distance + (0.5*tip_chord) - (0.5*root_chord)
    # Distance from 1/2 Root chord to 1/2 Tip Chord
    l = math.sqrt((semispan**2)+(cateto**2))

    CnAlpha_Fins_Barrwoman = (((4*Fins_Number)
                              * ((semispan / Bodytube_Diammeter)**2))
                              / (1 + math.sqrt(1+(((2*l)
                                                   / (root_chord + tip_chord))**2))))  # Barrowman

    # Fins Interference Factor for n = 3 or 4

    K = 1 + ((0.5*Bodytube_Diammeter)/((0.5*Bodytube_Diammeter)+(semispan)))

    # Cn alpha fins with interference factor
    CnAlpha_Fins_Barrwoman = CnAlpha_Fins_Barrwoman * K

    Leading_Edge_Position = Bodytube_Lenght + Nosecone_Length - \
        root_chord                # Fins Leading Edge Position

    x_cp_fins_barrowman = (Leading_Edge_Position
                           + ((sweep_distance*(root_chord+(2*tip_chord)))
                               / (3*(root_chord+tip_chord))) + ((1/6)
                                                                * (root_chord+tip_chord-((root_chord*tip_chord)
                                                                                         / (root_chord+tip_chord)))))                                   # CP position fins

    ################ Boattail #####################

    # Only if Boattail = True
    if Boattail:
        if Custom_part:
            LD_Ratio_Boattail = rocket.LD_ratio_boattail
            D1_Boattail = Bodytube_Diammeter
            Boattail_Length = LD_Ratio_Boattail * Bodytube_Diammeter
            D2_Boattail = rocket.boattail_Aft_Diammeter
        else:
            LD_Ratio_Boattail = rocket.LD_ratio_boattail
            D1_Boattail = Bodytube_Diammeter
            Boattail_Length = LD_Ratio_Boattail * Bodytube_Diammeter
            D2_Boattail = rocket.boattail_Aft_Diammeter

            CnAlpha_Boattail_Barrwoman = (2*(((D2_Boattail/Bodytube_Diammeter)**2)
                                             - ((D1_Boattail/Bodytube_Diammeter)**2)))            # Boattail Cn at 0 AoA

            # Boattail location
            x_Boattail = Nosecone_Length + Bodytube_Lenght

            x_cp_boattail_barrowman = (x_Boattail + ((Boattail_Length/3) *
                                                     (1 + ((1-(D1_Boattail/D2_Boattail))
                                                           / (1 - ((D1_Boattail/D2_Boattail)**2))))))            # CP Location Boattail
            if Angle_of_attack != 0:
                Area_Boattail = 0.5 * (D1_Boattail + D2_Boattail) * \
                    Boattail_Length              # Planform Area

                x_plane_boattail = (x_Boattail + ((Boattail_Length/3)
                                    * ((D1_Boattail + (2*D2_Boattail))
                                    / (D1_Boattail + D2_Boattail))))                              # Area Centroid

                CnAlpha_sq_Boattail = (4 * (K_Barrowman_Ext/pi)
                                       * (Area_Boattail/(Ref_Diammeter**2))
                                       * Angle_of_attack)                                         # Extended Barrowman

    else:                                                                                # If Boattail = False
        CnAlpha_Boattail_Barrwoman = 0
        x_cp_boattail_barrowman = 0
        CnAlpha_sq_Boattail = 0
        x_plane_boattail = 0

    ################# CP Position #####################

    CnAlpha_rocket = (CnAlpha_Nose_Barrwoman
                      + CnAlpha_Fins_Barrwoman
                      + CnAlpha_Boattail_Barrwoman)                                       # Cn of the rocket at 0 AoA
    # Nose + Body + Boattail CG with prop + Fins (Time = 0)
    # Sets the Rocket CG to the Calculated
    Actual_CG = rocket.cg(0)
    # If AoA != 0 --> Barrowman Extended
    if Angle_of_attack != 0:
        Cn_rocket = (CnAlpha_rocket + CnAlpha_sq_nose
                     + CnAlpha_sq_Body + CnAlpha_sq_Boattail)                               # Cn at AoA != 0
        x_cp_rocket = (((CnAlpha_Nose_Barrwoman*x_cp_nose_barrowman)
                        + (CnAlpha_Fins_Barrwoman*x_cp_fins_barrowman)
                        + (CnAlpha_Boattail_Barrwoman * x_cp_boattail_barrowman)
                        + (CnAlpha_sq_nose * x_plane_nose) +
                        (CnAlpha_sq_Body*x_plane_body)
                        + (x_plane_boattail * CnAlpha_sq_Boattail))/(Cn_rocket))              # CP Location
    else:                                                                                 # Pure Barrowman
        x_cp_rocket = (((CnAlpha_Nose_Barrwoman*x_cp_nose_barrowman)
                        + (CnAlpha_Fins_Barrwoman*x_cp_fins_barrowman)
                        + (CnAlpha_Boattail_Barrwoman * x_cp_boattail_barrowman))
                       / (CnAlpha_rocket))                                                  # CP Location

    # Static Margin For MAch < 0.3
    StaticMargin_sub = (x_cp_rocket - Actual_CG)/Bodytube_Diammeter

    # return StaticMargin_sub
    return round(StaticMargin_sub, 2)


class rocket:
    def __init__(self, dry_mass, initial_cg, diammeter, body_length, boattail_present=False, custom_part=False):
        # Rocket dry mass (Kg)
        self.dry_mass = dry_mass
        # Bodytube + nosecone + boattail CG
        self.initial_cg = initial_cg
        # Maximum rocket diammeter
        self.diammeter = diammeter
        # bo dytube length (cylinder)
        self.body_length = body_length
        # Boattail boolean
        self.boattail_present = boattail_present
        self.custom_part = custom_part

        # Default value for Thrust time
        self.Total_burn_time = None
        # Default value for Fins array
        self.fins_dimensions = None
        self.LD_ratio_boattail = None

    def nosecone(self, nosecone_type, LD_ratio_nose):
        self.nosecone_type = nosecone_type
        self.LD_ratio_nose = LD_ratio_nose
        self.nosecone_length = LD_ratio_nose * self.diammeter

    def fins(self, fins_number, fins_thickness, fins_dimensions=None):                # Define Fins
        self.fins_number = fins_number                                        # Fins number
        # Fins thickness (rectangular shape)
        self.fins_thickness = fins_thickness
        # Dimensions array
        self.fins_dimensions = fins_dimensions

    def boattail(self, LD_ratio_boattail, boattail_Aft_Diammeter):
        # Adding boattail
        self.boattail_present = True
        self.LD_ratio_boattail = LD_ratio_boattail
        # After diammeter in mm
        self.boattail_Aft_Diammeter = boattail_Aft_Diammeter
        self.boattail_length = LD_ratio_boattail * \
            self.diammeter             # Boattail length in mm

    def custom_boattail(custom_part):
        # Adding boattail
        self.boattail_present = True
        self.custom_part = True
        self.LD_ratio_boattail = custom_part.LD_ratio_boattail
        # After diammeter in mm
        self.boattail_Aft_Diammeter = custom_part.boattail_Aft_Diammeter
        self.boattail_length = custom_part.LD_ratio_boattail * \
            self.diammeter             # Boattail length in mm

    def engine(self, thrust_curve=None, prop_mass=0):                     # Define engine

        # Thrust curve file name (str)
        self.thrust_curve = thrust_curve
        # propellant mass
        self.prop_mass = prop_mass

        # Sets the initial Thrust to 0
        Initial_Thrust = numpy.array([[0, 0]])
        Burn = numpy.loadtxt(self.thrust_curve,
                             comments=('H', 'I', 'J', ';'), skiprows=1)            # Loads the points, skips the comments and the title
        # Skips the first rows if the engine name starts with H I J

        self.Thrust_curve = numpy.append(
            Initial_Thrust, Burn, axis=0)          # Full Thrust Curve

        shape = self.Thrust_curve.shape
        burn_time = self.Thrust_curve[shape[0]-1, 0]

        # Total Burn time
        self.Total_burn_time = burn_time

    # Gets the thrust in a given time
    def Thrust(self, t):
        try:                                                                  # Checks if an engine is defined
            if t <= self.Total_burn_time:
                # Interpolation
                Engine_Thrust = interp1d(
                    self.Thrust_curve[:, 0], self.Thrust_curve[:, 1])
                return Engine_Thrust(t)
            else:
                # Thrust = 0 if prop = 0
                return 0
        except Exception:                                                     # If no egine has been defined
            print('\n### No engine defined ###\n')

    def mass(self, t):
        if self.Total_burn_time is None:                                      # No engines
            print('\n### No engine defined ###\n')
            if self.fins_dimensions is None:                                  # No engines nor fins
                print('\n### No fins defined ###\n')
                return self.dry_mass
            else:                                                             # No engines but there are fins
                return self.dry_mass + fins_mass(self)
        elif self.fins_dimensions is None:                                    # No fins
            print('\n### No fins defined ###\n')
            if self.Total_burn_time is None:                                  # No fins nor engines
                print('\n### DONT PANIC ###\n ### No engine defined ###\n')
            else:                                                             # No fins but there are engines
                if t < self.Total_burn_time:
                    mass_data = numpy.array([self.dry_mass
                                             + self.prop_mass,
                                            self.dry_mass])                  # Mass Array

                    mass_time = numpy.array(
                        [0, self.Total_burn_time])         # Time array
                    # Mass interpolation
                    Mass = interp1d(mass_time, mass_data)
                    return Mass(t)

                else:                                                         # if prop = 0 mass == Dry mass
                    return self.dry_mass

        # Real time mass (Linear function)
        else:

            if t < self.Total_burn_time:
                # Calculates fins mass (WORK IN PROGRESS)
                self.fins_mass = fins_mass(self)
                mass_data = numpy.array([self.dry_mass
                                         + self.prop_mass
                                         + self.fins_mass,
                                        self.dry_mass
                                         + self.fins_mass])                    # Mass Array

                mass_time = numpy.array(
                    [0, self.Total_burn_time])         # Time array
                # Mass interpolation
                Mass = interp1d(mass_time, mass_data)
                return Mass(t)

            else:                                                             # if prop = 0 mass == Dry mass
                return self.dry_mass

    def cg(self, t):                                                           # CG calculations
        if self.fins_dimensions is not None:                                  # Check if there are fins
            self.fins_cg = (self.nosecone_length +
                            self.body_length - (0.5
                                                * self.fins_dimensions[0]))                       # Fins CG at 1/2 chord

            Actual_cg = (((self.fins_cg * fins_mass(self))
                          + (self.initial_cg * (self.mass(t)
                                                - fins_mass(self))))
                         / (self.mass(t)))

            return round(Actual_cg, 2)
        # Initial CG with propellant if there are no fins
        else:
            print('\n### No fins or nosecone defined ###\n')
            return self.initial_cg

    def static_margin(self, mach, aoa=0):
        if self.fins_dimensions is not None:                                  # Checks if fins are defined
            if mach < 0.3:                                                    # Barrowman Method
                sm = StaticMarginSubsonic(self, aoa)
                return sm
            else:                                                             # Stability for mach > 0.3
                pass

        else:                                                                 # if fins are not defined
            print('\n### No fins defined ####\n')

    def drag(self, v, altitude=0):
        drag = rocket_drag(self, atmo_data, altitude, v)
        return drag


class custom_boattail:
    def __init__(self,  LD_ratio_boattail, boattail_Aft_Diammeter, cp_array):
        self.LD_ratio_boattail = LD_ratio_boattail
        self.boattail_Aft_Diammeter = boattail_Aft_Diammeter
        self.cp_array = cp_array
