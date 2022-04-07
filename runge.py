from Foguete import Foguete
import numpy as np
import matpltlib.pyplt as plt

def runge():
    def f(x,y):
        return #formula 
        v0 = 0.01
        Fa = 0.5 * cd * rho * A * v**2
        P = m * g


        # PVI
        a0 = a = (empuxo - P - Fa) / m
        dH = dH = ((a * tempo**2) / 2) + (v0 * tempo)
        v = ((v0**2) - 2 * a * dH)**(1/2)

        # Intervalo
        a = 0        # Inicio intervalo
        b = 1.8923   # Final do intervalo
        dt = 0.01    # Passo

        # Calculando n (numero de divisões)
        n = (b-a)/dt

        # Vetores que irão acumular os pontos gerados
        ai = [a0]
        dHi = [dH]
        vi = [v0]     

        # Iterando com o método
        for i in range(n):
            k1v = a(x0, v0, t0)*dt
            k1x = vn*dt
            k2v = a(x0 + 0.5*k1, v0 + 0.5*k1v, t0 + 0.5*dt)*dt
            k2x = (v0 + 0.5*k1v)*dt
            k3v = a(x0 + 0.5*k2, v0 + 0.5*k2v, t0 + 0.5*dt)*dt
            k3x = (v0 + 0.5*k2v)*dt
            k4v = a(x0 + k3x, v0 + k3v, t0 + dt)*dt
            k4x = (v0 + k3x)*dt


            v2 = v0 + (1/6)*(k1v + 2*k2v + 2*k3v + k4v)

            x2 = x0 + (1/6)*(k1x + 2*k2x + 2*k3x + k4x)


            k2 = f(x0 + h, y0 + h*k1)
            yk = y0 + 0.5*(k1 + k2)*h
            x0 = x0 + h
            y0 = yk
            print(f"({round(x0,2)} , {round(y0,2)})")

            # Adicionando os pontos nos vetores
            xi.append(round(x0,2))
            yi.append(round(yk.2))


    print()
    print(xi)
    print(yi)

    
    plt.plot(xi,yi,"o")
    plt.show()