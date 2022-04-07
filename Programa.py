from tempfile import TemporaryDirectory
from Foguete import Foguete
import numpy as np
from runge import Runge

foguete = Foguete(4.97, 0.4, "MOTOR.txt")
foguete.plotEmpuxo()

foguete.interpolador()
    