import numpy as np
import funzioni as f
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import warnings
warnings.filterwarnings('ignore')

# Il tuo codice qui


#condizioni iniziali


n1 = input('Digitare il numero di strati al silicio desiderati: ') # numero di piani
n = int(n1)
E1 = input("Digitare l'energia desiderata del fascio di muoni, nell intervallo 10-10000 MeV: ") # energia del fascio 10-10000 MeV
E = float(E1)
d1 = input('Digitare la distanza fra gli strati desiderata, in cm: ')  #distanza tra i piani cm
d = float(d1)
s1 = input('Digitare la larghezza desiderata del pitch (elementi discreti di misura), in um: ') #dimensione dei pitch 0.005-0.04 cm
s = float(s1) /10000
N = 10000 #elementi del fascio

#condizioni iniziali
dev_x = np.empty(0)
err_x = np.empty(0)
dev_y = np.empty(0)
err_y = np.empty(0)

traiett = np.zeros((n+1, 3))
serr = s / np.sqrt(12)


#simulazione di N traiettorie di particelle
for j in range(N):

    x0 = 0
    y0 = 0
    z0 = 0
    
    r0 = np.array([x0, y0, z0])
    
    ox0 = 0
    oy0 = 0
    
    traiettoria = np.array([ r0 ])
    
    for i in range(n) :
        
        ox, dx, oy, dy = f.cscatter(ox0, oy0, E)
        z, x, y =  f.spostamento( ox, dx, oy, dy, d)
        z1 = z0 + z
        x1 = x0 + x
        y1 = y0 + y
        r1 = np.array([ x1, y1, z1 ])
        x2 = f.mstrip(x1, s)
        y2 = f.mstrip(y1, s)
        r = np.array([x2, y2, z1])
        traiettoria = np.vstack([traiettoria, r ])
        z0 = z1
        x0 = x1
        y0 = y1
        ox0 += ox
        oy0 += oy

    a, b, e = f.regrx(traiettoria, serr)
    c, l, g = f.regry(traiettoria, serr)
    dev_x = np.append(dev_x, b)
    err_x = np.append(err_x, e)
    dev_y = np.append(dev_y, l)
    err_y = np.append(err_y, g)
    traiett = traiettoria
    

print('====================================================================== ')
#calcolo angolo medio di uscita e relativo errore
dev, err = f.angolo(dev_x, err_x, dev_y, err_y)
dev1 = 'Angolo medio di deviazione, in radianti: {:.2e}'.format(dev)
print(dev1)
err1 = "Deviazione standard dell'angolo di deviazione, in radianti {:.2e}".format(err)
print(err1)
print('====================================================================== ')
print(' ')
print(' ')


#sezione grafici
f.grafico(traiett, n, d )
f.fitxy(traiett, s)
b1, b2 = f.infobins(dev_x, N)
f.istxy2(dev_x, dev_y, b1)
f.heatmap(dev_x, dev_y, b1)
f.ist3d(dev_x, dev_y, 2 * b1 , b2)






