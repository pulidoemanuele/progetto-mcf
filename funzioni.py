import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import linalg

h = 2 * 10**(-2) #cm
mc2 = 105.7 #MeV


def cscatter(ox1, oy1, E ) :
    """
    Genera angoli di deviazione (Ox, Oy) e spostamenti (dx, dy) a causa dello scattering coulombiano secondo le formule del PDG.
    ----------------------------------------------------------------------------------------------------------------------------
    Parametri:
    ox1 = float, angolo nel piano zx rispetto all'asse z con cui la particella entra nello strato di rivelazione;
    oy1 = float, angolo nel piano zy rispetto all'asse z con cui la particella entra nello strato di rivelazione;
    E = float, energia cinetica relativistica della particella.
    """
    
    pc = np.sqrt(( E**2) + 2*E*mc2 )
    b = pc / (E + mc2)

    
    ax = np.tan(ox1)

    ay = np.tan(oy1)

    h1= h * np.sqrt((ax)**2 + (ay)**2 + 1 )
    O = 13.6 / (b * pc ) * np.sqrt( h1 / 9.36 ) * (1 + 0.038*np.log(h1/(9.36 * (b**2))))
    z = np.random.normal(0, 1, 4)
    dx = z[0]*h1*O / (np.sqrt(12)) + z[1]*h1*O / 2
    Ox = z[1] * O
    dy = z[2]*h1*O / (np.sqrt(12)) + z[3]*h1*O / 2
    Oy = z[3] * O

    return Ox, dx, Oy, dy


def spostamento(ox, dx, oy, dy, d) :
    """
    Calcola e restituisce di quanto la particella si è spostata in ogni direzione fra uno strato ed il successivo.
    --------------------------------------------------------------------------------------------------------------
    Parametri:
    ox = float, angolo nel piano zx rispetto all'asse z;
    dx = float, spostamento che la particella ha subito lungo x nell'attraversare uno strato del tracciatore;
    oy = float, angolo nel piano zy rispetto all'asse z;
    dy = float, spostamento che la particella ha subito lungo y nell'attraversare uno strato del tracciatore;
    """
    deltaz = d
    deltax = dx + ((d-h) * np.tan(ox))
    deltay = dy + ((d-h) * np.tan(oy))    
    return deltaz, deltax, deltay


def mstrip(a, b):
    """
    Dato che gli elementi di misura sono discreti, arrotonda la posizione calcolata a multipli della larghezza del pitch.
    ---------------------------------------------------------------------------------------------------------------------
    Parametri:
    a = float, posizione x o y calcolata;
    b = float, larghezza del pitch.
    """
    return ( round(a / b) * b )






    
def retta(matrice, n, d):
    """
    Restituisce una retta, matrice nx3 dove le righe sono punti, che rappresenta una regressione lineare, molto approssimativa, della traiettoria della particella;
    Questa funzione viene usata esclusivamente per il grafico, non è di utilità per ottenere risultati numerici data la sua imprecizione.
    ---------------------------------------------------------------------------------------------------------------------------------------------------------------
    Parametri:
    matrice = matrice nx3 di float, rappresenta la traiettoria ricostruita dal traccia;
    n = int, numero di piani, quindi numero di punti registrati dal tracciatore;
    d = float, distanza fra i piani del tracciatore.
    """
    baricentro = np.mean(matrice, axis = 0)
    norma = np.linalg.norm(baricentro)
    direzione = baricentro / norma
    zmax = n * d
    tmax = zmax/ direzione[2]
    t= np.linspace(0, tmax, n * 10)
    puntiretta = np.empty((0, 3))
    for i in range( len(t)):
        punto = direzione * t[i]
        puntiretta = np.vstack([puntiretta, punto])
    return puntiretta





def grafico(matrice, m, d):
    """
    Produce un grafico tridimensionale che mostra la traiettoria della particella misurata dal tracciatore, la traiettoria approsimata dal fit e i vari piani di rivelazione.
    -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Parametri:
    matrice = matrice nx3 di float, rappresenta la traiettoria ricostruita dal traccia;
    n = int, numero di piani, quindi numero di punti registrati dal tracciatore;
    d = float, distanza fra i piani del tracciatore.
    """
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Estrai coordinate
    x_points = matrice[:, 0]
    y_points = matrice[:, 1]
    z_points = matrice[:, 2]
    
    # Calcola limiti per stessa scala x/y
    x_min, x_max = x_points.min(), x_points.max()
    y_min, y_max = y_points.min(), y_points.max()
    
    x_center = (x_min + x_max) / 2
    y_center = (y_min + y_max) / 2
    
    x_range_actual = x_max - x_min
    y_range_actual = y_max - y_min
    
    # Usa il range massimo per mantenere la stessa scala
    max_range_actual = max(x_range_actual, y_range_actual)
    
    # Se max_range_actual è molto piccolo (dell'ordine di 10^-6), 
    # mantieni i limiti reali dei dati con un piccolo margine
    if max_range_actual < 1e-5:  # Se i dati sono molto piccoli
        # Usa un margine percentuale invece di assoluto
        margin_factor = 0.1
        x_lim = [x_min - x_range_actual * margin_factor, 
                 x_max + x_range_actual * margin_factor]
        y_lim = [y_min - y_range_actual * margin_factor, 
                 y_max + y_range_actual * margin_factor]
    else:
        # Per dati più grandi, usa l'approccio originale
        x_lim = [x_center - max_range_actual/2 - max_range_actual*0.1, 
                 x_center + max_range_actual/2 + max_range_actual*0.1]
        y_lim = [y_center - max_range_actual/2 - max_range_actual*0.1, 
                 y_center + max_range_actual/2 + max_range_actual*0.1]
    
    # Crea griglia per piani
    X, Y = np.meshgrid(np.linspace(x_lim[0], x_lim[1], 50),
                       np.linspace(y_lim[0], y_lim[1], 50))
    
    # Disegna piani
    for i in range(m + 1):
        z = i * d
        Z = np.full_like(X, z)
        
        ax.plot_surface(X, Y, Z, alpha=0.1, color='lightblue', 
                       edgecolor='gray', linewidth=0.1)
    
    # Punti e linee
    ax.plot(x_points, y_points, z_points, 'b-', alpha=0.6)
    ax.scatter(x_points, y_points, z_points, c='red', s=30)
    
    # Fit lineare
    retta1 = retta(matrice, m, d)
    ax.plot(retta1[:,0], retta1[:,1], retta1[:,2], 
                'r-', linewidth=2)
    
    # Imposta limiti assi
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.set_zlim([0, m*d])  
    
    # Stessa scala x/y
    ax.set_box_aspect([1, 1, 0.5])
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'traiettoria ricostruita + fit 3d')
    
    plt.tight_layout()
    plt.show()



    




def regrx(matrice, s):
    """
    Tramite il metodo dei minimi quadrati restituisce la retta che meglio approssima la traiettoria nel piano zx, il coefficiente angolare di tale retta e l'errore associato a questo.
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Parametri:
    matrice = matrice nx3 di float, rappresenta la traiettoria ricostruita dal traccia;
    s = risoluzione spaziale.
    """
    delta = np.sum(matrice[:,2]**2)
    m = np.sum(matrice[:,2] * matrice[:,0]) / delta
    zx = m * matrice[:,2]
    errorex = s * np.sqrt(len(matrice[:,0]) / delta )
    
    return zx, m, errorex



def regry(matrice, s):
    """
    Tramite il metodo dei minimi quadrati restituisce la retta che meglio approssima la traiettoria nel piano zx, il coefficiente angolare di tale retta e l'errore associato a questo.
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Parametri:
    matrice = matrice nx3 di float, rappresenta la traiettoria ricostruita dal traccia;
    s = risoluzione spaziale.
    """
    delta = np.sum(matrice[:,2]**2)
    m = np.sum(matrice[:,2] * matrice[:,1]) / delta
    zy = m * matrice[:,2]
    errorey = s * np.sqrt(len(matrice[:,0]) / delta )
    return zy, m, errorey






def fitxy(matrice, s, titolo1="fit lineare per x", titolo2="fit lineare per y"):
    """
    Genera due grafici affiancati che mostrano la traiettoria registrata nel piani zx e zy con rispettivi fit.
    ----------------------------------------------------------------------------------------------------------
    Parametri:
    matrice = matrice nx3 di float, rappresenta la traiettoria ricostruita dal traccia;
    s = risoluzione spaziale.
    """
    a, b, e = regrx(matrice, s)
    c, d, f = regry(matrice, s)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Primo grafico
    ax1.plot(matrice[:,2], a, 'b-', linewidth=2, markersize=4)
    ax1.scatter( matrice[:,2], matrice[:,0], color='royalblue')
    ax1.set_xlabel('Z')
    ax1.set_ylabel('X')
    ax1.set_title('fit x')
    ax1.grid(True, alpha=0.3)
    
    # Secondo grafico
    ax2.plot(matrice[:,2], c, 'r-', linewidth=2, markersize=4)
    ax2.scatter( matrice[:,2], matrice[:,1], color='royalblue')
    ax2.set_xlabel('Z')
    ax2.set_ylabel('Y')
    ax2.set_title('fit y')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()



def pesata(x, xerr):
    """
    Restituisce media pesata e relativo errore di un array di dati di cui ognuno ha differente incertezza.
    ---------------------------------------------------------------
    Parametri:
    x = array di float
    xerr = incertezze dei dati.
    """
    pesi = 1 / xerr**2
    mpesata = np.average(x, weights = pesi)
    err = np.sqrt(1 / np.sum(pesi))
    return mpesata, err



def angolo(x, xerr, y, yerr):
    """
    Restituisce l'angolo rispetto all'asse z, con errore associato, con cui la particella esce dal tracciatore.
    ------------------------------------------------------------------------------------
    Parametri:
    x, y = array di float, deviazioni nei rispettivi piani della particella;
    xerr, yerr = array di float, errori associati ad ogni deviazione.
    """
    a, b = pesata(x, xerr)
    c, d = pesata(y, yerr)
    o  = np.sign(a) * np.sign( c ) * np.sqrt( a**2 + c**2 )
    err = np.sqrt( b**2 + d**2 )
    return o, err




def infobins(dati, n):
    """
    Restituisce numero e larghezza ottimale dei bins di un istogramma, a seconda dei dati che si vogliono utilizzare, usando le regole di Sturges e Scott.
    ------------------------------------------------------------------------------------------------------------------------------------------------------
    Parametri:
    dati = array di float, insieme di dati con cui si vuole realizzare l'istogramma;
    n = int, dimensione dell'array di dati.
    
    """
    sigma = np.std(dati, ddof=1)
    a = 3.49 * sigma * (n ** (-1/3))
    b = int(np.ceil(1 + np.log2(n)))
    return b, a




def heatmap(x, y, b1, cmap='viridis', title=None):
    """
    Produce una heatmap, paragonabile a un istogramma di due variabili, con scatterplot per visualizzare la distribuzione dei valori.
    ---------------------------------------------------------------------------------------------------------------------------------
    Parametri:
    x, y = float, coppia di dati su cui costruire i grafici;
    b1 = numero di bins;
    cmap = string, gradiente di colori da utilizzare, opzionale;
    title = string, titolo del grafico, opzionale.
    """
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))       
    
    # Calcola istogramma
    hist, xedges, yedges = np.histogram2d(x, y, bins= b1  )
    
    # Calcola informazioni sui bins
    if len(xedges) > 1:
        actual_bin_width_x = xedges[1] - xedges[0]
        bins_x_count = len(xedges) - 1
    else:
        actual_bin_width_x = 0
        bins_x_count = 0
        
    if len(yedges) > 1:
        actual_bin_width_y = yedges[1] - yedges[0]
        bins_y_count = len(yedges) - 1
    else:
        actual_bin_width_y = 0
        bins_y_count = 0
    
    # Heatmap
    im = ax1.imshow(hist.T, origin='lower',  extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cmap, aspect='auto')
    
    cbar = plt.colorbar(im, ax=ax1, label='Frequenza')
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    
    # Titolo
    if title is None:
        title = f'Heatmap 2D'
    ax1.set_title(title)
    
    # Informazioni sui bins nel grafico
    bins_info_text = (f'Bins: {bins_x_count}×{bins_y_count}\n'
                     f'Width: ({actual_bin_width_x:.3f}, {actual_bin_width_y:.3f})')
    
    ax1.text(0.02, 0.98, bins_info_text, transform=ax1.transAxes,
            fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
    
    # Scatter plot
    ax2.scatter(x, y, alpha=0.3, s=10, c='blue', edgecolors='none')
    ax2.set_xlabel('angolo x')
    ax2.set_ylabel('angolo y')
    ax2.set_title('distribuzione delle deviazioni')
    ax2.grid(True, alpha=0.3)
    
    
    plt.tight_layout()
    plt.show()

    

def ist3d(x, y, b1, b2):
    """
    Produce un istogramma tridimensionale per coppie di valori.
    -----------------------------------------------------------
    Parametri:
    x, y = float, coppia di dati su cui costruire i grafici;
    b1 = numero di bins;
    b2 = larghezza dei bins.
    """
    
    # Crea figura 3D
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Crea l'istogramma 2D
    hist, xedges, yedges = np.histogram2d(x, y, bins=b1)
    
    # Calcola le larghezze dei bin
    x_bin_widths = np.diff(xedges)
    y_bin_widths = np.diff(yedges)
    
    # Crea meshgrid delle posizioni dei bin (centri dei bin)
    xpos, ypos = np.meshgrid(xedges[:-1] + x_bin_widths/2, 
                             yedges[:-1] + y_bin_widths/2)
    
    # Appiattisci le matrici
    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros_like(xpos)
    
    dx = dy = b2 
    
    # Appiattisci l'istogramma
    dz = hist.flatten()
    
    # Plot dell'istogramma 3D
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, zsort='average',
             color='lightblue', edgecolor='black', alpha=0.8)

    # Etichette
    ax.set_xlabel('angolo x')
    ax.set_ylabel('angolo y')
    ax.set_zlabel('Frequenza')
    ax.set_title('Istogramma 3D di coppie (Ox,Oy)')
    
    # Imposta limiti degli assi
    ax.set_xlim(xedges[0], xedges[-1])
    ax.set_ylim(yedges[0], yedges[-1])
    
    plt.show()




def istxy2(x, y, n):
    """
    Produce due istogrammi affiancati.
    ----------------------------------
    Parametri:
    x, y = array di float, dati con cui realizzare il corrispettivo istogramma;
    n = int, dimensione degli array di dati.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    axes[0].hist(x, bins=n, 
                 color='skyblue', alpha=0.7, edgecolor='black')
    axes[0].set_title('Distribuzione Ox')
    
    axes[1].hist(y, bins=n, 
                 color='salmon', alpha=0.7, edgecolor='black')
    axes[1].set_title('Distribuzione Oy')
    
    plt.tight_layout()
    plt.show()

