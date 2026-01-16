Il codice per la simulazione di un tracciatore configurabile al silicio, per muoni, è stato strutturato in due script .py: funzioni.py contiene tutte le funzioni progettate ed utilizzate per la simulazione;
tracciatore.py è il codice effettivo da eseguire, in cui viene importato il primo script come modulo. In questo modo il codice "principale" risulta semplice e permette di essere modificato più agevolmente.
Per ogni funzione utilizzata è stata approntata una Docstring per eventuali chiarimenti.

All'esecuzione di tracciatore.py viene chiesto di inserire le variabili della simulazione, cioè il numero di strati del tracciatore, la distanza fra questi, la dimensione degli elementi di misura e l'energia del fascio di muoni. Il programma restituisce l'angolo medio di deviazione e relativa incertezza, come richiesto, e produce vari grafici informativi.

Per scelta di comodità, il numero di elementi del fascio di muoni è stato impostato a 10000; il valore è tuttavia modificabile entrando nello script tracciatore.py e intervenendo sulla variabile N.