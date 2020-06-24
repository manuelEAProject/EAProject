from scipy.interpolate import griddata
import numpy as np
import time
from matplotlib import pyplot


a=np.random.random((10, 33))


"""
var_range = 0.8
variation_start = [(1 - var_range + (var_range / (chromo_resolution / 2)) * gen_value) for gen_value in chromo[-3:-1:1]]

print (list(enumerate(chromo)))
chromo.sort()
print (list(enumerate(chromo)))

variation_start.append(55)
print (variation_start)
"""
"""
variation_start = 1 - var_range + (var_range / (chromo_resolution / 2)) * gen_value

print(variation_start)

def wechsleFarbe():
    lab1["bg"] = "#FFFF00"
    return
def gebeStr():
    lab1["text"]=eingabe.get()
    eingabe.delete(0, END)

#Main-Fenster öffnen
main = Tk()

#Label machen
lab1 = Label(main, text="Hallo!")
lab1.grid(row=0, column=0)

#Eingabe
eingabe = Entry(main, width=20)
eingabe.grid(row=0, column=1)
eingabe.insert(0, "Was geht")

#Knopf
b = Button(main, text="Knopf", command=wechsleFarbe)
b.grid(row=1, column=0)
c = Button(main, text="Übergeben", command=gebeStr)
c.grid(row=1, column=1)



main.mainloop()
"""