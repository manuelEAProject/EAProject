from scipy.interpolate import griddata
import numpy as np









def func(x, y):
     return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2
#on a grid in [0, 1]x[0, 1]


grid_x, grid_y = np.mgrid[0:1:100j, 0:1:200j]
#but we only know its values at 1000 data points:


l_list = [[0.12,0.23],[0.32,0.41],[0.12,0.41],[0.23,0.41],[0.45,0.41],[0.21,0.41],[0.89,0.41],[0.78,0.41],
          [0.14,0.22],[0.98,0.98]]
points = np.asarray(l_list, dtype=np.complex64)

#points = np.random.rand(1000, 2)
values = np.array([0.11,0.22,0.34,0.98,0.11,0.22,0.34,0.98,0.34,0.34], dtype=np.complex64)


#Thiss can be done with griddata – below, we try out all of the interpolation methods:


grid_z0 = griddata(points, values, (grid_x, grid_y), method='nearest')
grid_z1 = griddata(points, values, (grid_x, grid_y), method='linear')
grid_z2 = griddata(points, values, (grid_x, grid_y), method='cubic')

grid_z0= np.array(grid_z0, dtype=np.float32)
grid_z1= np.array(grid_z1, dtype=np.float32)
grid_z2= np.array(grid_z2, dtype=np.float32)

import matplotlib.pyplot as plt
#plt.subplot(221)
#plt.imshow(func(grid_x, grid_y).T, extent=(0,1,0,1), origin='lower')
#plt.plot(points[:,0], points[:,1], 'k.', ms=1)
#plt.title('Original')
plt.subplot(222)
plt.imshow(grid_z0.T, extent=(0,1,0,1), origin='lower')
plt.title('Nearest')
plt.subplot(223)
plt.imshow(grid_z1.T, extent=(0,1,0,1), origin='lower')
plt.title('Linear')
plt.subplot(224)
plt.imshow(grid_z2.T, extent=(0,1,0,1), origin='lower')
plt.title('Cubic')
plt.gcf().set_size_inches(6, 6)
plt.show()





#var_range = 1

#chromo_resolution=100


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