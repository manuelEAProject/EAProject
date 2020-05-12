from tkinter import *

for i in range(5): print(i)

l=list(range(5))*3
for i in l: print(i)
"""
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