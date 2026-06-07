from array import array
import matplotlib.pyplot as plt
import csv
import numpy as np

X = []
Uinit = []
Ufinal = []

# opening the CSV file
with open('initialCondition.csv', mode ='r')as file:
   
  # reading the CSV file
  csvFile = csv.reader(file)
 
  # displaying the contents of the CSV file
  for lines in csvFile:
        X.append(float(lines[0]))
        Uinit.append(float(lines[1]))


with open('final.csv', mode ='r')as file:
   
  # reading the CSV file
  csvFile = csv.reader(file)
 
  # displaying the contents of the CSV file
  for lines in csvFile:
        Ufinal.append(float(lines[1]))
        
#print(X, Uinit)
plt.plot(X, Uinit, '.-', color = 'k', label="Init")
plt.plot(X, Ufinal, '.-', color = 'r', label="Final")
plt.legend()
plt.grid()
plt.show()