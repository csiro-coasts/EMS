import matplotlib.pyplot as plt
import csv
x1=[]
y1=[]
with open('mass.txt','r') as csvfile:
    next(csvfile)
    plots = csv.reader(csvfile,delimiter=' ')
    for row in plots:
        x1.append(float(row[0]))
        y1.append(float(row[3]))
x2=[]
y2=[]
with open('mass.txt','r') as csvfile:
    next(csvfile)
    plots = csv.reader(csvfile,delimiter=' ')
    for row in plots:
        x2.append(float(row[0]))
        y2.append(float(row[4]))

x3=[]
y3=[]
with open('mass.txt','r') as f:
    next(f)
    plots = csv.reader(f,delimiter=' ')
    for row in plots:
        x3.append(float(row[0]))
        y3.append(float(row[7]))




#plt.plot(x1,y1,label="3")
#plt.plot(x2,y2,label="4")
plt.plot(x3,y3)
plt.xlabel('Seconds')
plt.ylabel('Mass (in some units)')
plt.title('Mud')
plt.legend()
plt.show()
