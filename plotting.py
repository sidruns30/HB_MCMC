import matplotlib.pyplot as plt

filename = "110602878.out"
data = []

with open(filename, "r") as f:
    for i,line in enumerate(f.readlines()):
        if (i==0): pass
        else: data.append(line.split(' '))

time, lc, model = [], [], [] 
for d in data:
    time.append(float(d[0]))
    lc.append(float(d[1]))
    model.append(float(d[2]))

plt.style.use("dark_background")
plt.figure(figsize=(10, 5))
plt.title("Computed vs model lightcurve", fontsize=30, font="serif")
plt.plot(time, lc, 'ro')
plt.plot(time, model, 'w--', linewidth=3)
plt.xlabel("Time", fontsize=15, font="serif")
plt.ylabel("Flux", fontsize=15, font="serif")
plt.xlim((1480, 1520))
plt.ylim((0.8, 1.2))
plt.show()