from geopy.geocoders import Nominatim
import matplotlib.pyplot as plt
from Cooper import cooper
import numpy as np
from MinMaxDist import p_center_prob_w

geolocator = Nominatim(user_agent="LocationPlanning", timeout=5)

with open('midt.txt', 'r', encoding='utf8') as dat:
    cities = []
    populations = []
    for line in dat:
        temp = line.split()
        if len(temp) == 2 or len(temp) == 3:
            *city, pop = temp
            city = " ".join(city)
            cities.append(city)
            if int(pop) == 0:  # Otherwise we might divide by 0
                pop = 1
            populations.append(int(pop))

locations = []

try:
    with open('locations.txt', 'r') as locdat:
        print("City locations loaded from file")
        for line in locdat:
            temp = line.split()
            locations.append((float(temp[0]), float(temp[1])))
except FileNotFoundError:
    with open('locations.txt', 'w') as locdat:
        c = 0
        print("Gathering city locations...")
        for city in cities:
            if city == "Over Hornbæk":
                city = "Kjølvejen 2, 8920"
            print(city)
            c += 1
            loc = geolocator.geocode(city + " Jylland")
            locdat.write('{} {}\n'.format(loc.latitude, loc.longitude))
            locations.append((loc.latitude, loc.longitude))
            print("{}/{} cities done".format(c, len(cities)))

# plt.plot(*zip(*locations), 'ro')
# plt.show()

lats, lons = list(zip(*locations))

lx = [min(lons)]
lx.extend(lons)
ly = [min(lats)]
ly.extend(lats)
R = 6371

x = [R * (lx_i - lx[0]) * np.pi / 180 * np.cos(np.pi * ly[0] / 180) for lx_i in lx[1:]]
y = [R * (ly_i - ly[0]) * np.pi / 180 for ly_i in ly[1:]]


points = np.array(list(zip(x, y)))
weights = np.array(populations)

sols, objs, max_ds = [], [], []
for p in range(5, 11):
    obj, sol, max_d = cooper(points, p, weights)
    sols.append(sol)
    objs.append(obj)
    max_ds.append(max_d)

# for obj, sol, max_d, p in zip(objs, sols, max_ds, range(5, 11)):
#     print("Solution with {} facilities:".format(p))
#     print("----------------------------------------")
#     print("Objective value:", obj)
#     print("Largest distance:", max_d)
    # cooper(points, p, weights, print_last=True, init_X=sol)

# Weighted distance
sols, max_ds, max_wds = [], [], []
for p in range(5, 11):
    print(p)
    sol, max_d, max_wd = p_center_prob_w(points, p, weights)
    sols.append(sol)
    max_ds.append(max_d)
    max_wds.append(max_wd)

print(max_ds)

# Unweighted distance
sols, max_ds, max_wds = [], [], []
for p in range(5, 11):
    print(p)
    no_weights = np.array([1] * len(weights))
    sol, max_d, max_wd = p_center_prob_w(points, p, no_weights)
    sols.append(sol)
    max_ds.append(max_d)
    max_wds.append(max_wd)

print(max_ds)
