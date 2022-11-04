from geopy.geocoders import Nominatim
import matplotlib.pyplot as plt
from Cooper import cooper
import numpy as np
from MinMaxDist import p_center_prob_w
import folium
import webbrowser


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
            loc = geolocator.geocode(city + " Midtjylland")
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
origin = (lx[0], ly[0])

x = [R * (lx_i - lx[0]) * np.pi / 180 * np.cos(np.pi * ly[0] / 180) for lx_i in lx[1:]]
y = [R * (ly_i - ly[0]) * np.pi / 180 for ly_i in ly[1:]]


def euc_to_lon_lat(point, origin):
    R = 6371
    lon = point[0]*180 / (R * np.pi * np.cos(np.pi * origin[1] / 180)) + origin[0]
    lat = point[1]*180 / (R * np.pi) + origin[1]
    return lat, lon


points = np.array(list(zip(x, y)))
weights = np.array(populations)

runs_of_each = 50

print("Cooper")
max_ds_cooper = [np.inf] * 6
best_obj_cooper = [np.inf] * 6
best_sols_cooper = [None] * 6
for _ in range(runs_of_each):
    print(_, end=" ")
    sols, objs, max_ds = [], [], []
    for p in range(5, 11):
        obj, sol, max_d = cooper(points, p, weights)
        sols.append(sol)
        objs.append(obj)
        max_ds.append(max_d)

        if obj < best_obj_cooper[p-5]:
            max_ds_cooper[p-5] = max_d
            best_obj_cooper[p-5] = obj
            best_sols_cooper[p-5] = sol
print()
print("Min max distance", max_ds_cooper)
print("Cooper obj", best_obj_cooper)

# Weighted distance
print("P-center: Weighted")
best_ds_pcw = [np.inf] * 6
best_wds_pcw = [np.inf] * 6
best_wdsols_pcw = [None] * 6
for _ in range(runs_of_each):
    print(_, end=" ")
    sols, max_ds, max_wds = [], [], []
    for p in range(5, 11):
        sol, max_d, max_wd = p_center_prob_w(points, p, weights)
        sols.append(sol)
        max_ds.append(max_d)
        max_wds.append(max_wd)

        if max_wd < best_wds_pcw[p-5]:
            best_ds_pcw[p-5] = max_d
            best_wds_pcw[p-5] = max_wd
            best_wdsols_pcw[p-5] = sol
print()
print("Min max distance", best_ds_pcw)
print("Min max weighted distance", best_wds_pcw)

# Unweighted distance
print("P-center: Unweighted")
best_ds_pc = [np.inf] * 6
best_wds_pc = [np.inf] * 6
best_dsols_pc = [None] * 6
for _ in range(runs_of_each):
    print(_, end=" ")
    sols, max_ds, max_wds = [], [], []
    for p in range(5, 11):
        sol, max_d, max_wd = p_center_prob_w(points, p, weights, unweighted=True)
        sols.append(sol)
        max_ds.append(max_d)
        max_wds.append(max_wd)

        if max_d < best_ds_pc[p-5]:
            best_ds_pc[p-5] = max_d
            best_wds_pc[p-5] = max_wd
            best_dsols_pc[p-5] = sol

print()
print("Min max distance", best_ds_pc)
print("Min max weighted distance", best_wds_pc)

# Transform solutions to lon, lat coordinates
cooper_lon_lats = []
for sol in best_sols_cooper:
    sol_translated = [euc_to_lon_lat(point, origin) for point in sol]
    cooper_lon_lats.append(np.array(sol_translated))

pcw_lon_lats = []
for sol in best_wdsols_pcw:
    sol_translated = [euc_to_lon_lat(point, origin) for point in sol]
    pcw_lon_lats.append(np.array(sol_translated))

pc_lon_lats = []
for sol in best_dsols_pc:
    sol_translated = [euc_to_lon_lat(point, origin) for point in sol]
    pc_lon_lats.append(np.array(sol_translated))


p = 5
for sol in cooper_lon_lats:
    cooper_map = folium.Map((56.302139, 9.502777), zoom_start=9)  # Create map
    for location in locations:  # Add cities to map
        folium.Marker(location).add_to(cooper_map)
    for point in sol:  # Add drone stations
        folium.Marker(list(point), icon=folium.Icon(icon="circle-dot", color="red")).add_to(cooper_map)
    cooper_map.save("../maps/cooper"+str(p)+".html")
    p += 1

p = 5
for sol in pcw_lon_lats:
    pcw_map = folium.Map((56.302139, 9.502777), zoom_start=9)  # Create map
    for location in locations:  # Add cities to map
        folium.Marker(location).add_to(pcw_map)
    for point in sol:  # Add drone stations
        folium.Marker(list(point), icon=folium.Icon(icon="circle-dot", color="red")).add_to(pcw_map)
    pcw_map.save("../maps/pcw"+str(p)+".html")
    p += 1

p = 5
for sol in pc_lon_lats:
    pc_map = folium.Map((56.302139, 9.502777), zoom_start=9)  # Create map
    for location in locations:  # Add cities to map
        folium.Marker(location).add_to(pc_map)
    for point in sol:  # Add drone stations
        folium.Marker(list(point), icon=folium.Icon(icon="circle-dot", color="red")).add_to(pc_map)
    pc_map.save("../maps/pc"+str(p)+".html")
    p += 1
