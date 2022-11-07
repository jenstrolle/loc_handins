# Handins from Applied Optimization: Location Planning

This repository is dedicated to the code and results from the take-home exercises presented in the course Applied Optimization: Location Planning in fall 2022 at the Department of Math at Aarhus University.

## Handin 1
In the folder ```h1``` one can find everything pertaining to the first handin. This exercise is based on a problem from the Danish Ministry of Defense regarding placement of drone stations.
The subfolder ```code``` contains all of the code, while ```maps``` contains the resulting maps from the solutions using the python module ```folium```. Finally ```plots``` contains some plots of the resulting key values from the best solutions we found.

The file ```report.pdf``` is a written report detailing our findings.

## Code Structure
The code is split into 4 files. ```main.py```, ```Cooper.py```, ```SingleWeber.py``` and ```MinMaxDist.py```. To run the code one should run the file ```main.py```, when then calls functions from the other three files. Data handling, plots and creation of maps are also found in this file. ```SingleWeber.py``` contains multiple methods for solving the single Weber problem, specifically the Ostresh method, which is called in ```Cooper.py``` to solve the Multi-Weber problem. The $p$-center problem is solved in the file ```MinMaxDist.py```, which can handle both the weighted case as well as the unweighted case by giving weights all equal to 1.

## Map Structure
All the maps can be found in the subfolder ```maps```. Each map has all cities marked by small blue dots, while the placed drone stations are marked with a red marker. The maps are named as ```prefix + number of drone stations.html```, where the prefix can be ```pc``` for the unweighted $p$-center problem, ```pcw``` for the weighted p-center problem and finally ```cooper``` for the results of the Cooper heuristic. As an example, when trying to place 8 drone stations, the locations found by Cooper's allocation algorithm can be found in the file ```cooper8.html```.

This repository was created by Jens Trolle, Maria Askholm Nielsen and Nikolaj Kjær-Sørensen.
