# Handins from Applied Optimization: Location Planning

This repository is dedicated to the code and results from the take-home exercises presented in the course Applied Optimization: Location Planning in fall 2022 at the Department of Math at Aarhus University.

## Handin 1
In the folder ```h1``` one can find everything pertaining to the first handin. This exercise is based on a problem from the Danish Ministry of Defense regarding placement of drone stations.
The subfolder ```code``` contains all of the code, while ```maps``` contains the resulting maps from the solutions using the python module ```folium```. Finally ```plots``` contains some plots of the resulting key values from the best solutions we found.

The file ```report.pdf``` is a written report detailing our findings.

## Code Structure
Niller skriver noget klogt her

## Map Structure
All the maps can be found in the subfolder "maps". Each map has all cities marked by small blue dots, while the placed drone stations are marked with a red marker. The maps are named as "prefix" + "number of drone stations", where the prefix can be "pc" for the unweighted $$p$$-center problem, "pcw" for the weighted p-center problem and finally "cooper" for the results of the Cooper heuristic. As an example, when trying to place 8 drone stations, the locations found by Cooper's allocation algorithm can be found in the file "cooper8.html".

This repository was created by Jens Trolle, Maria Askholm Nielsen and Nikolaj Kjær-Sørensen.
