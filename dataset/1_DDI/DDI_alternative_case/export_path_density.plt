set title "Dynamic Density Contour (Path 1 ) Unit: veh/mile/lane" 
set xlabel "Time Horizon"
set ylabel "Space (Node Sequence)"  offset -1
set xtics (" 6:00" 0 ," 6:10" 10 ," 6:20" 20 ," 6:30" 30 ," 6:40" 40 ," 6:50" 50 ," 7:00" 60 ," 7:10" 70 ," 7:20" 80 ," 7:30" 90 ," 7:40" 100 ," 7:50" 110 ," 8:00" 120 ) 
set ytics ("2560" 0, "2571" 801, "329" 1537, "789" 2239, "2299" 2910, "2353" 3518)
set xrange [0:121] 
set yrange [0:3518] 
set palette defined (0 "white", 10 "green", 30 "yellow", 50 "red")
set pm3d map
splot 'C:\Users\18016\Documents\GitHub\DLSim\dataset\1_DDI\DDI_alternative_case\export_path_density.txt' matrix notitle
