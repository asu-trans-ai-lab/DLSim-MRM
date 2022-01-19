set title "Space time trajectory diagram"
set xlabel "Time Horizon"
set ylabel "Space (distance)"  offset -1
set xtics (" 6:00" 0 ," 6:10" 10 ," 6:20" 20 ," 6:30" 30 ," 6:40" 40 ," 6:50" 50 ," 7:00" 60 ," 7:10" 70 ," 7:20" 80 ," 7:30" 90 ," 7:40" 100 ," 7:50" 110 ," 8:00" 120 ) 
set ytics (" " -2147483648)
set xrange [0:121] 
set yrange [0:351.14] 
plot "agent1.txt" using 1:2 title 'agent 1'  with lines,\
"agent2.txt" using 1:2 title 'agent 2'  with lines,\
"agent3.txt" using 1:2 title 'agent 3'  with lines,\
"agent4.txt" using 1:2 title 'agent 4'  with lines
