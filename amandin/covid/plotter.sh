# /bin/zsh

gcc -g model.c -lm
./a.out
gnuplot -e "plot 'new_cases.txt' with lines title 'Simulated Data', 'data_7_day_avg.csv' with lines title 'Real Data 7 Day Rolling Average'" -p
