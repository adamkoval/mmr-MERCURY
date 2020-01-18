import os

#g = 1
#gs = []
#for i in range(g, g+10):
#    print(g)
#    gs.append(g)
#    g+=1

runs_there = os.listdir('mcm/completed/21/info/')
runs_wanted = 6
run = runs_there + 1
while run < runs_wanted:
    print(run)
    run = runs_there + 1
    runs_there = os.listdir('mcm/completed/21/info/')
