import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors

steps = 10000
count = 0
statics = []
with open("{:}_data_export.dat".format(steps), 'r') as f:
    while True:
        line = f.readline()
        # if line is empty
        # end of file is reached
        if not line:
            break
        if 'static_{:}'.format(steps) not in line.strip().split():
            statics.append(count)
        count +=1
print(statics)

with open("{:}_data_export.dat".format(steps), 'r') as f:
    array = np.genfromtxt(f)

viridis = cm.get_cmap('viridis')
cmap = plt.get_cmap("viridis")


#567                      -1732                    0                        harmonic_100000          10                       700                      0.0203118399

ring_limits = array[:,0]
t_list = array[:,1]
i_list = array[:,4]
a_list = array[:,5]
e_list = array[:,6]



dict = {}
for i in statics:

    T = t_list[int(i)]
    if str(T) not in dict.keys():
        dict['{:}'.format(T)] = {}

    vals = i_list[int(i)]

    if str(vals) not in dict['{:}'.format(T)].keys() and int(float(ring_limits[int(i)])) == 567:
        dict['{:}'.format(T)]['{:}'.format(vals)] = {}
        dict['{:}'.format(T)]['{:}'.format(vals)]['Area'] = []
        dict['{:}'.format(T)]['{:}'.format(vals)]['Energy'] = []

    if int(float(ring_limits[int(i)])) == 567:
        dict['{:}'.format(T)]['{:}'.format(vals)]['Area'].append(a_list[i])
        dict['{:}'.format(T)]['{:}'.format(vals)]['Energy'].append(np.log10(e_list[i]))

print(dict)

norm_angle = matplotlib.colors.Normalize(vmin=0, vmax=100)

print(norm_angle(5.0))
for T in dict.keys():

    intercepts = sorted(dict['{:}'.format(T)].keys())
    list_intercepts = [float(intercept) for intercept in intercepts]
    list_intercepts = sorted(list_intercepts)
    print(intercepts)
    print(list_intercepts)
    for intercept in list_intercepts:
#        if int(float(intercept)) == 45:
            x = dict['{:}'.format(T)]['{:}'.format(intercept)]['Area']
            y = dict['{:}'.format(T)]['{:}'.format(intercept)]['Energy']
            y = [i + 0*float(intercept) for i in y]
            x, y = zip(*sorted(zip(x, y)))

            plt.plot(x,y,color=viridis(norm_angle(float(intercept))), label='{:}'.format(int(float(intercept))))
    plt.legend()
    plt.title('Temperature : {:}  at {:} steps'.format(T, steps))
    #plt.ylim(0,250)
    plt.show()

def calc_E_harm():

    return