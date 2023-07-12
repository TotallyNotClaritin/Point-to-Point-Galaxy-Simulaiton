import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

f  = open('nbodyleap-10000500100000000 - 8 Threads.dat', 'r')
arraydata = np.loadtxt(f, delimiter=",")

t = arraydata[:, 0]
x = arraydata[:, 1]
y = arraydata[:, 3]
z = arraydata[:, 5]
#energy = arraydata[:, -1]

#print(energy)

#plt.plot(t,energy)
#plt.title("Total Energy vs Time")
#plt.ylabel("Total Energy")
#plt.xlabel("Seconds")
n = 1



N_trajectories = 500



# Set up figure & 3D axis for animation
fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1], projection='3d')
ax.axis('off')

# choose a different color for each trajectory
colors = plt.cm.autumn(np.linspace(0, 1, N_trajectories))

# set up lines and points

pts = sum([ax.plot([], [], [], 'o', c=c, markersize = 1) 
           for c in colors], [])

ax.set_facecolor('xkcd:black')
# prepare the axes limits
ax.set_xlim((-100000, 100000))
ax.set_ylim((-100000, 100000))
ax.set_zlim((-100000, 100000))

# set point-of-view: specified by (altitude degrees, azimuth degrees)
ax.view_init(0, 90)

# initialization function: plot the background of each frame
def init():
    for pt in pts:
        
        pt.set_data([], [])
        pt.set_3d_properties([])
    return pts

i = 0
# animation function.  This will be called sequentially with the frame number
def animate(i):
    # we'll step two time-steps per frame.  This leads to nice results.
    i = i + 1
    h = 0
    for pt in pts:
            x1 = arraydata[:i, 1 + 6*h]
            y1 = arraydata[:i, 3 + 6*h]
            z1 = arraydata[:i, 5 + 6*h]
            pt.set_data(x1[-1:], y1[-1:])
            pt.set_3d_properties(z1[-1:])
            h = h + 1
            
    fig.canvas.draw()

    return pts

# instantiate the animator.
anim = animation.FuncAnimation(fig, animate, init_func=init, interval=1, frames=100001, blit=True, repeat = 'False', cache_frame_data=False)

anim.save("animated.mp4", fps=120)
f.close()
