import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import matplotlib.animation as animation

f = open('nbodyleap-10000500100000000 - 8 Threads.dat', 'r')
arraydata = np.loadtxt(f, delimiter=",")

N_trajectories = 500
timesimmed = 100000000
timestep = 10000
timesteps = int(timesimmed/timestep + 1)
print(timesteps)


# Set up figure & 3D axis for animation
fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1], projection='3d')
ax.axis('off')

colors = plt.cm.jet(np.linspace(0, 1, N_trajectories))

# set up lines and points

#pts = sum([ax.plot([], [], [], 'o', c=c)   for c in colors], [])

# prepare the axes limits
ax.set_xlim((-75000, 75000))
ax.set_ylim((-75000, 75000))
ax.set_zlim((-75000, 35000))
ax.set_facecolor('xkcd:black')
# set point-of-view: specified by (altitude degrees, azimuth degrees)
ax.view_init(30, 0)

scatters = []
xlist=np.zeros(N_trajectories)
ylist=np.zeros(N_trajectories)
zlist=np.zeros(N_trajectories)

extension = ".png"
name = "image"

for i in range(timesteps):
    for h in range(N_trajectories):
        np.put(xlist, h, arraydata[i, 1+6*h])
        np.put(ylist, h, arraydata[i, 3+6*h])
        np.put(zlist, h, arraydata[i, 5+6*h])
        
    scatter = ax.scatter(xlist,ylist,zlist, c="b")
    path = ("partials/%i.png" % i)
    plt.savefig(path)
    plt.cla()
    #scatters.append([scatter])

# instantiate the animator.
#anim = animation.ArtistAnimation(fig, scatters, blit=True, interval=50, repeat=False)

# Save as mp4. This requires mplayer or ffmpeg to be installed
#anim.save('lorentz_attractor.mp4', fps=15, extra_args=['-vcodec', 'libx264'])

#anim.save("attempt.mp4", dpi=300, fps=24)

print("Job Done")

f.close()
