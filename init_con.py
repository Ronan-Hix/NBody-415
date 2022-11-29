import numpy as np
import matplotlib.pyplot as plt


seed = 861621
np.random.seed(seed )
n_particles = 10000
n_sub_halos = 100
n_per_sub = int(n_particles * 0.50 / n_sub_halos)
n_central_halo = int(n_particles * 0.50)

box_size = 50  # Mpc
maximum_size_frac = 0.10  # halo can not be bigger than about 10 percent of the box
length_unit = 3e22  # 1 Mpc in meters
mass_unit = (1 / 6.674e-11) * (length_unit**3 / 3e17**2) * (1 / 2e30)  # msun
velocity_unit = 1e5  # m/s


# total mass in box assuming critical density
total_mass = 50**3 * 2.8 * 1e11 / 1000
# assume uniform dm masses, sets mass resolution
mass_per_particle = total_mass / mass_unit
masses = mass_per_particle * np.ones(n_particles)


# central halo
central_halo_size= np.random.uniform(1, box_size * maximum_size_frac)
central_halo_position = central_halo_size * np.random.randn(n_central_halo, 3)
central_halo_velocities = np.random.randn(n_central_halo, 3)

master_positions = [
    central_halo_position,
]
master_velocities = [
    central_halo_velocities,
]

# subhalos
for i in range(n_sub_halos):
    random_coord = np.random.uniform(-box_size, box_size, size=(1, 3))
    #random_size = np.random.uniform(1, box_size * maximum_size_frac)
    random_size = np.random.uniform(1, central_halo_size)
    postions = random_size * np.random.randn(n_per_sub, 3) - random_coord
    velocities = np.random.randn(n_per_sub, 3)

    master_positions.append(postions)
    master_velocities.append(velocities)

master_positions = np.concatenate(master_positions)/2
master_velocities = np.concatenate(master_velocities)


fig = plt.figure(figsize=(5, 5), dpi=300)
ax = plt.axes(projection="3d")
p = ax.scatter(
    master_positions[:, 0],
    master_positions[:, 1],
    master_positions[:, 2],
    c=np.linalg.norm(master_velocities, axis=1),
    s=0.2,
    alpha=0.8,
    cmap=plt.cm.viridis,
)

ax.set(xlim=(-50, 50), ylim=(-50, 50), zlim=(-50, 50))
fig.colorbar(p, ax=ax, orientation="horizontal", label="vmag")
plt.savefig("./init_test.png")
# plt.figure(dpi=200, figsize=(6, 6))
# plt.scatter(postions[:, 0], postions[:, 1], s=0.2)
# plt.xlim(-50,50)
# plt.ylim(-50,50)
masses = np.expand_dims(masses, axis=1)
master = np.concatenate((masses, master_positions, master_velocities), axis=1)
text_file_headr = "# mass x y z vx vy vz"
np.savetxt("init_{}_{}_{}.txt".format(n_particles,n_sub_halos, seed),master)