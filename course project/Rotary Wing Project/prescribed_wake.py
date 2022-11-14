##############################################################################
# Imports
##############################################################################
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import matplotlib.animation as animation

##############################################################################
# Code
##############################################################################
pi = np.pi

# Parameter initialization
R = 0.8255
C_T = 0.0050351
SIGMA = 0.0979
b = 4  # number of blades
theta_twist = -13
omega = 1 / 10  # in rads/sec


def get_vortex_loc(
    phi_v, theta_twist=theta_twist, C_T=C_T, SIGMA=SIGMA, b=b, R=R
):
    TWIST = theta_twist
    k1 = -0.25 * (C_T / SIGMA + 0.001 * TWIST)
    k2 = -((C_T) ** 0.5) - 0.01 * ((C_T) ** 0.5) * TWIST
    Zv = 0
    scaling = 1

    if phi_v < 2 * pi / b:
        Zv = R * (phi_v) * k1 * scaling
    elif phi_v > 2 * pi / b:
        Zv = R * scaling * (k1 * 2 * pi / b + k2 * (phi_v - 2 * pi / b))

    A = 0.78
    lamda = 0.145 + 27 * C_T
    Rv = (A + (1 - A) * np.exp(-lamda * phi_v)) * R

    # Xv = Rv * np.cos(phi_v)
    # Yv = Rv * np.sin(phi_v)

    return (Rv, Zv)


class Blade:
    def __init__(self, blade_id, num_of_blades, radius=R, omega=omega):
        self.blade_id = blade_id
        self.num_of_blades = num_of_blades
        self.root_loc = np.zeros(3)
        self.tip_theta = (blade_id / num_of_blades) * 2 * np.pi
        self.tip_loc = np.array(
            [
                self.root_loc[0] + radius * np.cos(self.tip_theta),
                self.root_loc[1] + radius * np.sin(self.tip_theta),
                self.root_loc[2],
            ]
        )
        self.omega = omega
        self.radius = radius

    def update_tip_loc(self, time):
        self.tip_theta = (
            self.blade_id / num_of_blades
        ) * 2 * np.pi - self.omega * time
        self.tip_loc = np.array(
            [
                self.root_loc[0] + self.radius * np.cos(self.tip_theta),
                self.root_loc[1] + self.radius * np.sin(self.tip_theta),
                self.root_loc[2],
            ]
        )

    def get_tip_theta(self):
        return self.tip_theta

    def get_tip_loc(self):
        return self.tip_loc

    def get_root_loc(self):
        return self.root_loc


class Vortex:
    def __init__(self, vortex_id, origin_blade, phi_v, phi_v_max):
        self.origin_blade = origin_blade
        self.vortex_id = vortex_id
        self.phi_v_0 = phi_v
        self.phi_v_curr = phi_v
        self.theta = phi_v + origin_blade.get_tip_theta()
        self.R, self.Z = get_vortex_loc(self.phi_v_curr)
        self.phi_v_max = phi_v_max

    def update_vortex_loc(self, time):
        self.phi_v_curr = self.phi_v_0 + omega * time
        if self.phi_v_curr > self.phi_v_max:
            self.phi_v_0 = -omega * time
            self.phi_v_curr = self.phi_v_0 + omega * time
            self.theta = self.origin_blade.get_tip_theta()
        self.R, self.Z = get_vortex_loc(self.phi_v_curr)

    def get_vortex_loc(self):
        vortex_loc = np.array(
            [self.R * np.cos(self.theta), self.R * np.sin(self.theta), self.Z]
        )
        return vortex_loc


def update_structure(time, blade, vortex_arr):
    blade.update_tip_loc(time)
    for vortex in vortex_arr:
        vortex.update_vortex_loc(time)


num_of_vortices = 60
phi_v_max = 3 * pi
phi_v_arr = np.linspace(0, phi_v_max, num_of_vortices)
num_time_steps = 120
time_arr = np.linspace(0, 120, num_time_steps)

num_of_blades = b

blade_arr = [Blade(i, num_of_blades) for i in range(num_of_blades)]
vortex_arr = [
    [
        Vortex(i, blade_arr[j], phi_v_arr[i], phi_v_max)
        for i in range(num_of_vortices)
    ]
    for j in range(num_of_blades)
]

##############################################################################
# Time evolution
##############################################################################

# Attaching 3D axis to the figure
fig = plt.figure()
ax = axes3d.Axes3D(fig)


# NOTE: Can't pass empty arrays into 3d version of plot(); for marker="o"
points = [
    [
        ax.plot(
            vortex_arr[blade_id][vortex_id].get_vortex_loc()[0],
            vortex_arr[blade_id][vortex_id].get_vortex_loc()[1],
            vortex_arr[blade_id][vortex_id].get_vortex_loc()[2],
            marker="o",
            c="r",
        )[0]
        for vortex_id in range(num_of_vortices)
    ]
    for blade_id in range(num_of_blades)
]

blade_root = ax.plot(
    [blade_arr[0].get_root_loc()[0],],
    [blade_arr[0].get_root_loc()[1],],
    [blade_arr[0].get_root_loc()[2],],
    "og",
)[0]

rotor_axis = ax.plot(
    [blade_arr[0].get_root_loc()[0], blade_arr[0].get_root_loc()[0]],
    [blade_arr[0].get_root_loc()[1], blade_arr[0].get_root_loc()[1]],
    [blade_arr[0].get_root_loc()[2], blade_arr[0].get_root_loc()[2] - 0.4],
    "--g",
)[0]

# blade_lines = [
#     ax.plot(
#         [
#             blade_arr[blade_id].get_root_loc()[0],
#             blade_arr[blade_id].get_tip_loc()[0],
#         ],
#         [
#             blade_arr[blade_id].get_root_loc()[1],
#             blade_arr[blade_id].get_tip_loc()[1],
#         ],
#         [
#             blade_arr[blade_id].get_root_loc()[2],
#             blade_arr[blade_id].get_tip_loc()[2],
#         ],
#         "b",
#     )[0]
#     for blade_id in range(num_of_blades)
# ]

# draw = []
# for i in range(len(points)):
#     for j in range(len(points[i])):
#         draw.append(points[i][j])
# for i in range(len(blade_lines)):
#     draw.append(blade_lines[i])


def update_plot(iter):
    for blade_id in range(num_of_blades):
        update_structure(
            time_arr[iter], blade_arr[blade_id], vortex_arr[blade_id]
        )

        for vortex_id in range(num_of_vortices):
            points[blade_id][vortex_id].set_data(
                vortex_arr[blade_id][vortex_id].get_vortex_loc()[0],
                vortex_arr[blade_id][vortex_id].get_vortex_loc()[1],
            )
            points[blade_id][vortex_id].set_3d_properties(
                vortex_arr[blade_id][vortex_id].get_vortex_loc()[2]
            )

        # blade_lines[blade_id].set_data(
        #     [
        #         blade_arr[blade_id].get_root_loc()[0],
        #         blade_arr[blade_id].get_tip_loc()[0],
        #     ],
        #     [
        #         blade_arr[blade_id].get_root_loc()[1],
        #         blade_arr[blade_id].get_tip_loc()[1],
        #     ],
        # )

        # blade_lines[blade_id].set_3d_properties(
        #     [
        #         blade_arr[blade_id].get_root_loc()[2],
        #         blade_arr[blade_id].get_tip_loc()[2],
        #     ]
        # )

        if iter > 0:
            blade_lines = [
                ax.plot(
                    [
                        blade_arr[blade_id].get_root_loc()[0],
                        blade_arr[blade_id].get_tip_loc()[0],
                    ],
                    [
                        blade_arr[blade_id].get_root_loc()[1],
                        blade_arr[blade_id].get_tip_loc()[1],
                    ],
                    [
                        blade_arr[blade_id].get_root_loc()[2],
                        blade_arr[blade_id].get_tip_loc()[2],
                    ],
                    "b",
                )[0]
                for blade_id in range(num_of_blades)
            ]
            # blade_lines[blade_id].set_data(
            #     [
            #         blade_arr[blade_id].get_root_loc()[0],
            #         blade_arr[blade_id].get_tip_loc()[0],
            #     ],
            #     [
            #         blade_arr[blade_id].get_root_loc()[1],
            #         blade_arr[blade_id].get_tip_loc()[1],
            #     ],
            # )

    draw = []
    for i in range(len(points)):
        for j in range(len(points[i])):
            draw.append(points[i][j])
    if iter > 0:
        for i in range(num_of_blades):
            draw.append(blade_lines[i])

    return draw


# Setting the axes properties
ax.set_xlim3d([-1.0, 1.0])
ax.set_xlabel("X")

ax.set_ylim3d([-1.0, 1.0])
ax.set_ylabel("Y")

ax.set_zlim3d([-0.35, 0.1])
ax.set_zlabel("Z")

ax.set_title("Prescribed Wake Theory")


def on_press(event):
    if event.key.isspace():
        if vortex_ani.running:
            vortex_ani.event_source.stop()
        else:
            vortex_ani.event_source.start()
        vortex_ani.running ^= True


fig.canvas.mpl_connect("key_press_event", on_press)

# Creating the Animation object
vortex_ani = animation.FuncAnimation(
    fig, update_plot, num_time_steps, interval=50, blit=True, repeat=False,
)

SAVE = False
if SAVE == True:
    writer = animation.FFMpegWriter(fps=20, metadata=dict(artist="Me"))
    vortex_ani.save(
        "D:/nakul-pavilion/Google Drive/IITB/SS/Sem 6/AE667 (Rotary Wing Aerodynamics)/Project/vortex_evolution.mov",
        writer=writer,
    )

vortex_ani.running = False

plt.show()

