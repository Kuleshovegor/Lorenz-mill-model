import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.constants as sc


def lorenz(x, y, z, s, b, r):
    x_diff = s * (y - x)
    y_diff = r * x - y - x * z
    z_diff = x * y - b * z
    return x_diff, y_diff, z_diff


def get_lorenz_attractor_data(dt, num_steps, start_x, start_y, start_z):
    x_array = np.empty(num_steps + 1)
    y_array = np.empty(num_steps + 1)
    z_array = np.empty(num_steps + 1)

    x_array[0], y_array[0], z_array[0] = (start_x, start_y, start_z)

    for i in range(num_steps):
        x_diff, y_diff, z_diff = lorenz(x_array[i], y_array[i], z_array[i], sigma, b, r)
        x_array[i + 1] = x_array[i] + (x_diff * dt)
        y_array[i + 1] = y_array[i] + (y_diff * dt)
        z_array[i + 1] = z_array[i] + (z_diff * dt)
    return x_array, y_array, z_array


def plot_lorenz_attractor(x_array, y_array, z_array):
    figure1 = plt.figure()
    attractor = figure1.gca(projection='3d')
    attractor.plot(x_array, y_array, z_array, lw=0.05)
    attractor.set_xlabel("X")
    attractor.set_ylabel("Y")
    attractor.set_zlabel("Z")
    attractor.set_title("Lorenz Attractor")
    plt.show()


def lorenz_mill(a, b, w, k, alpha, q, v, r, inerc):
    a_diff = w * b - k * a
    b_diff = -w * a - k * b + q
    w_diff = (-v / inerc) * w + ((math.pi * sc.g * math.sin(alpha) * r) / inerc) * a
    return a_diff, b_diff, w_diff


def get_lorenz_mill_data(dt, num_steps, k, alpha, q, v, r, inerc, start_a, start_b, start_w):
    a_array = np.empty(num_steps + 1)
    b_array = np.empty(num_steps + 1)
    w_array = np.empty(num_steps + 1)
    t_array = np.empty(num_steps + 1)

    a_array[0], b_array[0], w_array[0] = (start_a, start_b, start_w)
    t_array[0] = 0

    for i in range(num_steps):
        a_diff, b_diff, w_diff = lorenz_mill(a_array[i], b_array[i], w_array[i], k=k, alpha=alpha, q=q, v=v, r=r,
                                             inerc=inerc)
        a_array[i + 1] = a_array[i] + (a_diff * dt)
        b_array[i + 1] = b_array[i] + (b_diff * dt)
        w_array[i + 1] = w_array[i] + (w_diff * dt)
        t_array[i + 1] = t_array[i] + dt
    return a_array, b_array, w_array, t_array


def plot_lorenz_mill_graphs(a_array, b_array, w_array, t_array):
    fig1 = plt.figure()
    fig2 = plt.figure()
    fig3 = plt.figure()
    attractor = fig1.gca(projection='3d')
    angular_velocity = fig2.gca()
    velocity_distribution = fig3.gca()
    max_velocity = max(int(np.max(w_array) + 1), abs(int(np.min(w_array) + 1)))
    values = np.arange(-max_velocity, max_velocity, 0.5)
    num_values = np.empty(4 * max_velocity)
    for i in range(4 * max_velocity):
        num_values[i] = 0
    for i in range(len(w_array)):
        mn = max_velocity
        mn_value = 0
        for c in values:
            if abs(w_array[i] - c) < mn:
                mn = abs(w_array[i] - c)
                mn_value = c
        num_values[int(2 * (mn_value + max_velocity)) - 1] += 1

    velocity_distribution.plot(values, num_values)

    angular_velocity.plot(t_array, w_array, lw=0.1)
    attractor.plot(a_array, b_array, w_array, lw=0.1)
    attractor.set_xlabel("A")
    attractor.set_ylabel("B")
    attractor.set_zlabel("W")
    attractor.set_title("Lorenz Mill Attractor")
    velocity_distribution.set_xlabel("Angular velocity")
    velocity_distribution.set_ylabel("Number of repetitions")
    velocity_distribution.set_title("Velocity distribution")
    angular_velocity.set_xlabel("Time")
    angular_velocity.set_ylabel("Angular velocity")
    angular_velocity.set_title("Angular velocity")

    plt.show()


start_x = 1.
start_y = 0.
start_z = 0.

sigma = 10
r = 28
b = 8. / 3.

dt = 0.01
num_steps = 100000

x_array, y_array, z_array = get_lorenz_attractor_data(dt, num_steps, start_x, start_y, start_z)
plot_lorenz_attractor(x_array, y_array, z_array)

start_a = 0.
start_b = 0.
start_w = 0.1

k = 1
alpha = math.pi / 7
q = 2.5
v = 3.3
r = 1.4
inerc = 1

dt = 0.01
num_steps = 100000

x_array, y_array, z_array, t_array = get_lorenz_mill_data(dt, num_steps, k, alpha, q, v, r, inerc, start_a, start_b, start_w)
plot_lorenz_mill_graphs(x_array, y_array, z_array, t_array)

x_array, y_array, z_array, t_array = get_lorenz_mill_data(dt, num_steps, k, alpha, q, v, r, inerc, start_a, start_b, start_w + 0.5)
plot_lorenz_mill_graphs(x_array, y_array, z_array, t_array)
