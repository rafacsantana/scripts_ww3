import numpy as np
import matplotlib.pyplot as plt

plt.ion()

#dx/dt = cg * cos(theta)
#dy/dt = cg * sin(theta)
#dtheta/dt = -dk/dy * cg

def ray_trace_step(x, y, theta, dt, bathymetry_function):
    # Calculate group velocity and wavenumber based on bathymetry
    # Update x, y, and theta using the ray equations and dt
    return new_x, new_y, new_theta

# Example (simplified, no bathymetry or currents)
def deep_water_group_velocity(wave_period):
    g = 9.81  # Acceleration due to gravity
    return g * wave_period / (4 * np.pi)

def ray_trace_step(x, y, theta, dt, wave_period):
    cg = deep_water_group_velocity(wave_period)
    new_x = x + cg * np.cos(theta) * dt
    new_y = y + cg * np.sin(theta) * dt
    new_theta = theta # No change in direction for this simplified example
    return new_x, new_y, new_theta

# Parameters
wave_period = 10 #seconds
initial_thetas = np.linspace(0, 2 * np.pi, 10) # Radians
dt = 10 #seconds
num_steps = 100

# Initial positions
x0 = np.zeros_like(initial_thetas)
y0 = np.zeros_like(initial_thetas)

# Store ray paths
x_paths = np.zeros((len(initial_thetas), num_steps))
y_paths = np.zeros((len(initial_thetas), num_steps))

# Ray tracing loop
for i, initial_theta in enumerate(initial_thetas):
    x = x0[i]
    y = y0[i]
    theta = initial_theta
    for j in range(num_steps):
      x, y, theta = ray_trace_step(x, y, theta, dt, wave_period)
      x_paths[i,j] = x
      y_paths[i,j] = y

# Plotting
plt.figure()
for i in range(len(initial_thetas)):
    plt.plot(x_paths[i,:], y_paths[i,:])
plt.xlabel("x")
plt.ylabel("y")
plt.title("Ocean Wave Ray Tracing")
plt.grid(True)
plt.show()

