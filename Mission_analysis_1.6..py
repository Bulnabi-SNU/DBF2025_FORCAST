import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

#### Test ####
### Constants ###
rho = 1.2  # air density (kg/m^3)
g = 9.81  # gravity (m/s^2)
m_glider = 7  # glider mass (kg)
m_payload = 2.5  # payload mass (kg)
m_x1 = 0.2  # additional mass (kg)
W = (m_glider + m_payload + m_x1) * g  # total weight (N)
m = m_glider + m_payload + m_x1  # total mass (kg)
S = 0.6  # wing area (m^2)
AR = 7.2  # aspect ratio
lw = 0.2
lh = 1
CD0 = 0.02  # zero-lift drag coefficient
CL_0 = 0.0 * (lh - lw) / lh  # lift coefficient at zero angle of attack
CL_alpha = 0.086 * (lh - lw) / lh  # lift coefficient gradient per degree
e = 0.8  # Oswald efficiency factor
T_max = 6.6 * g  # maximum thrust (N)
V_stall = 15.7  # stall speed (m/s) 15.7

# # Climb constants
# theta_deg = 40  # pitch angle
# theta_rad = math.radians(theta_deg)

# Turn constants
phi_deg = 20  # bank angle
phi_rad = math.radians(phi_deg)
total_angle_rad = 0

### Helper Functions ###
def magnitude(vector):
    return math.sqrt(sum(x*x for x in vector))

def calculate_induced_drag(C_L):
    return (C_L**2) / (math.pi * AR * e)

def calculate_cruise_alpha_w(v):
    speed = magnitude(v)
    def equation(alpha_w):
        CL = CL_0 + CL_alpha * alpha_w
        L = 0.5 * rho * speed**2 * S * CL
        return L + T_max * 0.9 * math.sin(math.radians(alpha_w)) - W
    alpha_w_solution = fsolve(equation, 5, xtol=1e-8, maxfev=1000)
    return alpha_w_solution[0]

def calculate_turn_speed(desired_radius):
    """
    Calculate the appropriate speed for a given desired turn radius
    using the bank angle (phi_rad)
    """
    # From v^2 = g * R * tan(phi)
    return math.sqrt(g * desired_radius * math.tan(phi_rad))
    speed = magnitude(v)
    q = 0.5 * rho * speed**2
    K = 1 / (math.pi * AR * e)
    load_factor = math.sqrt(q / (K * (W / S)) * ((T_max * 0.55 / W) - (q * CD0 / (W / S))))
    return min(load_factor, 1.5)

def calculate_drag_coefficient(CL):
    return CD0 + calculate_induced_drag(CL)

### Result Lists ###
time_list = []
distance_list = []
load_factor_list = []
alpha_w_list = []
gamma_list = []
radius_list_turn = []
position_list = []
v_list = []
a_list = []
phase_index = []

### Acceleration Functions ###
def calculate_acceleration_ground(v):
    speed = magnitude(v)
    CL = 0.99 # flap on
    CDi = calculate_induced_drag(CL)
    CD = CD0 + CDi
    D = 0.5 * rho * speed**2 * S * CD
    L = 0.5 * rho * speed**2 * S * CL
    a_x = -g/W * (T_max*0.9 - D - 0.03*(W-L))
    return np.array([a_x, 0, 0])

def calculate_acceleration_cruise(v, alpha_w_deg):
    speed = magnitude(v)
    CL = CL_0 + CL_alpha * alpha_w_deg
    CDi = calculate_induced_drag(CL)
    CD = CD0 + CDi
    D = 0.5 * rho * speed**2 * S * CD
    a_x = (T_max * 0.9) / m * math.cos(math.radians(alpha_w_deg)) - D / m
    return np.array([a_x, 0, 0])

def calculate_acceleration_climb(v, alpha_w_deg, gamma_rad, theta_deg, z_pos):
    speed = magnitude(v)
    if (z_pos > 10) :
        CL = CL_0 + CL_alpha * alpha_w_deg
    else:
        CL = 0.99 + 0.06 * alpha_w_deg
    CDi = calculate_induced_drag(CL)
    CD = CD0 + CDi
    D = 0.5 * rho * speed**2 * S * CD
    L = 0.5 * rho * speed**2 * S * CL
    a_x = -(T_max * 0.9 * math.cos(math.radians(theta_deg)) - L * math.sin(math.radians(alpha_w_deg)) - D * math.cos(gamma_rad)) / m
    a_z = (T_max * 0.9 * math.sin(math.radians(theta_deg)) + L * math.cos(math.radians(alpha_w_deg)) - D * math.sin(gamma_rad) - W) / m
    return np.array([a_x, 0, a_z])

def calculate_acceleration_turn(v, turn_direction="right"):
    """Calculate acceleration during coordinated turn"""
    speed = magnitude(v)
    
    # Required lift coefficient for coordinated turn
    load_factor = 1 / math.cos(phi_rad)
    CL = (load_factor * W) / (0.5 * rho * speed**2 * S)
    
    # Calculate forces
    L = 0.5 * rho * speed**2 * S * CL
    CD = calculate_drag_coefficient(CL)
    D = 0.5 * rho * speed**2 * S * CD
    T = min(T_max * 0.55, D + W * math.sin(phi_rad))
    
    # Calculate centripetal and tangential accelerations
    a_centripetal = g * math.sqrt(load_factor**2 - 1)
    a_tangential = (T - D) / m
    
    return a_centripetal, a_tangential

### Main Simulation Functions ###
def takeoff_simulation():
    print("\nRunning Takeoff Simulation...")
    dt = 0.01
    v = np.array([0.0, 0.0, 0.0])
    position = np.array([0.0, 0.0, 0.0])
    t = 0.0
    
    # Ground roll until rotation speed
    while magnitude(v) < 1.3 * V_stall:
        t += dt
        time_list.append(t)
        
        a = calculate_acceleration_ground(v)
        v += a * dt
        position += v * dt
        
        L = 0.5 * rho * magnitude(v)**2 * S * CL_0
        load_factor_list.append(L / W)
        v_list.append(v.copy())
        alpha_w_list.append(0)
        a_list.append(a)
        position_list.append(tuple(position))

def climb_simulation(h_max):
    print("\nRunning Climb Simulation...")
    
    dt = 0.01
    n_steps = int(60 / dt)  # Max 60 seconds simulation
    v = v_list[-1].copy()
    d = distance_list[-1] if distance_list else 0
    t = time_list[-1]
    x_pos, y_pos, z_pos = position_list[-1]
    theta_deg = 0

    for step in range(n_steps):
        t += dt
        time_list.append(t)

        theta_deg += 0.3
        theta_deg = min(theta_deg, 15)

        # Calculate climb angle
        gamma_rad = math.atan2(abs(v[2]), abs(v[0]))
        if (step%10 == 0): print(math.degrees(gamma_rad))
        alpha_w_deg = theta_deg - math.degrees(gamma_rad)

        # Calculate load factor
        if (z_pos > 10) :
            CL = CL_0 + CL_alpha * alpha_w_deg
        else:
            CL = 0.99 + 0.06 * alpha_w_deg
        L = 0.5 * rho * magnitude(v)**2 * S * CL
        load_factor = L / W
        load_factor_list.append(load_factor)

        # RK4 integration
        a1 = calculate_acceleration_climb(v, alpha_w_deg, gamma_rad, theta_deg, z_pos)
        v1 = v + (a1*dt/2)
        a2 = calculate_acceleration_climb(v1, alpha_w_deg, gamma_rad, theta_deg, z_pos)
        v2 = v + (a2*dt/2)
        a3 = calculate_acceleration_climb(v2, alpha_w_deg, gamma_rad, theta_deg, z_pos)
        v3 = v + a3*dt
        a4 = calculate_acceleration_climb(v3, alpha_w_deg, gamma_rad, theta_deg, z_pos)
        
        a = (a1 + 2*a2 + 2*a3 + a4)/6
        v += a*dt

        # Speed limiting
        speed = magnitude(v)
        if speed > 50:
            v = v * (50 / speed)

        # Update position
        x_pos += v[0] * dt
        z_pos += v[2] * dt
        d += magnitude(v)*dt
        position_list.append((x_pos, y_pos, z_pos))

        # Store results
        v_list.append(v.copy())
        gamma_list.append(math.degrees(gamma_rad))
        alpha_w_list.append(alpha_w_deg)
        a_list.append(a)
        distance_list.append(d)

        if z_pos > h_max:
            break

def cruise_simulation(x_final, direction='+'):
    print("\nRunning Cruise Simulation...")
    dt = 0.1
    max_steps = 1000
    step = 0
    
    # Initialize vectors
    v = v_list[-1].copy()
    v[2] = 0  # Zero vertical velocity
    speed = magnitude(v)
    if direction == '+':
        v = np.array([speed, 0, 0])  # Align with x-axis
    else:
        v = np.array([-speed, 0, 0])
        
    t = time_list[-1]
    x_pos, y_pos, z_pos = position_list[-1]
    d = distance_list[-1] if distance_list else 0
    
    while step < max_steps:
        step += 1
        t += dt
        time_list.append(t)
        
        # Calculate alpha_w first
        alpha_w_deg = calculate_cruise_alpha_w(v)
        
        # RK4 integration
        a1 = calculate_acceleration_cruise(v, alpha_w_deg)
        v1 = v + a1 * dt / 2
        a2 = calculate_acceleration_cruise(v1, alpha_w_deg)
        v2 = v + a2 * dt / 2
        a3 = calculate_acceleration_cruise(v2, alpha_w_deg)
        v3 = v + a3 * dt
        a4 = calculate_acceleration_cruise(v3, alpha_w_deg)
        
        a = (a1 + 2 * a2 + 2 * a3 + a4) / 6
        if direction == '+': v += a * dt
        else: v -= a * dt
            
        # Speed limiting while maintaining direction
        speed = magnitude(v)
        if speed > 50:  # Original speed limit
            v = v * (50 / speed)
            
        # Update position
        dx = v[0] * dt
        dy = v[1] * dt
        x_pos += dx
        y_pos += dy
        d += math.sqrt(dx*dx + dy*dy)
        position_list.append((x_pos, y_pos, z_pos))
        
        # Calculate and store results
        CL = CL_0 + CL_alpha * alpha_w_deg
        L = 0.5 * rho * speed**2 * S * CL
        
        # Store results
        load_factor_list.append(L / W)
        v_list.append(v.copy())
        alpha_w_list.append(alpha_w_deg)
        a_list.append(a)
        distance_list.append(d)
        
        # Check if we've reached target x position
        if direction == '+':
            if x_pos >= x_final:
                break
        else:
            if x_pos <= x_final:
                break

def turn_simulation(target_heading, desired_radius, direction="right"):
    print("\nRunning Turn Simulation...")
    
    dt = 0.01
    v = v_list[-1].copy()
    v[2] = 0  # Zero vertical velocity
    
    # Define desired turn radius and calculate required speed
    turn_speed = calculate_turn_speed(desired_radius)
    print(f"Turn speed for {desired_radius}m radius: {turn_speed:.2f} m/s")

    # Gradually reduce speed to turn_speed
    speed = magnitude(v)
    
    t = time_list[-1]
    x_pos, y_pos, z_pos = position_list[-1]
    d = distance_list[-1]
    
    # Get initial heading and setup turn center
    initial_heading = math.atan2(v[1], v[0])
    current_heading = initial_heading

    load_factor = 1 / math.cos(phi_rad)  # For coordinated turn
    
    # Calculate actual turn radius based on final speed
    R = turn_speed**2 / (g * math.tan(phi_rad))
    omega = turn_speed / R  # Angular velocity

    print(f"Actual turn radius: {R:.2f} m")

    # Calculate turn center
    if direction == "right":
        center_x = x_pos - R * math.sin(initial_heading)
        center_y = y_pos + R * math.cos(initial_heading)
    else:
        center_x = x_pos + R * math.sin(initial_heading)
        center_y = y_pos - R * math.cos(initial_heading)
    
    # Initialize turn tracking
    target_angle = math.radians(target_heading)
    turned_angle = 0

    # Execute turn at constant speed
    while abs(turned_angle) < abs(target_angle):
        t += dt
        time_list.append(t)
        
        # Calculate accelerations
        a_centripetal, a_tangential = calculate_acceleration_turn(v, direction)

        # Update heading based on angular velocity
        if direction == "right":
            current_heading += omega * dt
            turned_angle += omega * dt
        else:
            current_heading -= omega * dt
            turned_angle -= omega * dt
        
        # Calculate new position relative to turn center
        if direction == "right":
            x_pos = center_x + R * math.sin(current_heading)
            y_pos = center_y - R * math.cos(current_heading)
        else:
            x_pos = center_x - R * math.sin(current_heading)
            y_pos = center_y + R * math.cos(current_heading)

        # Update velocity direction (tangent to the circular path)
        v = np.array([
            turn_speed * math.cos(current_heading),
            turn_speed * math.sin(current_heading),
            0
        ])

        a = np.array([a_tangential * math.cos(current_heading) - a_centripetal * math.sin(current_heading),
                     a_tangential * math.sin(current_heading) + a_centripetal * math.cos(current_heading),
                     0])
        
        # Store results
        a_list.append(a)
        d += turn_speed * dt
        position_list.append((x_pos, y_pos, z_pos))
        v_list.append(v.copy())
        distance_list.append(d)
        load_factor_list.append(load_factor)
        
        # Calculate angle of attack needed for the required lift
        CL = load_factor * W / (0.5 * rho * turn_speed**2 * S)
        alpha_w = (CL - CL_0) / CL_alpha
        alpha_w_list.append(alpha_w)

def run_mission():
    phase_index.append(0)

    # Phase 1: Takeoff
    takeoff_simulation()
    print(f"Takeoff Complete at position: {position_list[-1]}")
    phase_index.append(len(time_list))
    
    # Phase 2: Climb to 30m
    climb_simulation(30)
    print(f"Climb Complete at position: {position_list[-1]}")
    phase_index.append(len(time_list))
    
    # Phase 3: Initial cruise
    cruise_simulation(-152, direction="-")
    print(f"First Cruise Complete at position: {position_list[-1]}")
    phase_index.append(len(time_list))
    
    # Phase 4: First turn (180 degrees)
    turn_simulation(180, 100, direction="right")
    print(f"First Turn Complete at position: {position_list[-1]}")
    phase_index.append(len(time_list))
    
    # Phase 5: Return cruise
    cruise_simulation(0, direction="+")
    print(f"Second Cruise Complete at position: {position_list[-1]}")
    phase_index.append(len(time_list))
    
    # Phase 6: Full loop (360 degrees)
    turn_simulation(360, 80, direction="left")
    print(f"Loop Complete at position: {position_list[-1]}")
    phase_index.append(len(time_list))
    
    # Phase 7: Outbound cruise
    cruise_simulation(152, direction="+")
    print(f"Third Cruise Complete at position: {position_list[-1]}")
    phase_index.append(len(time_list))
    
    # Phase 8: Final turn (180 degrees)
    turn_simulation(180, 100, direction="right")
    print(f"Final Turn Complete at position: {position_list[-1]}")
    phase_index.append(len(time_list))
    
    # Phase 9: Return cruise
    cruise_simulation(0, direction="-")
    print(f"Final Cruise Complete at position: {position_list[-1]}")
    phase_index.append(len(time_list))

def plot_results():
    x_coords = [pos[0] for pos in position_list]
    y_coords = [pos[1] for pos in position_list]
    z_coords = [pos[2] for pos in position_list]
    speeds = [magnitude(v) for v in v_list]

    # Separate acceleration components
    acc_x = [a[0] for a in a_list]
    acc_y = [a[1] for a in a_list]
    acc_z = [a[2] for a in a_list]
    
    plt.figure(figsize=(20, 10))

    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'purple']  # Define colors for phases

    # 3D trajectory
    ax1 = plt.subplot(221, projection='3d')
    for i in range(len(phase_index) - 1):
        start, end = phase_index[i], phase_index[i + 1]
        ax1.plot(x_coords[start:end], y_coords[start:end], z_coords[start:end], color=colors[i % len(colors)], label=f"Phase {i+1}")
    ax1.set_box_aspect([max(x_coords)-min(x_coords), max(y_coords)-min(y_coords), max(z_coords)-min(z_coords)])
    ax1.set_title('3D Flight Path')
    ax1.set_xlabel('X (m)')
    ax1.set_ylabel('Y (m)')
    ax1.set_zlabel('Z (m)')

    # Speed profile
    ax2 = plt.subplot(222)
    for i in range(len(phase_index) - 1):
        start, end = phase_index[i], phase_index[i + 1]
        ax2.plot(time_list[start:end], speeds[start:end], color=colors[i % len(colors)], label=f"Phase {i+1}")
    ax2.set_title('Speed vs Time')
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Speed (m/s)')
    ax2.grid(True)

    # AOA profile
    ax3 = plt.subplot(223)
    for i in range(len(phase_index) - 1):
        start, end = phase_index[i], phase_index[i + 1]
        ax3.plot(time_list[start:end], alpha_w_list[start:end], color=colors[i % len(colors)], label=f"Phase {i+1}")
    ax3.set_title('AOA vs Time')
    ax3.set_xlabel('Time (s)')
    ax3.set_ylabel('AOA (deg)')
    ax3.grid(True)

    # Acc profile
    ax4 = plt.subplot(224)
    for i in range(len(phase_index) - 1):
        start, end = phase_index[i], phase_index[i + 1]
        ax4.plot(time_list[start:end], acc_x[start:end], color='r', label='acc_x')  # X acceleration
        ax4.plot(time_list[start:end], acc_y[start:end], color='g', label='acc_y')  # Y acceleration
        ax4.plot(time_list[start:end], acc_z[start:end], color='b', label='acc_z')  # Z acceleration
    ax4.set_title('Acc vs Time')
    ax4.set_xlabel('Time (s)')
    ax4.set_ylabel('ACC (m/s^2)')
    ax4.grid(True)
    ax4.legend()  # Legend for acc_x, acc_y, acc_z

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    run_mission()
    plot_results()
