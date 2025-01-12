import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

### Constants ###
rho = 1.2  # air density (kg/m^3)
g = 9.81  # gravity (m/s^2)
m_glider = 4 # glider mass (kg)
m_payload = 3  # payload mass (kg)
m_x1 = 0.2  # additional mass (kg)
W = (m_glider + m_payload + m_x1) * g  # total weight (N)
m = m_glider + m_payload + m_x1  # total mass (kg)
S = 0.6  # wing area (m^2)
AR = 7.2  # aspect ratio
lw = 0.2
lh = 1
CD0 = 0.1  # zero-lift drag coefficient
CL0 = 0.0  # lift coefficient at zero angle of attack , OpenVSP 결과과
CL_alpha = 0.086  # lift coefficient gradient per degree , OpenVSP 결과
e = 0.8  # Oswald efficiency factor
T_max = 6.6 * g  # maximum thrust (N)
alpha_stall = 10
V_stall = math.sqrt(W/(alpha_stall * CL_alpha * 0.5 * rho * S))
print(V_stall)
h_flap = 5

### Helper Functions ###
def magnitude(vector):
    return math.sqrt(sum(x*x for x in vector))

def calculate_induced_drag(C_L):
    return (C_L**2) / (math.pi * AR * e)

def calculate_cruise_alpha_w(v):
    speed = magnitude(v)
    def equation(alpha_w):
        CL = CL0 + CL_alpha * alpha_w
        L = 0.5 * rho * speed**2 * S * CL
        return L + T_max * 0.9 * math.sin(math.radians(alpha_w)) - W
    alpha_w_solution = fsolve(equation, 5, xtol=1e-8, maxfev=1000)
    return alpha_w_solution[0]

### Result Lists ###
time_list = []
distance_list = []
load_factor_list = []
AOA_list = []
position_list = []
v_list = []
a_list = []
phase_index = []
bank_angle_list = []

### Acceleration Functions ###
def calculate_acceleration_ground(v):
    speed = magnitude(v)
    CL = 0.99
    CDi = calculate_induced_drag(CL)
    CD = CD0 + CDi
    D = 0.5 * rho * speed**2 * S * CD
    L = 0.5 * rho * speed**2 * S * CL
    a_x = -g/W * (T_max*0.9 - D - 0.03*(W-L))
    return np.array([a_x, 0, 0])

def calculate_acceleration_cruise(v, alpha_w_deg):
    speed = magnitude(v)
    CL = CL0 + CL_alpha * alpha_w_deg
    CDi = calculate_induced_drag(CL)
    CD = CD0 + CDi
    D = 0.5 * rho * speed**2 * S * CD
    a_x = (T_max * 0.9) / m * math.cos(math.radians(alpha_w_deg)) - D / m
    return np.array([a_x, 0, 0])

def calculate_acceleration_climb(v, alpha_w_deg, gamma_rad, theta_deg, z_pos):
    speed = magnitude(v)
    # if (z_pos > h_flap) :
    #     CL = CL0 + CL_alpha * alpha_w_deg
    # else:
    #     CL = 0.99 + 0.06 * alpha_w_deg
    CL = CL0 + CL_alpha * alpha_w_deg
    CDi = calculate_induced_drag(CL)
    CD = CD0 + CDi
    D = 0.5 * rho * speed**2 * S * CD
    L = 0.5 * rho * speed**2 * S * CL
    a_x = -(T_max * 0.9 * math.cos(math.radians(theta_deg)) - L * math.sin(math.radians(alpha_w_deg)) - D * math.cos(gamma_rad)) / m
    a_z = (T_max * 0.9 * math.sin(math.radians(theta_deg)) + L * math.cos(math.radians(alpha_w_deg)) - D * math.sin(gamma_rad) - W) / m
    return np.array([a_x, 0, a_z])

### Simulation Functions ###
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
        
        L = 0.5 * rho * magnitude(v)**2 * S * CL0
        load_factor_list.append(L / W)
        v_list.append(v.copy())
        AOA_list.append(0)
        a_list.append(a)
        position_list.append(tuple(position))
        bank_angle_list.append(math.degrees(0))

def climb_simulation(h_max):
    print("\nRunning Climb Simulation...")
    
    dt = 0.01
    n_steps = int(60 / dt)  # Max 60 seconds simulation
    v = v_list[-1].copy()
    d = distance_list[-1] if distance_list else 0
    t = time_list[-1]
    x_pos, y_pos, z_pos = position_list[-1]
    theta_deg = 0
    bank_angle_list.append(math.degrees(0))

    for step in range(n_steps):
        t += dt
        time_list.append(t)

        if (abs(z_pos - h_max) < 10):
            theta_deg -= 0.5
            theta_deg = max(theta_deg, 10)
        else:
            theta_deg += 0.5
        theta_deg = min(theta_deg, 40) 
        if (theta_deg == 50): print(step)

        # Calculate climb angle
        gamma_rad = math.atan2(abs(v[2]), abs(v[0]))
        alpha_w_deg = theta_deg - math.degrees(gamma_rad)

        # Calculate load factor
        # if (z_pos > h_flap) :
        #     CL = CL0 + CL_alpha * alpha_w_deg
        # else:
        #     CL = 0.99 + 0.06 * alpha_w_deg
        CL = CL0 + CL_alpha * alpha_w_deg
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
        AOA_list.append(alpha_w_deg)
        a_list.append(a)
        distance_list.append(d)
        bank_angle_list.append(math.degrees(0))

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
        CL = CL0 + CL_alpha * alpha_w_deg
        L = 0.5 * rho * speed**2 * S * CL
        
        # Store results
        load_factor_list.append(L / W)
        v_list.append(v.copy())
        AOA_list.append(alpha_w_deg)
        a_list.append(a)
        distance_list.append(d)
        bank_angle_list.append(math.degrees(0))
        
        # Check if we've reached target x position
        if direction == '+':
            if x_pos >= x_final:
                break
        else:
            if x_pos <= x_final:
                break

def turn_simulation(target_angle_deg, direction="right"):
    print("\nRunning Turn Simulation...")
    
    # 초기 설정
    dt = 0.01
    v = v_list[-1].copy()
    d = distance_list[-1] if distance_list else 0
    t = time_list[-1]
    x_pos, y_pos, z_pos = position_list[-1]
    speed = magnitude(v)

    # Initialize turn tracking
    target_angle_rad = math.radians(target_angle_deg)
    turned_angle_rad = 0

    # Get initial heading and setup turn center
    initial_angle_rad = math.atan2(v[1], v[0])
    current_angle_rad = initial_angle_rad

    step = 0

    # Turn
    while abs(turned_angle_rad) < abs(target_angle_rad):
        t += dt
        time_list.append(t)
        
        CL = CL0 + CL_alpha * alpha_stall
        L = CL * (0.5 * rho * speed**2) * S
        phi_rad = math.acos(W/L)
        a_centripetal = (L * math.sin(phi_rad)) / m
        R = (m * speed**2)/(L * math.sin(phi_rad))
        omega = speed / R
        load_factor = 1 / math.cos(phi_rad)

        CD = CD0 + calculate_induced_drag(CL)
        D = CD * (0.5 * rho * speed**2) * S
        a_tangential = (T_max * 0.55 - D) / m
        speed += a_tangential * dt

        # Calculate turn center
        if direction == "right":
            center_x = x_pos - R * math.sin(current_angle_rad)
            center_y = y_pos + R * math.cos(current_angle_rad)
        else:
            center_x = x_pos + R * math.sin(current_angle_rad)
            center_y = y_pos - R * math.cos(current_angle_rad)

        # Update heading based on angular velocity
        if direction == "right":
            current_angle_rad += omega * dt
            turned_angle_rad += omega * dt
        else:
            current_angle_rad -= omega * dt
            turned_angle_rad -= omega * dt
        
        # Calculate new position relative to turn center
        if direction == "right":
            x_pos = center_x + R * math.sin(current_angle_rad)
            y_pos = center_y - R * math.cos(current_angle_rad)
        else:
            x_pos = center_x - R * math.sin(current_angle_rad)
            y_pos = center_y + R * math.cos(current_angle_rad)

        # Update velocity direction (tangent to the circular path)
        v = np.array([
            speed * math.cos(current_angle_rad),
            speed * math.sin(current_angle_rad),
            0
        ])

        a = np.array([a_tangential * math.cos(current_angle_rad) - a_centripetal * math.sin(current_angle_rad),
                     a_tangential * math.sin(current_angle_rad) + a_centripetal * math.cos(current_angle_rad),
                     0])
        
        # Store results
        a_list.append(a)
        d += speed * dt
        position_list.append((x_pos, y_pos, z_pos))
        v_list.append(v.copy())
        distance_list.append(d)
        load_factor_list.append(load_factor)
        AOA_list.append(alpha_stall)
        bank_angle_list.append(math.degrees(phi_rad))

        # if (step%500 == 0):
        #     print(f"CL: {CL:.2f}")
        #     print(f"L: {L:.2f} N")
        #     print(f"phi_deg: {math.degrees(phi_rad):.2f} deg")
        #     # print(f"radius: {R:.2f} m")
        #     print(f"speed: {speed:.2f} m/s\n")
        #     # print(f"omega: {omega:.2f} rad/s")

### Mission Function & Plotting ###
def run_mission3():
    phase_index.append(0)

    # Phase 1: Takeoff
    takeoff_simulation()
    print(f"Takeoff Complete at position: {position_list[-1]}")
    phase_index.append(len(time_list))
    
    # Phase 2: Climb to 30m
    climb_simulation(70)
    print(f"Climb Complete at position: {position_list[-1]}")
    phase_index.append(len(time_list))
    
    # Phase 3: Initial cruise
    cruise_simulation(-152, direction="-")
    print(f"First Cruise Complete at position: {position_list[-1]}")
    phase_index.append(len(time_list))
    
    # Phase 4: First turn (180 degrees)
    turn_simulation(180, direction="right")
    print(f"First Turn Complete at position: {position_list[-1]}")
    phase_index.append(len(time_list))
    
    # Phase 5: Return cruise
    cruise_simulation(0, direction="+")
    print(f"Second Cruise Complete at position: {position_list[-1]}")
    phase_index.append(len(time_list))
    
    # Phase 6: Full loop (360 degrees)
    turn_simulation(360, direction="left")
    print(f"Loop Complete at position: {position_list[-1]}")
    phase_index.append(len(time_list))
    
    # Phase 7: Outbound cruise
    cruise_simulation(152, direction="+")
    print(f"Third Cruise Complete at position: {position_list[-1]}")
    phase_index.append(len(time_list))
    
    # Phase 8: Final turn (180 degrees)
    turn_simulation(180, direction="right")
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
    
    plt.figure(figsize=(20, 10))

    gridspec = plt.GridSpec(2, 3, width_ratios=[1, 1, 1], height_ratios=[1, 1])

    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'purple']  # Define colors for phases

    # 3D trajectory
    ax1 = plt.subplot(gridspec[:, 0], projection='3d')
    for i in range(len(phase_index) - 1):
        start, end = phase_index[i], phase_index[i + 1]
        ax1.plot(x_coords[start:end], y_coords[start:end], z_coords[start:end], color=colors[i % len(colors)], label=f"Phase {i+1}")

    # Get current axis limits
    x_limits = ax1.get_xlim()
    y_limits = ax1.get_ylim()
    z_limits = ax1.get_zlim()

    # Find the max range for uniform scaling
    x_range = x_limits[1] - x_limits[0]
    y_range = y_limits[1] - y_limits[0]
    z_range = z_limits[1] - z_limits[0]
    max_range = max(x_range, y_range, z_range)

    # Set the new limits
    mid_x = 0.5 * (x_limits[0] + x_limits[1])
    mid_y = 0.5 * (y_limits[0] + y_limits[1])
    mid_z = 0.5 * (z_limits[0] + z_limits[1])
    ax1.set_xlim(mid_x - max_range / 2, mid_x + max_range / 2)
    ax1.set_ylim(mid_y - max_range / 2, mid_y + max_range / 2)
    ax1.set_zlim(mid_z - max_range / 2, mid_z + max_range / 2)

    ax1.set_title('3D Flight Path')
    ax1.set_xlabel('X (m)')
    ax1.set_ylabel('Y (m)')
    ax1.set_zlabel('Z (m)')


    # Speed profile
    ax2 = plt.subplot(gridspec[0, 1])
    for i in range(len(phase_index) - 1):
        start, end = phase_index[i], phase_index[i + 1]
        ax2.plot(time_list[start:end], speeds[start:end], color=colors[i % len(colors)], label=f"Phase {i+1}")
    ax2.set_title('Speed vs Time')
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Speed (m/s)')
    ax2.grid(True)

    # AOA profile
    ax3 = plt.subplot(gridspec[1, 1])
    for i in range(len(phase_index) - 1):
        start, end = phase_index[i], phase_index[i + 1]
        ax3.plot(time_list[start:end], AOA_list[start:end], color=colors[i % len(colors)], label=f"Phase {i+1}")
    ax3.set_title('AOA vs Time')
    ax3.set_xlabel('Time (s)')
    ax3.set_ylabel('AOA (deg)')
    ax3.grid(True)

    # Bank angle profile
    ax4 = plt.subplot(gridspec[0, 2])
    for i in range(len(phase_index) - 1):
        start, end = phase_index[i], phase_index[i + 1]
        ax4.plot(time_list[start:end], bank_angle_list[start:end], color=colors[i % len(colors)], label=f"Phase {i+1}")
    ax4.set_title('Bank Angle vs Time')
    ax4.set_xlabel('Time (s)')
    ax4.set_ylabel('Bank Angle (deg)')
    ax4.grid(True)

    # Load Factor profile
    ax5 = plt.subplot(gridspec[1, 2])
    for i in range(len(phase_index) - 1):
        start, end = phase_index[i], phase_index[i + 1]
        ax5.plot(time_list[start:end], load_factor_list[start:end], color=colors[i % len(colors)], label=f"Phase {i+1}")
    ax5.set_title('Load Factor vs Time')
    ax5.set_xlabel('Time (s)')
    ax5.set_ylabel('Load Factor')
    ax5.grid(True)

    plt.tight_layout()
    plt.show()

def save_results():
    import os

    # Create a result directory if it does not exist
    results = "results"
    if not os.path.exists(results):
        os.makedirs(results)

    # Save data to a .npz file
    np.savez(
        os.path.join(results, "mission3.npz"), 
        time_list=time_list, 
        distance_list=distance_list, 
        load_factor_list=load_factor_list, 
        AOA_list=AOA_list, 
        position_list=position_list, 
        v_list=v_list, 
        a_list=a_list, 
        phase_index=phase_index,
        bank_angle_list=bank_angle_list
    )

    print(f"\nData saved to {os.path.join(results, 'mission3.npz')}\n")

if __name__ == "__main__":
    run_mission3()
    plot_results()
    # save_results()