#
# Physics 1321 Final Project
# University of Pittsburgh
# 
# Connor McPartland
# Half-Car Model Simulation

import math
import numpy as np
import time
from scipy.misc import derivative

def initialize_shell(shell, road_choice, car_plot, com_plot, messages, y_car_front_box, y_car_rear_box, y_com_box, y_front_tire_box, y_rear_tire_box, 
                    A_box, omega_box, phi_box, **kwargs):
    shell.interact(kwargs.copy())
    clear(car_plot, com_plot, messages, y_car_front_box, y_car_rear_box, y_com_box, y_front_tire_box, y_rear_tire_box, **kwargs)
    open_inputs(road_choice, A_box, omega_box, phi_box, **kwargs)

def clear(car_plot, com_plot, messages, y_car_front_box, y_car_rear_box, y_com_box, y_front_tire_box, y_rear_tire_box, **kwargs):
    y_car_front_box.value = '0'
    y_car_rear_box.value = '0'
    y_front_tire_box.value = '0'
    y_rear_tire_box.value = '0'
    y_com_box.value = '0'
    car_plot.clear()
    com_plot.clear()
    messages.clear()

n_depvars = 8
f_return = np.zeros((n_depvars), dtype=np.float)
y_temp = np.zeros((n_depvars), dtype=np.float)
k1 = np.zeros((n_depvars), dtype=np.float)
k2 = np.zeros((n_depvars), dtype=np.float)
k3 = np.zeros((n_depvars), dtype=np.float)
k4 = np.zeros((n_depvars), dtype=np.float)
t=0.

def RK4(f, dt, t, y):
    f(t, y, f_return)
    k1[:] = dt*f_return 
    y_temp[:] = y + k1/2. 
    f(t + dt/2., y_temp, f_return) 
    k2[:] = dt*f_return 
    y_temp[:] = y + k2/2. 
    f(t + dt/2., y_temp, f_return)
    k3[:] = dt*f_return  
    y_temp[:] = y + k3 
    f(t+dt, y_temp, f_return) 
    k4[:] = dt*f_return  
    y += (k1+2.*(k2+k3)+k4)/6.

def open_inputs(road_choice, A_box, omega_box, phi_box, **kwargs):
    if road_choice.value == 'Sinusoid (Asin(wx + phi))' or road_choice.value=='Square Wave':
        A_box.enabled=True
        omega_box.enabled=True
        phi_box.enabled=True
    else:
        A_box.enabled=False
        omega_box.enabled=False
        phi_box.enabled=False
    
def run(v_0_box, road_choice, A_box, omega_box, phi_box, height_box,
        k_fs_box, k_rs_box, b_fs_box, b_rs_box, fl_rs_box, fl_fs_box, 
        r_R_box, r_F_box, m_rt_box, m_ft_box, k_ft_box, k_rt_box, b_rt_box, b_ft_box, 
        y_car_front_box, y_car_rear_box, y_com_box, y_front_tire_box, y_rear_tire_box,
        Lf_box, Lr_box, m_c_box, com_plot, car_plot, plot_com_choice, stop, messages, **kwargs):

    # Define Suspension Variables
    k_fs = k_fs_box.value       # spring constant of the Front Spring
    k_rs = k_rs_box.value       # spring constant of the Rear Spring
    fl_rs = fl_rs_box.value     # free length of Rear Spring
    fl_fs = fl_rs_box.value     # free length of Front Spring
    b_fs = b_fs_box.value       # damping coefficient of the Front Damper
    b_rs = b_rs_box.value       # damping coefficient of the Rear Damper
    m_ft = m_ft_box.value       # mass of the Front Tire
    m_rt = m_rt_box.value       # mass of the Rear Tire
    k_ft = k_ft_box.value       # spring constant of the Front Tire
    k_rt = k_rt_box.value       # spring constant of the Rear Tire
    r_R = r_R_box.value         # radius of the Rear Tire
    r_F = r_F_box.value         # radius of the Front Tire
    m_rt = m_rt_box.value       # mass of the Rear Tire in kg
    m_ft = m_ft_box.value       # mass of the Front Tire in kg
    b_ft = b_ft_box.value       # damping coefficient of the Front Damper
    b_rt = b_rt_box.value       # damping coefficient of the Front Damper
    Lf = Lf_box.value           # the distance from the c.o.m. to the front suspension
    Lr = Lr_box.value           # the distance from the c.o.m. to the rear suspension          
    m_c = m_c_box.value         # the mass of the vehicle, located at the c.o.m.
        
    L = Lf+Lr
    g = 9.81
    Ic = m_c*(1./3.)*(Lr*Lf**2 + Lf*Lr**2)/L # the moment of the inertia of the vehicle
    
    fl_rs += m_rt*g/k_rs        # add the extenions of the springs due to the weight of the tires
    fl_fs += m_ft*g/k_fs
    

    # Define the profile of the road.
    amp = float(A_box.value)
    w = float(omega_box.value)
    phi = float(phi_box.value)
    if road_choice.value == 'Sinusoid (Asin(wx + phi))':
        def y_road(x):
            return amp*np.sin(w*x + phi)
        def v_road(x):
            return amp*w*np.cos(w*x + phi)
    if road_choice.value == 'Flat Surface':
        def y_road(x):
            return 0*x
        def v_road(x):
            return 0*x
    if road_choice.value=='Square Wave':
        def y_road(x):
            return amp*(np.int_(x*w/5)%2==0)
        def v_road(x):
            return derivative(y_road, x)

    # Set up plotting
    car_plot.set_plot_properties(
        title='Half-Car Model Simulated Motion',
        x_label='x (m)',
        y_label='y (m)',
        x_scale='linear',
        y_scale='linear',
        aspect_ratio='auto')
    car_plot.new_curve('tire_rear', memory='growable', animated=False, 
                     line_style='-', line_width=0.5, line_color='black',
                     marker_style='o', marker_color='grey', marker_edge_color='black', marker_width=25, marker_edge_width=5)
    car_plot.new_curve('tire_front', memory='growable', animated=False, 
                     line_style='-', line_width=0.5, line_color='black',
                     marker_style='o', marker_color='grey', marker_edge_color='black', marker_width=25, marker_edge_width=5)
    car_plot.new_curve('car_outline', memory='growable', animated=False, 
                     line_style='-', line_width=5.0, line_color='red',
                     marker_style='o', marker_color='black', marker_width=5.0)
    car_plot.new_curve('road', memory='growable', animated=False, 
                     line_style='-', line_width=0.5, line_color='black', marker_width=3.0)
                     
    com_plot.set_plot_properties(
        title='Position of Center of Mass',
        x_label = 't',
        y_label = 'Height (m)',
        aspect_ratio='auto')
    com_plot.new_curve('com', memory='growable', animated=False, 
                     line_style='-', line_width=1.0, line_color='green', marker_width=3.0)
                     
    # Find initial conditions by solving Ay = b
    A = np.array([[-k_fs, -k_rs, k_fs, k_rs], 
                  [-k_fs*L, 0, k_fs*L, 0],
                  [0, k_rs, 0, -(k_rt+k_rs)],
                  [k_fs, 0, -(k_ft+k_fs), 0]])
    b = np.array([m_c*g - k_rs*fl_rs - k_fs*fl_fs,
                  m_c*g*Lr - k_fs*L*fl_fs,
                  m_rt*g - k_rt*r_R + k_rs*fl_rs - k_rt*y_road(0),
                  m_ft*g - k_ft*r_F + k_fs*fl_fs - k_ft*y_road(L)])
    temp = np.linalg.solve(A,b)
    y_car_front0 = temp[0]
    y_car_rear0 = temp[1]
    y_tire_front0 = temp[2]
    y_tire_rear0 = temp[3]
    theta0 = math.asin((y_car_front0-y_car_rear0)/L)
    messages.write('Initial Equilibrium Positions\n \tFront Car: %.5f\n \tRear Car: %.5f\n \tFront Tire: %.5f\n \tRear Tire: %.5f\n' % tuple(temp))
    initial_height = height_box.value
    y = np.array([y_car_front0 - Lf*math.sin(theta0) + initial_height, 0., theta0, 0, y_tire_front0 + initial_height, 0., y_tire_rear0+initial_height, 0.])
    v_0 = v_0_box.value*1000/60./60. # convert to m/s

    def f(t, y, f_return):
        y_com = y[0]
        v_com = y[1]
        theta = y[2]
        omega = y[3]
        y_tire_front = y[4]
        v_tire_front = y[5]
        y_tire_rear = y[6]
        v_tire_rear = y[7]
        
        delta_fs = (y_com + Lf*math.sin(theta)) - y_tire_front - fl_fs
        delta_rs = (y_com - Lr*math.sin(theta))- y_tire_rear - fl_rs
        # The clipping cuts off interaction between the road and the tire if the tire goes airborne.
        delta_ft = np.clip(np.array([y_tire_front - y_road(v_0*t + L) - r_F]), a_min=None, a_max=0)[0]
        delta_rt = np.clip(np.array([y_tire_rear - y_road(v_0*t) - r_R]), a_min=None, a_max=0)[0]
        v_fs = Lf*omega + v_com - v_tire_front 
        v_rs = -Lr*omega + v_com - v_tire_rear
        if delta_ft == 0:
            v_ft = 0
        else:
            v_ft = v_tire_front - v_road(v_0*t+L)
        if delta_rt == 0:
            v_rt = 0
        else:
            v_rt = v_tire_rear - v_road(v_0*t)
        
        # The clipping and if-cases basically make it so that if the car tire is airborn, the tire-spring and tire-damper produce no forces on the tire. 
        f_return[0] = v_com
        f_return[1] = -g - (k_fs/m_c)*delta_fs - (b_fs/m_c)*v_fs - (k_rs/m_c)*delta_rs - (b_rs/m_c)*v_rs
        f_return[2] = omega
        f_return[3] = (-Lf*math.cos(theta)*(k_fs*delta_fs + b_fs*v_fs) + Lr*math.cos(theta)*(k_rs*delta_rs + b_rs*v_rs))/Ic 
        f_return[4] = v_tire_front
        f_return[5] = -g - (k_ft/m_ft)*delta_ft - (b_ft/m_ft)*v_ft + (k_fs/m_ft)*delta_fs + (b_fs/m_ft)*v_fs
        f_return[6] = v_tire_rear
        f_return[7] = -g - (k_rt/m_rt)*delta_rt - (b_rt/m_rt)*v_rt + (k_rs/m_rt)*delta_rs + (b_rs/m_rt)*v_rs
    
    t = 0.01
    dt = .002
    xs = np.arange(-5, 5*v_0/dt, v_0*dt)
    road = y_road(xs)
    n=0
    
    y_com = y[0]
    theta = y[2]        
    y_tire_front = y[4]
    y_tire_rear = y[6]
    y_car_front = y_com + Lf*math.sin(theta)
    y_car_rear = y_com - Lr*math.sin(theta)
    car_plot.set_plot_properties(x_limits=(v_0*t-.25*L, v_0*(t+2) + 2*L), y_limits=(-A_box.value, 4))
    car_plot.set_data('road', np.column_stack((xs, road)), redraw=False)
    car_plot.set_data('car_outline', [[v_0*t+L, y_tire_front], [v_0*t + L, y_car_front], [v_0*t+Lr, y_com],  [v_0*t, y_car_rear], [v_0*t, y_tire_rear]], redraw=False)
    car_plot.set_data('tire_front', [[v_0*t + L, y_tire_front]], redraw=False)
    car_plot.set_data('tire_rear', [[v_0*t, y_tire_rear]], redraw=True)
    y_car_front_box.value = '%.5f' % y_car_front
    y_car_rear_box.value = '%.5f' % y_car_rear
    y_front_tire_box.value = '%.5f' % y_tire_front
    y_rear_tire_box.value = '%.5f' % y_tire_rear
    y_com_box.value = '%.5f' % y_com
    if plot_com_choice.value:
        com_plot.append_data('com', (t, y_com))
    while True:
        # Rename variables for convenience.
        y_com = y[0]
        theta = y[2]        
        y_tire_front = y[4]
        y_tire_rear = y[6]
        y_car_front = y_com + Lf*math.sin(theta)
        y_car_rear = y_com - Lr*math.sin(theta)
        
        if n%8 == 0:
            car_plot.set_data('car_outline', [[v_0*t+L, y_tire_front], [v_0*t + L, y_car_front], [v_0*t+Lr, y_com],  [v_0*t, y_car_rear], [v_0*t, y_tire_rear]], redraw=False)
            car_plot.set_data('tire_front', [[v_0*t + L, y_tire_front]], redraw=False)
            car_plot.set_data('tire_rear', [[v_0*t, y_tire_rear]], redraw=True)
            y_car_front_box.value = '%.5f' % y_car_front
            y_car_rear_box.value = '%.5f' % y_car_rear
            y_front_tire_box.value = '%.5f' % y_tire_front
            y_rear_tire_box.value = '%.5f' % y_tire_rear
            y_com_box.value = '%.5f' % y_com
            if plot_com_choice.value:
                com_plot.append_data('com', (t, y_com))
        if int(v_0*t+1) % int(1.9*v_0+2*L) == 0:
            car_plot.set_plot_properties(x_limits=(v_0*t-.25*L, v_0*(t+2) + 2*L))
        
        n+=1
        RK4(f,dt,t,y)        
        t+=dt
        if stop.value: break
    stop.value=False
  
   
