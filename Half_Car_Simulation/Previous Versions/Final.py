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

n_depvars = 14
f_return = np.zeros((n_depvars), dtype=np.float)
y_temp = np.zeros((n_depvars), dtype=np.float);
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
        k_frs_box, k_fls_box, k_rrs_box, k_rls_box, b_frs_box, b_fls_box, b_rrs_box, b_rls_box, fl_frs_box, fl_fls_box, fl_rrs_box, fl_rls_box,
        r_Fr_box, r_Fl_box, r_R_box,  m_frt_box, m_flt_box, m_rrt_box, m_rlt_box, 
        k_frt_box, k_flt_box, k_rrt_box, k_rlt_box, b_frt_box, b_flt_box, b_rrt_box, b_rlt_box, 
        y_car_front_box, y_car_rear_box, y_com_box, y_front_tire_box, y_rear_tire_box,
        Lf_box, Lr_box, Wr_box, Wl_box, m_c_box, com_plot, car_plot, plot_com_choice, stop, messages, **kwargs):

    # Define Suspension Variables
    k_frs = k_frs_box.value       # spring constant of the Front Right Spring
    k_fls = k_fls_box.value       # spring constant of the Front Left Spring
    k_rrs = k_rrs_box.value       # spring constant of the Rear Right Spring
    k_rls = k_rls_box.value       # spring constant of the Rear Left Spring
    fl_frs = fl_frs_box.value     # free length of Front Right Spring
    fl_fls = fl_fls_box.value     # free length of Front Left Spring
    fl_rrs = fl_rrs_box.value     # free length of Right Right Spring
    fl_rls = fl_rls_box.value     # free length of Right Left Spring
    b_frs = b_frs_box.value       # damping coefficient of the Front Right Damper
    b_fls = b_fls_box.value       # damping coefficient of the Front Left Damper
    b_rrs = b_rrs_box.value       # damping coefficient of the Rear Right Damper
    b_rls = b_rls_box.value       # damping coefficient of the Rear Left Damper
    m_frt = m_frt_box.value       # mass of the Front Right Tire
    m_flt = m_flt_box.value       # mass of the Front Left Tire
    m_rrt = m_rrt_box.value       # mass of the Rear Right Tire
    m_rlt = m_rlt_box.value       # mass of the Rear Left Tire
    k_frt = k_frt_box.value       # spring constant of the Front Right Tire
    k_flt = k_flt_box.value       # spring constant of the Front Left Tire
    k_rrt = k_rrt_box.value       # spring constant of the Rear Right Tire
    k_rlt = k_rlt_box.value       # spring constant of the Rear Left Tire
    r_Fr = r_Fr_box.value         # radius of the Front Right Tire = free length of the Front Right Tire
    r_Fl = r_Fl_box.value         # radius of the Front Left Tire = free length of the Front Left Tire
    r_Rr = r_Rr_box.value         # radius of the Rear Right Tire = free length of the Rear Right Tire
    r_Rl = r_Rl_box.value         # radius of the Rear Left Tire = free length of the Rear Left Tire
    m_frt = m_frt_box.value       # mass of the Front Right Tire in kg    
    m_flt = m_flt_box.value       # mass of the Front Left Tire in kg    
    m_rrt = m_rrt_box.value       # mass of the Rear Right Tire in kg
    m_rlt = m_rlt_box.value       # mass of the Rear Left Tire in kg

    b_frt = b_frt_box.value       # damping coefficient of the Front Right Tire
    b_flt = b_flt_box.value       # damping coefficient of the Front Left Tire
    b_rrt = b_rrt_box.value       # damping coefficient of the Rear Right Tire
    b_rlt = b_rlt_box.value       # damping coefficient of the Rear Left Tire
    Lfront = Lf_box.value           # the distance from the c.o.m. to the front suspension
    Lrear = Lr_box.value           # the distance from the c.o.m. to the rear suspension   
    Wright = Wr_box.value           # the distance from c.o.m. to right of car
    Wleft = Wl_box.value            # the distance from c.o.m. to left of car
    m_c = m_c_box.value         # the mass of the vehicle, located at the c.o.m.
        
    L = Lfront+Lrear
    W = Wright+Wleft
    g = 9.81
    Icl = m_c*(1./3.)*(Lrear*Lfront**2 + Lfront*Lrear**2)/L # the moment of the inertia of the vehicle
    Icw = m_c*(1./3.)*(Wright*Wleft**2 + Wleft*Wright**2)/W
    #
    #fl_rs += m_rt*g/k_rs        # add the extenions of the springs due to the weight of the tires
    #fl_fs += m_ft*g/k_fs
    

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
        theta1 = y[2]
        omega1 = y[3]
        theta2 = y[4]
        omega2 = y[5]
        z_frt = y[6]
        v_frt = y[7]
        z_flt = y[8]
        v_flt = y[9]
        z_rrt = y[10]
        v_rrt = y[11]
        z_rlt = y[12]
        v_rlt = y[13]
        
        delta_frs = (y_com + Lfront*math.sin(theta1) + Wright*math.sin(theta2)) - z_frt- fl_frs
        delta_fls = (y_com + Lfront*math.sin(theta1) - Wleft*math.sin(theta2)) - z_flt- fl_fls
        delta_rrs = (y_com - Lrear*math.sin(theta) + Wright*math.sin(theta2))- z_rrt - fl_rrs
        delta_rls = (y_com - Lrear*math.sin(theta) - Wleft*math.sin(theta2))- z_rlt - fl_rls
        # The clipping cuts off interaction between the road and the tire if the tire goes airborne.
        delta_frt = np.clip(np.array([z_frt - y_road(v_0*t + L) - r_Fr]), a_min=None, a_max=0)[0]
        delta_flt = np.clip(np.array([z_flt - y_road(v_0*t + L) - r_Fl]), a_min=None, a_max=0)[0]
        delta_rrt = np.clip(np.array([z_rrt - y_road(v_0*t + L) - r_Rr]), a_min=None, a_max=0)[0]
        delta_rlt = np.clip(np.array([z_rlt - y_road(v_0*t + L) - r_Rl]), a_min=None, a_max=0)[0]
        
        v_frs = Lfront*omega1 + Wright*omega2 + v_com - v_frt
        v_fls = Lfront*omega1 - Wright*omega2 + v_com - v_flt
        v_rrs = -Lr*omega + Wright*omega2 + v_com - v_rrt
        v_rls = -Lr*omega - Wright*omega2 + v_com - v_rlt
        if delta_frt == 0:
            dv_frt = 0
        else:
            dv_frt = v_frt - v_road(v_0*t+L)
        if delta_flt == 0:
            dv_flt = 0
        else:
            dv_flt = v_flt - v_road(v_0*t+L)
        if delta_rrt == 0:
            dv_rrt = 0
        else:
            dv_rrt = v_rrt - v_road(v_0*t)
        if delta_rlt == 0:
            dv_rlt = 0
        else:
            dv_rlt = v_rlt - v_road(v_0*t)
        
        # The clipping and if-cases basically make it so that if the car tire is airborn, the tire-spring and tire-damper produce no forces on the tire. 
        f_return[0] = v_com
        f_return[1] = -g - (k_frs/m_c)*delta_frs - (b_frs/m_c)*v_frs - (k_fls/m_c)*delta_fls - (b_fls/m_c)*v_fls - (k_rrs/m_c)*delta_rrs - (b_rrs/m_c)*v_rrs- (k_rls/m_c)*delta_rls - (b_rls/m_c)*v_rls
        f_return[2] = omega1
        f_return[3] = (-Lfront*math.cos(theta1)*(k_frs*delta_frs + b_frs*v_frs) + (-Lfront*math.cos(theta1)*(k_fls*delta_fls + b_fls*v_fls)  
                        + Lrear*math.cos(theta1)*(k_rrs*delta_rrs + b_rrs*v_rrs)+ Lrear*math.cos(theta1)*(k_rls*delta_rls + b_rls*v_rls))/Icl
        f_return[4] = omega2
        f_return[5] = (-Wright*math.cos(theta2)*(k_frs*delta_frs + b_frs*v_frs) + (-Wleft*math.cos(theta2)*(k_fls*delta_fls + b_fls*v_fls)  
                        + Wright*math.cos(theta2)*(k_rrs*delta_rrs + b_rrs*v_rrs)+ Wleft*math.cos(theta2)*(k_rls*delta_rls + b_rls*v_rls))/Icw
        f_return[6] = v_frt
        f_return[7] = -g - (k_frt/m_frt)*delta_frt - (b_frt/m_frt)*dv_frt + (k_frs/m_frt)*delta_frs + (b_frs/m_frt)*v_frs
        f_return[8] = v_flt
        f_return[9] = -g - (k_flt/m_flt)*delta_flt - (b_flt/m_flt)*dv_flt + (k_fls/m_frt)*delta_fls + (b_frs/m_flt)*v_fls
        f_return[10] = v_rrt
        f_return[11] = -g - (k_rrt/m_rrt)*delta_rrt - (b_rrt/m_rrt)*dv_rrt + (k_rrs/m_rrt)*delta_rrs + (b_rrs/m_rrt)*v_rrs
        f_return[12] = v_rlt
        f_return[13] = -g - (k_rlt/m_rlt)*delta_rlt - (b_rlt/m_rlt)*dv_rrt + (k_rls/m_rlt)*delta_rls + (b_rls/m_rlt)*v_rls
    
    t = 0.01
    dt = .002
    xs = np.arange(-5, v_0/dt, v_0*dt)
    road = y_road(xs)
    n=0
    
    z_com = z[0]
    theta1 = z[2] 
    theta2 = z[4]
    z_frt = z[6]
    z_flt = z[8]
    z_rrt = z[10]
    z_rlt = z[12]
    z_frc = z_com + Lfront*math.sin(theta1) + Wright*math.sin(theta2)
    z_flc = z_com + Lfront*math.sin(theta1) - WLeft*math.sin(theta2)
    z_rrc = z_com - Lrear*math.sin(theta1) + Wright*math.sin(theta2)
    z_rlc = z_com - Lrear*math.sin(theta1) - Wleft*math.sin(theta2)
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
  
   
