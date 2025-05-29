from vehicle_dynamics import *

gvw = 1000
gear_ratio = 12
wheel_radius = 0.27
frontal_area=2.5
cd = 0.45
rolling_resistance_coeff=0.015
gear_efficiency=95
rated_speed_kmph=30
continuous_gradient=5
max_speed=60
gradeability=12.4
time_to_cross_grade_from_rest=12
length_of_grade=20
acceleration_from_rest_to_speed=50 
time_to_accelerate=12
v_dc=48


vd = VehicleDynamics(gvw, gear_ratio, wheel_radius, frontal_area, cd, rolling_resistance_coeff, gear_efficiency, rated_speed_kmph, continuous_gradient, max_speed, gradeability, time_to_cross_grade_from_rest, length_of_grade, acceleration_from_rest_to_speed, time_to_accelerate, v_dc)

vd.run(verbose=False)

a = [0,	18,	25,	30,	32,	42,	52,	65,	90,	108]
b=  [0, 0,	25,	15,	15,	33,	20,	42,	42,	0]




c,d,e = vd.calculate_drive_cycle(a, b)

print(c,d*1000,e)

print(vd.calculate_equivalent_operating_points(c, d*1000, e))

