
# process
x4=[20.892,0.267376,198.898,0.498318,16.7562] # AM
x1=[9.04191, 0.122593,83.2673,0.197334,6.92663] # AR
x3=[26.4805,0.321682,313.643,0.699551,25.4554] # PS
x2=[17.6139,0.237026,84.5844,0.217731,8.18291] # SS

base=[8.32518,0.113716,83.2135, 0.196499,6.81607]
print([ i-j for i,j in zip(x4,base)])
print([ i-j for i,j in zip(x1,base)])
print([ i-j for i,j in zip(x3,base)])
print([ i-j for i,j in zip(x2,base)])

12.56682, 0.15366000000000002, 115.6845, 0.30181899999999995, 9.94013
0.7167300000000001, 0.008876999999999996, 0.05380000000000962, 0.0008350000000000024, 0.11056000000000044
18.15532, 0.20796600000000004, 230.42949999999996, 0.503052, 18.63933
9.288720000000001, 0.12330999999999999, 1.370900000000006, 0.021232, 1.3668399999999998

20.892,0.267376,198.898,0.498318,16.7562
9.04191, 0.122593,83.2673,0.197334,6.92663
26.4805,0.321682,313.643,0.699551,25.4554
17.6139,0.237026,84.5844,0.217731,8.18291