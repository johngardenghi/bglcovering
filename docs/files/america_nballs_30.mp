prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(30);
u := 20cm;
path pols[];
pair p[];
R := .0566320690u;
D := 2*R;

pols[1] = 
(        0.2250000000u,        1.2000000000u)--
(        0.1750000000u,        1.1900000000u)--
(        0.1350000000u,        1.1500000000u)--
(        0.1375000000u,        1.1075000000u)--
(        0.1500000000u,        1.0750000000u)--
(        0.2000000000u,        1.1000000000u)--
(        0.2250000000u,        1.2000000000u)--
cycle;
pols[2] = 
(        0.2250000000u,        1.2000000000u)--
(        0.2000000000u,        1.1000000000u)--
(        0.2750000000u,        1.1000000000u)--
(        0.2900000000u,        1.1900000000u)--
(        0.2250000000u,        1.2000000000u)--
cycle;
pols[3] = 
(        0.3000000000u,        1.0500000000u)--
(        0.3200000000u,        1.0000000000u)--
(        0.5000000000u,        1.0000000000u)--
(        0.3000000000u,        1.0500000000u)--
cycle;
pols[4] = 
(        0.2750000000u,        1.1000000000u)--
(        0.3000000000u,        1.0500000000u)--
(        0.5000000000u,        1.0000000000u)--
(        0.5000000000u,        1.0750000000u)--
(        0.4500000000u,        1.1750000000u)--
(        0.3650000000u,        1.1850000000u)--
(        0.2900000000u,        1.1900000000u)--
(        0.2750000000u,        1.1000000000u)--
cycle;
pols[5] = 
(        0.5000000000u,        1.0000000000u)--
(        0.5500000000u,        0.9550000000u)--
(        0.5600000000u,        1.0100000000u)--
(        0.5000000000u,        1.0750000000u)--
(        0.5000000000u,        1.0000000000u)--
cycle;
pols[6] = 
(        0.5000000000u,        1.0750000000u)--
(        0.5000000000u,        1.1100000000u)--
(        0.4500000000u,        1.1750000000u)--
(        0.5000000000u,        1.0750000000u)--
cycle;
pols[7] = 
(        0.5000000000u,        1.1100000000u)--
(        0.5100000000u,        1.1650000000u)--
(        0.4500000000u,        1.1750000000u)--
(        0.5000000000u,        1.1100000000u)--
cycle;
pols[8] = 
(        0.5000000000u,        1.1100000000u)--
(        0.5750000000u,        1.1500000000u)--
(        0.5500000000u,        1.2300000000u)--
(        0.5100000000u,        1.1650000000u)--
(        0.5000000000u,        1.1100000000u)--
cycle;
pols[9] = 
(        0.5500000000u,        0.9550000000u)--
(        0.5700000000u,        0.9200000000u)--
(        0.6250000000u,        0.9750000000u)--
(        0.6200000000u,        0.9950000000u)--
(        0.5600000000u,        1.0100000000u)--
(        0.5500000000u,        0.9550000000u)--
cycle;
pols[10] = 
(        0.6200000000u,        0.9950000000u)--
(        0.6900000000u,        1.0300000000u)--
(        0.5900000000u,        1.1200000000u)--
(        0.5600000000u,        1.0250000000u)--
(        0.5600000000u,        1.0100000000u)--
(        0.6200000000u,        0.9950000000u)--
cycle;
pols[11] = 
(        0.6250000000u,        0.9750000000u)--
(        0.6550000000u,        0.9750000000u)--
(        0.6200000000u,        0.9950000000u)--
(        0.6250000000u,        0.9750000000u)--
cycle;
pols[12] = 
(        0.3200000000u,        1.0000000000u)--
(        0.3050000000u,        0.9750000000u)--
(        0.3000000000u,        0.9350000000u)--
(        0.3100000000u,        0.9100000000u)--
(        0.3300000000u,        0.8800000000u)--
(        0.3400000000u,        0.8750000000u)--
(        0.3450000000u,        0.8750000000u)--
(        0.5650000000u,        0.8900000000u)--
(        0.5700000000u,        0.9200000000u)--
(        0.5500000000u,        0.9550000000u)--
(        0.5000000000u,        1.0000000000u)--
(        0.3200000000u,        1.0000000000u)--
cycle;
pols[13] = 
(        0.3450000000u,        0.8750000000u)--
(        0.5350000000u,        0.8700000000u)--
(        0.5650000000u,        0.8900000000u)--
(        0.3450000000u,        0.8750000000u)--
cycle;
pols[14] = 
(        0.5200000000u,        0.8600000000u)--
(        0.5250000000u,        0.8300000000u)--
(        0.5300000000u,        0.8300000000u)--
(        0.5350000000u,        0.8700000000u)--
(        0.5200000000u,        0.8600000000u)--
cycle;
pols[15] = 
(        0.3450000000u,        0.8750000000u)--
(        0.4650000000u,        0.8600000000u)--
(        0.5200000000u,        0.8600000000u)--
(        0.5350000000u,        0.8700000000u)--
(        0.3450000000u,        0.8750000000u)--
cycle;
pols[16] = 
(        0.3450000000u,        0.8750000000u)--
(        0.4200000000u,        0.8300000000u)--
(        0.4650000000u,        0.8600000000u)--
(        0.3450000000u,        0.8750000000u)--
cycle;
pols[17] = 
(        0.3450000000u,        0.8750000000u)--
(        0.3700000000u,        0.8300000000u)--
(        0.3900000000u,        0.7950000000u)--
(        0.4250000000u,        0.8000000000u)--
(        0.4200000000u,        0.8300000000u)--
(        0.3450000000u,        0.8750000000u)--
cycle;
pols[18] = 
(        0.3900000000u,        0.7950000000u)--
(        0.3850000000u,        0.7900000000u)--
(        0.4250000000u,        0.7650000000u)--
(        0.4450000000u,        0.7650000000u)--
(        0.4500000000u,        0.7800000000u)--
(        0.4250000000u,        0.8000000000u)--
(        0.3900000000u,        0.7950000000u)--
cycle;
pols[19] = 
(        0.4450000000u,        0.7650000000u)--
(        0.4600000000u,        0.7500000000u)--
(        0.4700000000u,        0.7650000000u)--
(        0.4650000000u,        0.7750000000u)--
(        0.4500000000u,        0.7800000000u)--
(        0.4450000000u,        0.7650000000u)--
cycle;
pols[20] = 
(        0.4650000000u,        0.7750000000u)--
(        0.4850000000u,        0.7800000000u)--
(        0.4950000000u,        0.8000000000u)--
(        0.4750000000u,        0.8000000000u)--
(        0.4500000000u,        0.7800000000u)--
(        0.4650000000u,        0.7750000000u)--
cycle;
pols[21] = 
(        0.3300000000u,        0.8800000000u)--
(        0.3400000000u,        0.8400000000u)--
(        0.3500000000u,        0.8400000000u)--
(        0.3400000000u,        0.8750000000u)--
(        0.3300000000u,        0.8800000000u)--
cycle;
pols[22] = 
(        0.3400000000u,        0.8400000000u)--
(        0.3550000000u,        0.8150000000u)--
(        0.3600000000u,        0.8150000000u)--
(        0.3500000000u,        0.8400000000u)--
(        0.3400000000u,        0.8400000000u)--
cycle;
pols[23] = 
(        0.4600000000u,        0.7500000000u)--
(        0.4850000000u,        0.7350000000u)--
(        0.5100000000u,        0.7250000000u)--
(        0.5100000000u,        0.7650000000u)--
(        0.4700000000u,        0.7650000000u)--
(        0.4600000000u,        0.7500000000u)--
cycle;
pols[24] = 
(        0.4850000000u,        0.7350000000u)--
(        0.5000000000u,        0.7200000000u)--
(        0.5400000000u,        0.7050000000u)--
(        0.5450000000u,        0.7100000000u)--
(        0.5100000000u,        0.7250000000u)--
(        0.4850000000u,        0.7350000000u)--
cycle;
pols[25] = 
(        0.5200000000u,        0.8100000000u)--
(        0.5500000000u,        0.7900000000u)--
(        0.5650000000u,        0.8000000000u)--
(        0.5200000000u,        0.8150000000u)--
(        0.5200000000u,        0.8100000000u)--
cycle;
pols[26] = 
(        0.5350000000u,        0.6600000000u)--
(        0.5250000000u,        0.6250000000u)--
(        0.5350000000u,        0.5625000000u)--
(        0.5700000000u,        0.5300000000u)--
(        0.7100000000u,        0.4850000000u)--
(        0.7500000000u,        0.5000000000u)--
(        0.7650000000u,        0.5400000000u)--
(        0.7650000000u,        0.5650000000u)--
(        0.5350000000u,        0.6600000000u)--
cycle;
pols[27] = 
(        0.6100000000u,        0.2700000000u)--
(        0.5950000000u,        0.2650000000u)--
(        0.6100000000u,        0.2600000000u)--
(        0.6250000000u,        0.2650000000u)--
(        0.6100000000u,        0.2700000000u)--
cycle;
pols[28] = 
(        0.7650000000u,        0.5650000000u)--
(        0.7850000000u,        0.6100000000u)--
(        0.7300000000u,        0.6400000000u)--
(        0.5450000000u,        0.7100000000u)--
(        0.5400000000u,        0.7050000000u)--
(        0.5350000000u,        0.6600000000u)--
(        0.7650000000u,        0.5650000000u)--
cycle;
pols[29] = 
(        0.7300000000u,        0.6400000000u)--
(        0.6900000000u,        0.6750000000u)--
(        0.6450000000u,        0.7050000000u)--
(        0.6050000000u,        0.7250000000u)--
(        0.5800000000u,        0.7300000000u)--
(        0.5450000000u,        0.7100000000u)--
(        0.7300000000u,        0.6400000000u)--
cycle;
pols[30] = 
(        0.6450000000u,        0.7050000000u)--
(        0.6250000000u,        0.7250000000u)--
(        0.6050000000u,        0.7250000000u)--
(        0.6450000000u,        0.7050000000u)--
cycle;
pols[31] = 
(        0.5700000000u,        0.5300000000u)--
(        0.5700000000u,        0.3750000000u)--
(        0.5750000000u,        0.2850000000u)--
(        0.5900000000u,        0.2750000000u)--
(        0.6050000000u,        0.2800000000u)--
(        0.6150000000u,        0.3350000000u)--
(        0.5700000000u,        0.5300000000u)--
cycle;
pols[32] = 
(        0.6150000000u,        0.3350000000u)--
(        0.6300000000u,        0.3850000000u)--
(        0.5700000000u,        0.5300000000u)--
(        0.6150000000u,        0.3350000000u)--
cycle;
pols[33] = 
(        0.6300000000u,        0.3850000000u)--
(        0.6600000000u,        0.3850000000u)--
(        0.6550000000u,        0.4200000000u)--
(        0.5700000000u,        0.5300000000u)--
(        0.6300000000u,        0.3850000000u)--
cycle;
pols[34] = 
(        0.6550000000u,        0.4200000000u)--
(        0.6750000000u,        0.4150000000u)--
(        0.6850000000u,        0.4300000000u)--
(        0.7100000000u,        0.4850000000u)--
(        0.5700000000u,        0.5300000000u)--
(        0.6550000000u,        0.4200000000u)--
cycle;
fill pols[1] withcolor 0.9Dandelion;
fill pols[2] withcolor 0.9Dandelion;
fill pols[3] withcolor 0.9Dandelion;
fill pols[4] withcolor 0.9Dandelion;
fill pols[5] withcolor 0.9Dandelion;
fill pols[6] withcolor 0.9Dandelion;
fill pols[7] withcolor 0.9Dandelion;
fill pols[8] withcolor 0.9Dandelion;
fill pols[9] withcolor 0.9Dandelion;
fill pols[10] withcolor 0.9Dandelion;
fill pols[11] withcolor 0.9Dandelion;
fill pols[12] withcolor 0.9Dandelion;
fill pols[13] withcolor 0.9Dandelion;
fill pols[14] withcolor 0.9Dandelion;
fill pols[15] withcolor 0.9Dandelion;
fill pols[16] withcolor 0.9Dandelion;
fill pols[17] withcolor 0.9Dandelion;
fill pols[18] withcolor 0.9Dandelion;
fill pols[19] withcolor 0.9Dandelion;
fill pols[20] withcolor 0.9Dandelion;
fill pols[21] withcolor 0.9Dandelion;
fill pols[22] withcolor 0.9Dandelion;
fill pols[23] withcolor 0.9Dandelion;
fill pols[24] withcolor 0.9Dandelion;
fill pols[25] withcolor 0.9Dandelion;
fill pols[26] withcolor 0.9Dandelion;
fill pols[27] withcolor 0.9Dandelion;
fill pols[28] withcolor 0.9Dandelion;
fill pols[29] withcolor 0.9Dandelion;
fill pols[30] withcolor 0.9Dandelion;
fill pols[31] withcolor 0.9Dandelion;
fill pols[32] withcolor 0.9Dandelion;
fill pols[33] withcolor 0.9Dandelion;
fill pols[34] withcolor 0.9Dandelion;
draw (        0.2000000000u,        1.1000000000u)--(        0.2250000000u,        1.2000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2250000000u,        1.2000000000u)--(        0.2000000000u,        1.1000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2750000000u,        1.1000000000u)--(        0.2900000000u,        1.1900000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3200000000u,        1.0000000000u)--(        0.5000000000u,        1.0000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5000000000u,        1.0000000000u)--(        0.3000000000u,        1.0500000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3000000000u,        1.0500000000u)--(        0.5000000000u,        1.0000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5000000000u,        1.0000000000u)--(        0.5000000000u,        1.0750000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5000000000u,        1.0750000000u)--(        0.4500000000u,        1.1750000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2900000000u,        1.1900000000u)--(        0.2750000000u,        1.1000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5000000000u,        1.0000000000u)--(        0.5500000000u,        0.9550000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5500000000u,        0.9550000000u)--(        0.5600000000u,        1.0100000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5000000000u,        1.0750000000u)--(        0.5000000000u,        1.0000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5000000000u,        1.1100000000u)--(        0.4500000000u,        1.1750000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4500000000u,        1.1750000000u)--(        0.5000000000u,        1.0750000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5000000000u,        1.1100000000u)--(        0.5100000000u,        1.1650000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4500000000u,        1.1750000000u)--(        0.5000000000u,        1.1100000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5100000000u,        1.1650000000u)--(        0.5000000000u,        1.1100000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5500000000u,        0.9550000000u)--(        0.5700000000u,        0.9200000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.6250000000u,        0.9750000000u)--(        0.6200000000u,        0.9950000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.6200000000u,        0.9950000000u)--(        0.5600000000u,        1.0100000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5600000000u,        1.0100000000u)--(        0.5500000000u,        0.9550000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5600000000u,        1.0100000000u)--(        0.6200000000u,        0.9950000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.6200000000u,        0.9950000000u)--(        0.6250000000u,        0.9750000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3300000000u,        0.8800000000u)--(        0.3400000000u,        0.8750000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3450000000u,        0.8750000000u)--(        0.5650000000u,        0.8900000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5700000000u,        0.9200000000u)--(        0.5500000000u,        0.9550000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5500000000u,        0.9550000000u)--(        0.5000000000u,        1.0000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5000000000u,        1.0000000000u)--(        0.3200000000u,        1.0000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3450000000u,        0.8750000000u)--(        0.5350000000u,        0.8700000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5650000000u,        0.8900000000u)--(        0.3450000000u,        0.8750000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5350000000u,        0.8700000000u)--(        0.5200000000u,        0.8600000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3450000000u,        0.8750000000u)--(        0.4650000000u,        0.8600000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5200000000u,        0.8600000000u)--(        0.5350000000u,        0.8700000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5350000000u,        0.8700000000u)--(        0.3450000000u,        0.8750000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3450000000u,        0.8750000000u)--(        0.4200000000u,        0.8300000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4650000000u,        0.8600000000u)--(        0.3450000000u,        0.8750000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3900000000u,        0.7950000000u)--(        0.4250000000u,        0.8000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4200000000u,        0.8300000000u)--(        0.3450000000u,        0.8750000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4450000000u,        0.7650000000u)--(        0.4500000000u,        0.7800000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4250000000u,        0.8000000000u)--(        0.3900000000u,        0.7950000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4600000000u,        0.7500000000u)--(        0.4700000000u,        0.7650000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4650000000u,        0.7750000000u)--(        0.4500000000u,        0.7800000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4500000000u,        0.7800000000u)--(        0.4450000000u,        0.7650000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4500000000u,        0.7800000000u)--(        0.4650000000u,        0.7750000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3400000000u,        0.8400000000u)--(        0.3500000000u,        0.8400000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3400000000u,        0.8750000000u)--(        0.3300000000u,        0.8800000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3500000000u,        0.8400000000u)--(        0.3400000000u,        0.8400000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4850000000u,        0.7350000000u)--(        0.5100000000u,        0.7250000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4700000000u,        0.7650000000u)--(        0.4600000000u,        0.7500000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5400000000u,        0.7050000000u)--(        0.5450000000u,        0.7100000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5100000000u,        0.7250000000u)--(        0.4850000000u,        0.7350000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5700000000u,        0.5300000000u)--(        0.7100000000u,        0.4850000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.7650000000u,        0.5650000000u)--(        0.5350000000u,        0.6600000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.7300000000u,        0.6400000000u)--(        0.5450000000u,        0.7100000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5450000000u,        0.7100000000u)--(        0.5400000000u,        0.7050000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5350000000u,        0.6600000000u)--(        0.7650000000u,        0.5650000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.6450000000u,        0.7050000000u)--(        0.6050000000u,        0.7250000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5450000000u,        0.7100000000u)--(        0.7300000000u,        0.6400000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.6050000000u,        0.7250000000u)--(        0.6450000000u,        0.7050000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.6150000000u,        0.3350000000u)--(        0.5700000000u,        0.5300000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.6300000000u,        0.3850000000u)--(        0.5700000000u,        0.5300000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5700000000u,        0.5300000000u)--(        0.6150000000u,        0.3350000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.6550000000u,        0.4200000000u)--(        0.5700000000u,        0.5300000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5700000000u,        0.5300000000u)--(        0.6300000000u,        0.3850000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.7100000000u,        0.4850000000u)--(        0.5700000000u,        0.5300000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5700000000u,        0.5300000000u)--(        0.6550000000u,        0.4200000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2250000000u,        1.2000000000u)--(        0.1750000000u,        1.1900000000u) withpen pencircle scaled 1.0;
draw (        0.1750000000u,        1.1900000000u)--(        0.1350000000u,        1.1500000000u) withpen pencircle scaled 1.0;
draw (        0.1350000000u,        1.1500000000u)--(        0.1375000000u,        1.1075000000u) withpen pencircle scaled 1.0;
draw (        0.1375000000u,        1.1075000000u)--(        0.1500000000u,        1.0750000000u) withpen pencircle scaled 1.0;
draw (        0.1500000000u,        1.0750000000u)--(        0.2000000000u,        1.1000000000u) withpen pencircle scaled 1.0;
draw (        0.2000000000u,        1.1000000000u)--(        0.2750000000u,        1.1000000000u) withpen pencircle scaled 1.0;
draw (        0.2900000000u,        1.1900000000u)--(        0.2250000000u,        1.2000000000u) withpen pencircle scaled 1.0;
draw (        0.3000000000u,        1.0500000000u)--(        0.3200000000u,        1.0000000000u) withpen pencircle scaled 1.0;
draw (        0.2750000000u,        1.1000000000u)--(        0.3000000000u,        1.0500000000u) withpen pencircle scaled 1.0;
draw (        0.4500000000u,        1.1750000000u)--(        0.3650000000u,        1.1850000000u) withpen pencircle scaled 1.0;
draw (        0.3650000000u,        1.1850000000u)--(        0.2900000000u,        1.1900000000u) withpen pencircle scaled 1.0;
draw (        0.5600000000u,        1.0100000000u)--(        0.5000000000u,        1.0750000000u) withpen pencircle scaled 1.0;
draw (        0.5000000000u,        1.0750000000u)--(        0.5000000000u,        1.1100000000u) withpen pencircle scaled 1.0;
draw (        0.5100000000u,        1.1650000000u)--(        0.4500000000u,        1.1750000000u) withpen pencircle scaled 1.0;
draw (        0.5000000000u,        1.1100000000u)--(        0.5750000000u,        1.1500000000u) withpen pencircle scaled 1.0;
draw (        0.5750000000u,        1.1500000000u)--(        0.5500000000u,        1.2300000000u) withpen pencircle scaled 1.0;
draw (        0.5500000000u,        1.2300000000u)--(        0.5100000000u,        1.1650000000u) withpen pencircle scaled 1.0;
draw (        0.5700000000u,        0.9200000000u)--(        0.6250000000u,        0.9750000000u) withpen pencircle scaled 1.0;
draw (        0.6200000000u,        0.9950000000u)--(        0.6900000000u,        1.0300000000u) withpen pencircle scaled 1.0;
draw (        0.6900000000u,        1.0300000000u)--(        0.5900000000u,        1.1200000000u) withpen pencircle scaled 1.0;
draw (        0.5900000000u,        1.1200000000u)--(        0.5600000000u,        1.0250000000u) withpen pencircle scaled 1.0;
draw (        0.5600000000u,        1.0250000000u)--(        0.5600000000u,        1.0100000000u) withpen pencircle scaled 1.0;
draw (        0.6250000000u,        0.9750000000u)--(        0.6550000000u,        0.9750000000u) withpen pencircle scaled 1.0;
draw (        0.6550000000u,        0.9750000000u)--(        0.6200000000u,        0.9950000000u) withpen pencircle scaled 1.0;
draw (        0.3200000000u,        1.0000000000u)--(        0.3050000000u,        0.9750000000u) withpen pencircle scaled 1.0;
draw (        0.3050000000u,        0.9750000000u)--(        0.3000000000u,        0.9350000000u) withpen pencircle scaled 1.0;
draw (        0.3000000000u,        0.9350000000u)--(        0.3100000000u,        0.9100000000u) withpen pencircle scaled 1.0;
draw (        0.3100000000u,        0.9100000000u)--(        0.3300000000u,        0.8800000000u) withpen pencircle scaled 1.0;
draw (        0.3400000000u,        0.8750000000u)--(        0.3450000000u,        0.8750000000u) withpen pencircle scaled 1.0;
draw (        0.5650000000u,        0.8900000000u)--(        0.5700000000u,        0.9200000000u) withpen pencircle scaled 1.0;
draw (        0.5350000000u,        0.8700000000u)--(        0.5650000000u,        0.8900000000u) withpen pencircle scaled 1.0;
draw (        0.5200000000u,        0.8600000000u)--(        0.5250000000u,        0.8300000000u) withpen pencircle scaled 1.0;
draw (        0.5250000000u,        0.8300000000u)--(        0.5300000000u,        0.8300000000u) withpen pencircle scaled 1.0;
draw (        0.5300000000u,        0.8300000000u)--(        0.5350000000u,        0.8700000000u) withpen pencircle scaled 1.0;
draw (        0.4650000000u,        0.8600000000u)--(        0.5200000000u,        0.8600000000u) withpen pencircle scaled 1.0;
draw (        0.4200000000u,        0.8300000000u)--(        0.4650000000u,        0.8600000000u) withpen pencircle scaled 1.0;
draw (        0.3450000000u,        0.8750000000u)--(        0.3700000000u,        0.8300000000u) withpen pencircle scaled 1.0;
draw (        0.3700000000u,        0.8300000000u)--(        0.3900000000u,        0.7950000000u) withpen pencircle scaled 1.0;
draw (        0.4250000000u,        0.8000000000u)--(        0.4200000000u,        0.8300000000u) withpen pencircle scaled 1.0;
draw (        0.3900000000u,        0.7950000000u)--(        0.3850000000u,        0.7900000000u) withpen pencircle scaled 1.0;
draw (        0.3850000000u,        0.7900000000u)--(        0.4250000000u,        0.7650000000u) withpen pencircle scaled 1.0;
draw (        0.4250000000u,        0.7650000000u)--(        0.4450000000u,        0.7650000000u) withpen pencircle scaled 1.0;
draw (        0.4500000000u,        0.7800000000u)--(        0.4250000000u,        0.8000000000u) withpen pencircle scaled 1.0;
draw (        0.4450000000u,        0.7650000000u)--(        0.4600000000u,        0.7500000000u) withpen pencircle scaled 1.0;
draw (        0.4700000000u,        0.7650000000u)--(        0.4650000000u,        0.7750000000u) withpen pencircle scaled 1.0;
draw (        0.4650000000u,        0.7750000000u)--(        0.4850000000u,        0.7800000000u) withpen pencircle scaled 1.0;
draw (        0.4850000000u,        0.7800000000u)--(        0.4950000000u,        0.8000000000u) withpen pencircle scaled 1.0;
draw (        0.4950000000u,        0.8000000000u)--(        0.4750000000u,        0.8000000000u) withpen pencircle scaled 1.0;
draw (        0.4750000000u,        0.8000000000u)--(        0.4500000000u,        0.7800000000u) withpen pencircle scaled 1.0;
draw (        0.3300000000u,        0.8800000000u)--(        0.3400000000u,        0.8400000000u) withpen pencircle scaled 1.0;
draw (        0.3500000000u,        0.8400000000u)--(        0.3400000000u,        0.8750000000u) withpen pencircle scaled 1.0;
draw (        0.3400000000u,        0.8400000000u)--(        0.3550000000u,        0.8150000000u) withpen pencircle scaled 1.0;
draw (        0.3550000000u,        0.8150000000u)--(        0.3600000000u,        0.8150000000u) withpen pencircle scaled 1.0;
draw (        0.3600000000u,        0.8150000000u)--(        0.3500000000u,        0.8400000000u) withpen pencircle scaled 1.0;
draw (        0.4600000000u,        0.7500000000u)--(        0.4850000000u,        0.7350000000u) withpen pencircle scaled 1.0;
draw (        0.5100000000u,        0.7250000000u)--(        0.5100000000u,        0.7650000000u) withpen pencircle scaled 1.0;
draw (        0.5100000000u,        0.7650000000u)--(        0.4700000000u,        0.7650000000u) withpen pencircle scaled 1.0;
draw (        0.4850000000u,        0.7350000000u)--(        0.5000000000u,        0.7200000000u) withpen pencircle scaled 1.0;
draw (        0.5000000000u,        0.7200000000u)--(        0.5400000000u,        0.7050000000u) withpen pencircle scaled 1.0;
draw (        0.5450000000u,        0.7100000000u)--(        0.5100000000u,        0.7250000000u) withpen pencircle scaled 1.0;
draw (        0.5200000000u,        0.8100000000u)--(        0.5500000000u,        0.7900000000u) withpen pencircle scaled 1.0;
draw (        0.5500000000u,        0.7900000000u)--(        0.5650000000u,        0.8000000000u) withpen pencircle scaled 1.0;
draw (        0.5650000000u,        0.8000000000u)--(        0.5200000000u,        0.8150000000u) withpen pencircle scaled 1.0;
draw (        0.5200000000u,        0.8150000000u)--(        0.5200000000u,        0.8100000000u) withpen pencircle scaled 1.0;
draw (        0.5350000000u,        0.6600000000u)--(        0.5250000000u,        0.6250000000u) withpen pencircle scaled 1.0;
draw (        0.5250000000u,        0.6250000000u)--(        0.5350000000u,        0.5625000000u) withpen pencircle scaled 1.0;
draw (        0.5350000000u,        0.5625000000u)--(        0.5700000000u,        0.5300000000u) withpen pencircle scaled 1.0;
draw (        0.7100000000u,        0.4850000000u)--(        0.7500000000u,        0.5000000000u) withpen pencircle scaled 1.0;
draw (        0.7500000000u,        0.5000000000u)--(        0.7650000000u,        0.5400000000u) withpen pencircle scaled 1.0;
draw (        0.7650000000u,        0.5400000000u)--(        0.7650000000u,        0.5650000000u) withpen pencircle scaled 1.0;
draw (        0.6100000000u,        0.2700000000u)--(        0.5950000000u,        0.2650000000u) withpen pencircle scaled 1.0;
draw (        0.5950000000u,        0.2650000000u)--(        0.6100000000u,        0.2600000000u) withpen pencircle scaled 1.0;
draw (        0.6100000000u,        0.2600000000u)--(        0.6250000000u,        0.2650000000u) withpen pencircle scaled 1.0;
draw (        0.6250000000u,        0.2650000000u)--(        0.6100000000u,        0.2700000000u) withpen pencircle scaled 1.0;
draw (        0.7650000000u,        0.5650000000u)--(        0.7850000000u,        0.6100000000u) withpen pencircle scaled 1.0;
draw (        0.7850000000u,        0.6100000000u)--(        0.7300000000u,        0.6400000000u) withpen pencircle scaled 1.0;
draw (        0.5400000000u,        0.7050000000u)--(        0.5350000000u,        0.6600000000u) withpen pencircle scaled 1.0;
draw (        0.7300000000u,        0.6400000000u)--(        0.6900000000u,        0.6750000000u) withpen pencircle scaled 1.0;
draw (        0.6900000000u,        0.6750000000u)--(        0.6450000000u,        0.7050000000u) withpen pencircle scaled 1.0;
draw (        0.6050000000u,        0.7250000000u)--(        0.5800000000u,        0.7300000000u) withpen pencircle scaled 1.0;
draw (        0.5800000000u,        0.7300000000u)--(        0.5450000000u,        0.7100000000u) withpen pencircle scaled 1.0;
draw (        0.6450000000u,        0.7050000000u)--(        0.6250000000u,        0.7250000000u) withpen pencircle scaled 1.0;
draw (        0.6250000000u,        0.7250000000u)--(        0.6050000000u,        0.7250000000u) withpen pencircle scaled 1.0;
draw (        0.5700000000u,        0.5300000000u)--(        0.5700000000u,        0.3750000000u) withpen pencircle scaled 1.0;
draw (        0.5700000000u,        0.3750000000u)--(        0.5750000000u,        0.2850000000u) withpen pencircle scaled 1.0;
draw (        0.5750000000u,        0.2850000000u)--(        0.5900000000u,        0.2750000000u) withpen pencircle scaled 1.0;
draw (        0.5900000000u,        0.2750000000u)--(        0.6050000000u,        0.2800000000u) withpen pencircle scaled 1.0;
draw (        0.6050000000u,        0.2800000000u)--(        0.6150000000u,        0.3350000000u) withpen pencircle scaled 1.0;
draw (        0.6150000000u,        0.3350000000u)--(        0.6300000000u,        0.3850000000u) withpen pencircle scaled 1.0;
draw (        0.6300000000u,        0.3850000000u)--(        0.6600000000u,        0.3850000000u) withpen pencircle scaled 1.0;
draw (        0.6600000000u,        0.3850000000u)--(        0.6550000000u,        0.4200000000u) withpen pencircle scaled 1.0;
draw (        0.6550000000u,        0.4200000000u)--(        0.6750000000u,        0.4150000000u) withpen pencircle scaled 1.0;
draw (        0.6750000000u,        0.4150000000u)--(        0.6850000000u,        0.4300000000u) withpen pencircle scaled 1.0;
draw (        0.6850000000u,        0.4300000000u)--(        0.7100000000u,        0.4850000000u) withpen pencircle scaled 1.0;

x1 = .5977959712u;
y1 = .3146704404u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .6078216380u;
y2 = .4070168758u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .6603715785u;
y3 = .4697108661u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .7241203526u;
y4 = .5424471453u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .5730310481u;
y5 = .5057205288u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .6383198021u;
y6 = .5723913625u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .7285593817u;
y7 = .6146534470u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .5564860967u;
y8 = .6094756513u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .6453383737u;
y9 = .6678577438u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = .5750341528u;
y10 = .7022389274u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;

x11 = .5986442999u;
y11 = .9690559255u;
draw fullcircle scaled D shifted z11 withcolor RoyalBlue withpen pencircle scaled 1.0;

x12 = .6413140162u;
y12 = 1.0589366935u;
draw fullcircle scaled D shifted z12 withcolor RoyalBlue withpen pencircle scaled 1.0;

x13 = .4658867848u;
y13 = 1.1369830188u;
draw fullcircle scaled D shifted z13 withcolor RoyalBlue withpen pencircle scaled 1.0;

x14 = .5342012240u;
y14 = .8858843934u;
draw fullcircle scaled D shifted z14 withcolor RoyalBlue withpen pencircle scaled 1.0;

x15 = .5159478190u;
y15 = .9857030072u;
draw fullcircle scaled D shifted z15 withcolor RoyalBlue withpen pencircle scaled 1.0;

x16 = .5616404133u;
y16 = 1.0709773000u;
draw fullcircle scaled D shifted z16 withcolor RoyalBlue withpen pencircle scaled 1.0;

x17 = .5173512273u;
y17 = .7693300003u;
draw fullcircle scaled D shifted z17 withcolor RoyalBlue withpen pencircle scaled 1.0;

x18 = .4540169597u;
y18 = .8986622596u;
draw fullcircle scaled D shifted z18 withcolor RoyalBlue withpen pencircle scaled 1.0;

x19 = .4153400960u;
y19 = .9691972157u;
draw fullcircle scaled D shifted z19 withcolor RoyalBlue withpen pencircle scaled 1.0;

x20 = .4574861055u;
y20 = 1.0543604725u;
draw fullcircle scaled D shifted z20 withcolor RoyalBlue withpen pencircle scaled 1.0;

x21 = .4331603000u;
y21 = .7933476654u;
draw fullcircle scaled D shifted z21 withcolor RoyalBlue withpen pencircle scaled 1.0;

x22 = .3726317461u;
y22 = .8688178962u;
draw fullcircle scaled D shifted z22 withcolor RoyalBlue withpen pencircle scaled 1.0;

x23 = .3562482186u;
y23 = 1.0466060922u;
draw fullcircle scaled D shifted z23 withcolor RoyalBlue withpen pencircle scaled 1.0;

x24 = .5350701699u;
y24 = 1.1753705256u;
draw fullcircle scaled D shifted z24 withcolor RoyalBlue withpen pencircle scaled 1.0;

x25 = .3287131465u;
y25 = .9490902755u;
draw fullcircle scaled D shifted z25 withcolor RoyalBlue withpen pencircle scaled 1.0;

x26 = .3943743593u;
y26 = 1.1311694234u;
draw fullcircle scaled D shifted z26 withcolor RoyalBlue withpen pencircle scaled 1.0;

x27 = .2766607852u;
y27 = 1.1022111932u;
draw fullcircle scaled D shifted z27 withcolor RoyalBlue withpen pencircle scaled 1.0;

x28 = .3255095011u;
y28 = 1.1544365501u;
draw fullcircle scaled D shifted z28 withcolor RoyalBlue withpen pencircle scaled 1.0;

x29 = .1706887017u;
y29 = 1.1277268819u;
draw fullcircle scaled D shifted z29 withcolor RoyalBlue withpen pencircle scaled 1.0;

x30 = .2259489389u;
y30 = 1.1865507976u;
draw fullcircle scaled D shifted z30 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.5977959712u,        0.3146704404u)--(        0.6250000000u,        0.2650000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5977959712u,        0.3146704404u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5977959712u,        0.3146704404u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6250000000u,        0.2650000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6250000000u,        0.2650000000u) withcolor Yellow;
draw (        0.5977959712u,        0.3146704404u)--(        0.5705929097u,        0.3643412032u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6078216380u,        0.4070168758u)--(        0.5705929097u,        0.3643412032u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5977959712u,        0.3146704404u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5977959712u,        0.3146704404u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6078216380u,        0.4070168758u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6078216380u,        0.4070168758u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5705929097u,        0.3643412032u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5705929097u,        0.3643412032u) withcolor Green;
draw (        0.6078216380u,        0.4070168758u)--(        0.6600000000u,        0.3850000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6078216380u,        0.4070168758u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6078216380u,        0.4070168758u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6600000000u,        0.3850000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6600000000u,        0.3850000000u) withcolor Yellow;
draw (        0.6078216380u,        0.4070168758u)--(        0.5700017450u,        0.4491695380u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5730310481u,        0.5057205288u)--(        0.5700017450u,        0.4491695380u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6078216380u,        0.4070168758u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6078216380u,        0.4070168758u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5730310481u,        0.5057205288u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5730310481u,        0.5057205288u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5700017450u,        0.4491695380u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5700017450u,        0.4491695380u) withcolor Green;
draw (        0.6603715785u,        0.4697108661u)--(        0.6750000000u,        0.4150000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6603715785u,        0.4697108661u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6603715785u,        0.4697108661u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6750000000u,        0.4150000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6750000000u,        0.4150000000u) withcolor Yellow;
draw (        0.6603715785u,        0.4697108661u)--(        0.6700821452u,        0.5255042003u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7241203526u,        0.5424471453u)--(        0.6700821452u,        0.5255042003u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6603715785u,        0.4697108661u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6603715785u,        0.4697108661u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7241203526u,        0.5424471453u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7241203526u,        0.5424471453u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6700821452u,        0.5255042003u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6700821452u,        0.5255042003u) withcolor MidnightBlue;
draw (        0.6603715785u,        0.4697108661u)--(        0.7144097859u,        0.4866538112u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7241203526u,        0.5424471453u)--(        0.7144097859u,        0.4866538112u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6603715785u,        0.4697108661u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6603715785u,        0.4697108661u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7241203526u,        0.5424471453u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7241203526u,        0.5424471453u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7144097859u,        0.4866538112u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7144097859u,        0.4866538112u) withcolor Green;
draw (        0.6603715785u,        0.4697108661u)--(        0.6286088425u,        0.5165971097u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5730310481u,        0.5057205288u)--(        0.6286088425u,        0.5165971097u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6603715785u,        0.4697108661u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6603715785u,        0.4697108661u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5730310481u,        0.5057205288u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5730310481u,        0.5057205288u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6286088425u,        0.5165971097u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6286088425u,        0.5165971097u) withcolor MidnightBlue;
draw (        0.6603715785u,        0.4697108661u)--(        0.6286100016u,        0.5165978949u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6383198021u,        0.5723913625u)--(        0.6286100016u,        0.5165978949u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6603715785u,        0.4697108661u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6603715785u,        0.4697108661u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6383198021u,        0.5723913625u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6383198021u,        0.5723913625u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6286100016u,        0.5165978949u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6286100016u,        0.5165978949u) withcolor MidnightBlue;
draw (        0.6603715785u,        0.4697108661u)--(        0.6700813790u,        0.5255043336u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6383198021u,        0.5723913625u)--(        0.6700813790u,        0.5255043336u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6603715785u,        0.4697108661u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6603715785u,        0.4697108661u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6383198021u,        0.5723913625u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6383198021u,        0.5723913625u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6700813790u,        0.5255043336u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6700813790u,        0.5255043336u) withcolor MidnightBlue;
draw (        0.7241203526u,        0.5424471453u)--(        0.6700819765u,        0.5255047383u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6383198021u,        0.5723913625u)--(        0.6700819765u,        0.5255047383u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7241203526u,        0.5424471453u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7241203526u,        0.5424471453u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6383198021u,        0.5723913625u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6383198021u,        0.5723913625u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6700819765u,        0.5255047383u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6700819765u,        0.5255047383u) withcolor MidnightBlue;
draw (        0.7241203526u,        0.5424471453u)--(        0.7698332924u,        0.5758764495u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7285593817u,        0.6146534470u)--(        0.7698332924u,        0.5758764495u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7241203526u,        0.5424471453u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7241203526u,        0.5424471453u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7285593817u,        0.6146534470u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7285593817u,        0.6146534470u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7698332924u,        0.5758764495u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7698332924u,        0.5758764495u) withcolor Green;
draw (        0.5730310481u,        0.5057205288u)--(        0.6286086426u,        0.5165981315u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6383198021u,        0.5723913625u)--(        0.6286086426u,        0.5165981315u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5730310481u,        0.5057205288u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5730310481u,        0.5057205288u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6383198021u,        0.5723913625u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6383198021u,        0.5723913625u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6286086426u,        0.5165981315u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6286086426u,        0.5165981315u) withcolor MidnightBlue;
draw (        0.5730310481u,        0.5057205288u)--(        0.5438688675u,        0.5542669856u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5564860967u,        0.6094756513u)--(        0.5438688675u,        0.5542669856u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5730310481u,        0.5057205288u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5730310481u,        0.5057205288u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5564860967u,        0.6094756513u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5564860967u,        0.6094756513u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5438688675u,        0.5542669856u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5438688675u,        0.5542669856u) withcolor Green;
draw (        0.6383198021u,        0.5723913625u)--(        0.6720207023u,        0.6179044438u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7285593817u,        0.6146534470u)--(        0.6720207023u,        0.6179044438u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6383198021u,        0.5723913625u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6383198021u,        0.5723913625u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7285593817u,        0.6146534470u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7285593817u,        0.6146534470u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6720207023u,        0.6179044438u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6720207023u,        0.6179044438u) withcolor MidnightBlue;
draw (        0.6383198021u,        0.5723913625u)--(        0.6116368931u,        0.6223434756u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5564860967u,        0.6094756513u)--(        0.6116368931u,        0.6223434756u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6383198021u,        0.5723913625u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6383198021u,        0.5723913625u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5564860967u,        0.6094756513u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5564860967u,        0.6094756513u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6116368931u,        0.6223434756u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6116368931u,        0.6223434756u) withcolor MidnightBlue;
draw (        0.6383198021u,        0.5723913625u)--(        0.6116381614u,        0.6223441531u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6453383737u,        0.6678577438u)--(        0.6116381614u,        0.6223441531u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6383198021u,        0.5723913625u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6383198021u,        0.5723913625u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6453383737u,        0.6678577438u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6453383737u,        0.6678577438u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6116381614u,        0.6223441531u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6116381614u,        0.6223441531u) withcolor MidnightBlue;
draw (        0.6383198021u,        0.5723913625u)--(        0.6720200144u,        0.6179049532u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6453383737u,        0.6678577438u)--(        0.6720200144u,        0.6179049532u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6383198021u,        0.5723913625u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6383198021u,        0.5723913625u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6453383737u,        0.6678577438u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6453383737u,        0.6678577438u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6720200144u,        0.6179049532u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6720200144u,        0.6179049532u) withcolor MidnightBlue;
draw (        0.7285593817u,        0.6146534470u)--(        0.7850000000u,        0.6100000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7285593817u,        0.6146534470u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7285593817u,        0.6146534470u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7850000000u,        0.6100000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7850000000u,        0.6100000000u) withcolor Yellow;
draw (        0.7285593817u,        0.6146534470u)--(        0.6720207543u,        0.6179053484u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6453383737u,        0.6678577438u)--(        0.6720207543u,        0.6179053484u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7285593817u,        0.6146534470u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7285593817u,        0.6146534470u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6453383737u,        0.6678577438u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6453383737u,        0.6678577438u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6720207543u,        0.6179053484u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6720207543u,        0.6179053484u) withcolor MidnightBlue;
draw (        0.7285593817u,        0.6146534470u)--(        0.7018770010u,        0.6646058424u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6453383737u,        0.6678577438u)--(        0.7018770010u,        0.6646058424u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7285593817u,        0.6146534470u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7285593817u,        0.6146534470u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6453383737u,        0.6678577438u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6453383737u,        0.6678577438u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7018770010u,        0.6646058424u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7018770010u,        0.6646058424u) withcolor Green;
draw (        0.5564860967u,        0.6094756513u)--(        0.6116364370u,        0.6223454299u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6453383737u,        0.6678577438u)--(        0.6116364370u,        0.6223454299u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5564860967u,        0.6094756513u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5564860967u,        0.6094756513u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6453383737u,        0.6678577438u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6453383737u,        0.6678577438u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6116364370u,        0.6223454299u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6116364370u,        0.6223454299u) withcolor MidnightBlue;
draw (        0.5564860967u,        0.6094756513u)--(        0.5352209665u,        0.6619636069u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5750341528u,        0.7022389274u)--(        0.5352209665u,        0.6619636069u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5564860967u,        0.6094756513u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5564860967u,        0.6094756513u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5750341528u,        0.7022389274u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5750341528u,        0.7022389274u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5352209665u,        0.6619636069u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5352209665u,        0.6619636069u) withcolor Green;
draw (        0.6453383737u,        0.6678577438u)--(        0.6281714182u,        0.7218252051u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5750341528u,        0.7022389274u)--(        0.6281714182u,        0.7218252051u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6453383737u,        0.6678577438u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6453383737u,        0.6678577438u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5750341528u,        0.7022389274u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5750341528u,        0.7022389274u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6281714182u,        0.7218252051u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6281714182u,        0.7218252051u) withcolor Green;
draw (        0.5750341528u,        0.7022389274u)--(        0.5193831293u,        0.7127343943u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5173512273u,        0.7693300003u)--(        0.5193831293u,        0.7127343943u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5750341528u,        0.7022389274u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5750341528u,        0.7022389274u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5173512273u,        0.7693300003u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5173512273u,        0.7693300003u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5193831293u,        0.7127343943u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5193831293u,        0.7127343943u) withcolor Green;
draw (        0.5986442999u,        0.9690559255u)--(        0.6550000000u,        0.9750000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5986442999u,        0.9690559255u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5986442999u,        0.9690559255u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6550000000u,        0.9750000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6550000000u,        0.9750000000u) withcolor Yellow;
draw (        0.5986442999u,        0.9690559255u)--(        0.5955317507u,        1.0256023957u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6413140162u,        1.0589366935u)--(        0.5955317507u,        1.0256023957u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5986442999u,        0.9690559255u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5986442999u,        0.9690559255u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6413140162u,        1.0589366935u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6413140162u,        1.0589366935u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5955317507u,        1.0256023957u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5955317507u,        1.0256023957u) withcolor MidnightBlue;
draw (        0.5986442999u,        0.9690559255u)--(        0.5498492102u,        0.9403117006u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5342012240u,        0.8858843934u)--(        0.5498492102u,        0.9403117006u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5986442999u,        0.9690559255u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5986442999u,        0.9690559255u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5342012240u,        0.8858843934u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5342012240u,        0.8858843934u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5498492102u,        0.9403117006u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5498492102u,        0.9403117006u) withcolor MidnightBlue;
draw (        0.5986442999u,        0.9690559255u)--(        0.5498380188u,        0.9403307072u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5159478190u,        0.9857030072u)--(        0.5498380188u,        0.9403307072u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5986442999u,        0.9690559255u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5986442999u,        0.9690559255u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5159478190u,        0.9857030072u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5159478190u,        0.9857030072u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5498380188u,        0.9403307072u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5498380188u,        0.9403307072u) withcolor MidnightBlue;
draw (        0.5986442999u,        0.9690559255u)--(        0.5647541001u,        1.0144282255u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5159478190u,        0.9857030072u)--(        0.5647541001u,        1.0144282255u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5986442999u,        0.9690559255u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5986442999u,        0.9690559255u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5159478190u,        0.9857030072u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5159478190u,        0.9857030072u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5647541001u,        1.0144282255u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5647541001u,        1.0144282255u) withcolor MidnightBlue;
draw (        0.5986442999u,        0.9690559255u)--(        0.5647579554u,        1.0144311049u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5616404133u,        1.0709773000u)--(        0.5647579554u,        1.0144311049u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5986442999u,        0.9690559255u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5986442999u,        0.9690559255u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5616404133u,        1.0709773000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5616404133u,        1.0709773000u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5647579554u,        1.0144311049u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5647579554u,        1.0144311049u) withcolor MidnightBlue;
draw (        0.5986442999u,        0.9690559255u)--(        0.5955267577u,        1.0256021206u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5616404133u,        1.0709773000u)--(        0.5955267577u,        1.0256021206u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5986442999u,        0.9690559255u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5986442999u,        0.9690559255u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5616404133u,        1.0709773000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5616404133u,        1.0709773000u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5955267577u,        1.0256021206u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5955267577u,        1.0256021206u) withcolor MidnightBlue;
draw (        0.6413140162u,        1.0589366935u)--(        0.6900000000u,        1.0300000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6413140162u,        1.0589366935u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6413140162u,        1.0589366935u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6900000000u,        1.0300000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6900000000u,        1.0300000000u) withcolor Yellow;
draw (        0.6413140162u,        1.0589366935u)--(        0.5955301218u,        1.0256046331u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5616404133u,        1.0709773000u)--(        0.5955301218u,        1.0256046331u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6413140162u,        1.0589366935u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6413140162u,        1.0589366935u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5616404133u,        1.0709773000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5616404133u,        1.0709773000u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5955301218u,        1.0256046331u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5955301218u,        1.0256046331u) withcolor MidnightBlue;
draw (        0.6413140162u,        1.0589366935u)--(        0.6074243077u,        1.1043093604u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5616404133u,        1.0709773000u)--(        0.6074243077u,        1.1043093604u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6413140162u,        1.0589366935u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6413140162u,        1.0589366935u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5616404133u,        1.0709773000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5616404133u,        1.0709773000u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6074243077u,        1.1043093604u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6074243077u,        1.1043093604u) withcolor Green;
draw (        0.4658867848u,        1.1369830188u)--(        0.4999977181u,        1.0917764322u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4574861055u,        1.0543604725u)--(        0.4999977181u,        1.0917764322u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4658867848u,        1.1369830188u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4658867848u,        1.1369830188u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4574861055u,        1.0543604725u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4574861055u,        1.0543604725u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4999977181u,        1.0917764322u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4999977181u,        1.0917764322u) withcolor Green;
draw (        0.4658867848u,        1.1369830188u)--(        0.5201400645u,        1.1207419324u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5350701699u,        1.1753705256u)--(        0.5201400645u,        1.1207419324u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4658867848u,        1.1369830188u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4658867848u,        1.1369830188u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5350701699u,        1.1753705256u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5350701699u,        1.1753705256u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5201400645u,        1.1207419324u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5201400645u,        1.1207419324u) withcolor Green;
draw (        0.4658867848u,        1.1369830188u)--(        0.4265799037u,        1.1777526195u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3943743593u,        1.1311694234u)--(        0.4265799037u,        1.1777526195u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4658867848u,        1.1369830188u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4658867848u,        1.1369830188u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3943743593u,        1.1311694234u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3943743593u,        1.1311694234u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4265799037u,        1.1777526195u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4265799037u,        1.1777526195u) withcolor Green;
draw (        0.5342012240u,        0.8858843934u)--(        0.5250000000u,        0.8300000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5342012240u,        0.8858843934u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5342012240u,        0.8858843934u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5250000000u,        0.8300000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5250000000u,        0.8300000000u) withcolor Yellow;
draw (        0.5342012240u,        0.8858843934u)--(        0.5003262982u,        0.9312680980u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5159478190u,        0.9857030072u)--(        0.5003262982u,        0.9312680980u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5342012240u,        0.8858843934u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5342012240u,        0.8858843934u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5159478190u,        0.9857030072u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5159478190u,        0.9857030072u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5003262982u,        0.9312680980u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5003262982u,        0.9312680980u) withcolor MidnightBlue;
draw (        0.5342012240u,        0.8858843934u)--(        0.5498227448u,        0.9403193025u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5159478190u,        0.9857030072u)--(        0.5498227448u,        0.9403193025u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5342012240u,        0.8858843934u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5342012240u,        0.8858843934u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5159478190u,        0.9857030072u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5159478190u,        0.9857030072u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5498227448u,        0.9403193025u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5498227448u,        0.9403193025u) withcolor MidnightBlue;
draw (        0.5342012240u,        0.8858843934u)--(        0.5003227291u,        0.9312654338u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4540169597u,        0.8986622596u)--(        0.5003227291u,        0.9312654338u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5342012240u,        0.8858843934u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5342012240u,        0.8858843934u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4540169597u,        0.8986622596u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4540169597u,        0.8986622596u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5003227291u,        0.9312654338u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5003227291u,        0.9312654338u) withcolor MidnightBlue;
draw (        0.5159478190u,        0.9857030072u)--(        0.5128356527u,        1.0422494983u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5616404133u,        1.0709773000u)--(        0.5128356527u,        1.0422494983u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5159478190u,        0.9857030072u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5159478190u,        0.9857030072u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5616404133u,        1.0709773000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5616404133u,        1.0709773000u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5128356527u,        1.0422494983u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5128356527u,        1.0422494983u) withcolor MidnightBlue;
draw (        0.5159478190u,        0.9857030072u)--(        0.5647525796u,        1.0144308088u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5616404133u,        1.0709773000u)--(        0.5647525796u,        1.0144308088u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5159478190u,        0.9857030072u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5159478190u,        0.9857030072u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5616404133u,        1.0709773000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5616404133u,        1.0709773000u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5647525796u,        1.0144308088u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5647525796u,        1.0144308088u) withcolor MidnightBlue;
draw (        0.5159478190u,        0.9857030072u)--(        0.5003194736u,        0.9312700570u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4540169597u,        0.8986622596u)--(        0.5003194736u,        0.9312700570u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5159478190u,        0.9857030072u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5159478190u,        0.9857030072u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4540169597u,        0.8986622596u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4540169597u,        0.8986622596u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5003194736u,        0.9312700570u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5003194736u,        0.9312700570u) withcolor MidnightBlue;
draw (        0.5159478190u,        0.9857030072u)--(        0.4696453052u,        0.9530952097u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4540169597u,        0.8986622596u)--(        0.4696453052u,        0.9530952097u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5159478190u,        0.9857030072u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5159478190u,        0.9857030072u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4540169597u,        0.8986622596u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4540169597u,        0.8986622596u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4696453052u,        0.9530952097u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4696453052u,        0.9530952097u) withcolor MidnightBlue;
draw (        0.5159478190u,        0.9857030072u)--(        0.4696379242u,        0.9531056932u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4153400960u,        0.9691972157u)--(        0.4696379242u,        0.9531056932u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5159478190u,        0.9857030072u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5159478190u,        0.9857030072u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4153400960u,        0.9691972157u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4153400960u,        0.9691972157u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4696379242u,        0.9531056932u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4696379242u,        0.9531056932u) withcolor MidnightBlue;
draw (        0.5159478190u,        0.9857030072u)--(        0.5128076935u,        1.0422479526u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4574861055u,        1.0543604725u)--(        0.5128076935u,        1.0422479526u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5159478190u,        0.9857030072u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5159478190u,        0.9857030072u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4574861055u,        1.0543604725u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4574861055u,        1.0543604725u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5128076935u,        1.0422479526u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5128076935u,        1.0422479526u) withcolor MidnightBlue;
draw (        0.5616404133u,        1.0709773000u)--(        0.5900000000u,        1.1200000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5616404133u,        1.0709773000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5616404133u,        1.0709773000u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5900000000u,        1.1200000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5900000000u,        1.1200000000u) withcolor Yellow;
draw (        0.5616404133u,        1.0709773000u)--(        0.5128155159u,        1.0422837355u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4574861055u,        1.0543604725u)--(        0.5128155159u,        1.0422837355u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5616404133u,        1.0709773000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5616404133u,        1.0709773000u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4574861055u,        1.0543604725u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4574861055u,        1.0543604725u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5128155159u,        1.0422837355u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5128155159u,        1.0422837355u) withcolor MidnightBlue;
draw (        0.5173512273u,        0.7693300003u)--(        0.5650000000u,        0.8000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5173512273u,        0.7693300003u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5173512273u,        0.7693300003u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5650000000u,        0.8000000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5650000000u,        0.8000000000u) withcolor Yellow;
draw (        0.4540169597u,        0.8986622596u)--(        0.4696356389u,        0.9530979841u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4153400960u,        0.9691972157u)--(        0.4696356389u,        0.9530979841u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4540169597u,        0.8986622596u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4540169597u,        0.8986622596u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4153400960u,        0.9691972157u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4153400960u,        0.9691972157u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4696356389u,        0.9530979841u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4696356389u,        0.9530979841u) withcolor MidnightBlue;
draw (        0.4540169597u,        0.8986622596u)--(        0.4258866491u,        0.8495106892u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4331603000u,        0.7933476654u)--(        0.4258866491u,        0.8495106892u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4540169597u,        0.8986622596u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4540169597u,        0.8986622596u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4331603000u,        0.7933476654u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4331603000u,        0.7933476654u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4258866491u,        0.8495106892u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4258866491u,        0.8495106892u) withcolor MidnightBlue;
draw (        0.4540169597u,        0.8986622596u)--(        0.4258737276u,        0.8495180867u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3726317461u,        0.8688178962u)--(        0.4258737276u,        0.8495180867u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4540169597u,        0.8986622596u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4540169597u,        0.8986622596u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3726317461u,        0.8688178962u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3726317461u,        0.8688178962u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4258737276u,        0.8495180867u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4258737276u,        0.8495180867u) withcolor MidnightBlue;
draw (        0.4153400960u,        0.9691972157u)--(        0.4087898101u,        1.0254491935u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4574861055u,        1.0543604725u)--(        0.4087898101u,        1.0254491935u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4153400960u,        0.9691972157u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4153400960u,        0.9691972157u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4574861055u,        1.0543604725u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4574861055u,        1.0543604725u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4087898101u,        1.0254491935u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4087898101u,        1.0254491935u) withcolor MidnightBlue;
draw (        0.4153400960u,        0.9691972157u)--(        0.3799645356u,        0.9249732293u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3726317461u,        0.8688178962u)--(        0.3799645356u,        0.9249732293u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4153400960u,        0.9691972157u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4153400960u,        0.9691972157u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3726317461u,        0.8688178962u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3726317461u,        0.8688178962u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3799645356u,        0.9249732293u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3799645356u,        0.9249732293u) withcolor MidnightBlue;
draw (        0.4153400960u,        0.9691972157u)--(        0.4087794052u,        1.0254479809u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3562482186u,        1.0466060922u)--(        0.4087794052u,        1.0254479809u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4153400960u,        0.9691972157u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4153400960u,        0.9691972157u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3562482186u,        1.0466060922u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3562482186u,        1.0466060922u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4087794052u,        1.0254479809u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4087794052u,        1.0254479809u) withcolor MidnightBlue;
draw (        0.4153400960u,        0.9691972157u)--(        0.3799564157u,        0.9249797258u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3287131465u,        0.9490902755u)--(        0.3799564157u,        0.9249797258u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4153400960u,        0.9691972157u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4153400960u,        0.9691972157u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3287131465u,        0.9490902755u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3287131465u,        0.9490902755u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3799564157u,        0.9249797258u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3799564157u,        0.9249797258u) withcolor MidnightBlue;
draw (        0.4574861055u,        1.0543604725u)--(        0.4087839015u,        1.0254591477u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3562482186u,        1.0466060922u)--(        0.4087839015u,        1.0254591477u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4574861055u,        1.0543604725u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4574861055u,        1.0543604725u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3562482186u,        1.0466060922u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3562482186u,        1.0466060922u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4087839015u,        1.0254591477u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4087839015u,        1.0254591477u) withcolor MidnightBlue;
draw (        0.4574861055u,        1.0543604725u)--(        0.4049504225u,        1.0755074170u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3562482186u,        1.0466060922u)--(        0.4049504225u,        1.0755074170u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4574861055u,        1.0543604725u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4574861055u,        1.0543604725u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3562482186u,        1.0466060922u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3562482186u,        1.0466060922u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4049504225u,        1.0755074170u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4049504225u,        1.0755074170u) withcolor MidnightBlue;
draw (        0.4574861055u,        1.0543604725u)--(        0.4049618716u,        1.0755358380u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3943743593u,        1.1311694234u)--(        0.4049618716u,        1.0755358380u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4574861055u,        1.0543604725u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4574861055u,        1.0543604725u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3943743593u,        1.1311694234u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3943743593u,        1.1311694234u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4049618716u,        1.0755358380u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4049618716u,        1.0755358380u) withcolor MidnightBlue;
draw (        0.4331603000u,        0.7933476654u)--(        0.4258702756u,        0.8495085663u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3726317461u,        0.8688178962u)--(        0.4258702756u,        0.8495085663u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4331603000u,        0.7933476654u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4331603000u,        0.7933476654u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3726317461u,        0.8688178962u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3726317461u,        0.8688178962u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4258702756u,        0.8495085663u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4258702756u,        0.8495085663u) withcolor MidnightBlue;
draw (        0.3726317461u,        0.8688178962u)--(        0.3550000000u,        0.8150000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3726317461u,        0.8688178962u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3726317461u,        0.8688178962u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3550000000u,        0.8150000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3550000000u,        0.8150000000u) withcolor Yellow;
draw (        0.3726317461u,        0.8688178962u)--(        0.3799540053u,        0.9249746034u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3287131465u,        0.9490902755u)--(        0.3799540053u,        0.9249746034u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3726317461u,        0.8688178962u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3726317461u,        0.8688178962u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3287131465u,        0.9490902755u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3287131465u,        0.9490902755u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3799540053u,        0.9249746034u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3799540053u,        0.9249746034u) withcolor MidnightBlue;
draw (        0.3562482186u,        1.0466060922u)--(        0.4049364198u,        1.0755310004u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3943743593u,        1.1311694234u)--(        0.4049364198u,        1.0755310004u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3562482186u,        1.0466060922u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3562482186u,        1.0466060922u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3943743593u,        1.1311694234u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3943743593u,        1.1311694234u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4049364198u,        1.0755310004u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4049364198u,        1.0755310004u) withcolor MidnightBlue;
draw (        0.3562482186u,        1.0466060922u)--(        0.2997501605u,        1.0504997631u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2766607852u,        1.1022111932u)--(        0.2997501605u,        1.0504997631u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3562482186u,        1.0466060922u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3562482186u,        1.0466060922u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2766607852u,        1.1022111932u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2766607852u,        1.1022111932u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2997501605u,        1.0504997631u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2997501605u,        1.0504997631u) withcolor Green;
draw (        0.3562482186u,        1.0466060922u)--(        0.3331588434u,        1.0983175222u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2766607852u,        1.1022111932u)--(        0.3331588434u,        1.0983175222u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3562482186u,        1.0466060922u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3562482186u,        1.0466060922u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2766607852u,        1.1022111932u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2766607852u,        1.1022111932u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3331588434u,        1.0983175222u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3331588434u,        1.0983175222u) withcolor MidnightBlue;
draw (        0.3562482186u,        1.0466060922u)--(        0.3331780126u,        1.0983260771u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3255095011u,        1.1544365501u)--(        0.3331780126u,        1.0983260771u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3562482186u,        1.0466060922u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3562482186u,        1.0466060922u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3255095011u,        1.1544365501u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3255095011u,        1.1544365501u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3331780126u,        1.0983260771u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3331780126u,        1.0983260771u) withcolor MidnightBlue;
draw (        0.5350701699u,        1.1753705256u)--(        0.5500000000u,        1.2300000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5350701699u,        1.1753705256u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5350701699u,        1.1753705256u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5500000000u,        1.2300000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5500000000u,        1.2300000000u) withcolor Yellow;
draw (        0.2766607852u,        1.1022111932u)--(        0.3331592562u,        1.0983235169u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3255095011u,        1.1544365501u)--(        0.3331592562u,        1.0983235169u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2766607852u,        1.1022111932u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2766607852u,        1.1022111932u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3255095011u,        1.1544365501u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3255095011u,        1.1544365501u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3331592562u,        1.0983235169u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3331592562u,        1.0983235169u) withcolor MidnightBlue;
draw (        0.2766607852u,        1.1022111932u)--(        0.2200717103u,        1.1000048744u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1706887017u,        1.1277268819u)--(        0.2200717103u,        1.1000048744u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2766607852u,        1.1022111932u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2766607852u,        1.1022111932u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1706887017u,        1.1277268819u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1706887017u,        1.1277268819u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2200717103u,        1.1000048744u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2200717103u,        1.1000048744u) withcolor Green;
draw (        0.2766607852u,        1.1022111932u)--(        0.2272777766u,        1.1299332006u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1706887017u,        1.1277268819u)--(        0.2272777766u,        1.1299332006u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2766607852u,        1.1022111932u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2766607852u,        1.1022111932u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1706887017u,        1.1277268819u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1706887017u,        1.1277268819u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2272777766u,        1.1299332006u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2272777766u,        1.1299332006u) withcolor MidnightBlue;
draw (        0.2766607852u,        1.1022111932u)--(        0.2272784139u,        1.1299343359u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2259489389u,        1.1865507976u)--(        0.2272784139u,        1.1299343359u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2766607852u,        1.1022111932u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2766607852u,        1.1022111932u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2259489389u,        1.1865507976u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2259489389u,        1.1865507976u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2272784139u,        1.1299343359u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2272784139u,        1.1299343359u) withcolor MidnightBlue;
draw (        0.1706887017u,        1.1277268819u)--(        0.1500000000u,        1.0750000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.1706887017u,        1.1277268819u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1706887017u,        1.1277268819u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1500000000u,        1.0750000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1500000000u,        1.0750000000u) withcolor Yellow;
draw (        0.1706887017u,        1.1277268819u)--(        0.2272777329u,        1.1299343199u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2259489389u,        1.1865507976u)--(        0.2272777329u,        1.1299343199u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.1706887017u,        1.1277268819u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1706887017u,        1.1277268819u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2259489389u,        1.1865507976u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2259489389u,        1.1865507976u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2272777329u,        1.1299343199u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2272777329u,        1.1299343199u) withcolor MidnightBlue;
% r =         0.0566320690
endfig;
end;
