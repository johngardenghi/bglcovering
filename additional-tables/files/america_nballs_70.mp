prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(70);
u := 20cm;
path pols[];
pair p[];
R := .0350628371u;
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

x1 = .7477047925u;
y1 = .5493615527u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .7502440084u;
y2 = .6053104094u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .7005450825u;
y3 = .6410791074u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .6266584787u;
y4 = .9954769032u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .6591878639u;
y5 = 1.0465328070u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .6076157576u;
y6 = 1.0897917912u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .7206992929u;
y7 = .5172116103u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .6943759482u;
y8 = .5812474804u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .6485346317u;
y9 = .6152423693u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = .6556387253u;
y10 = .6793256353u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;

x11 = .5793751086u;
y11 = .9543037800u;
draw fullcircle scaled D shifted z11 withcolor RoyalBlue withpen pencircle scaled 1.0;

x12 = .5626538694u;
y12 = 1.0042615380u;
draw fullcircle scaled D shifted z12 withcolor RoyalBlue withpen pencircle scaled 1.0;

x13 = .5979704306u;
y13 = 1.0480931569u;
draw fullcircle scaled D shifted z13 withcolor RoyalBlue withpen pencircle scaled 1.0;

x14 = .5400277005u;
y14 = 1.1497207711u;
draw fullcircle scaled D shifted z14 withcolor RoyalBlue withpen pencircle scaled 1.0;

x15 = .5454565025u;
y15 = 1.1953242446u;
draw fullcircle scaled D shifted z15 withcolor RoyalBlue withpen pencircle scaled 1.0;

x16 = .6852954712u;
y16 = .4586378910u;
draw fullcircle scaled D shifted z16 withcolor RoyalBlue withpen pencircle scaled 1.0;

x17 = .6781399816u;
y17 = .5217832117u;
draw fullcircle scaled D shifted z17 withcolor RoyalBlue withpen pencircle scaled 1.0;

x18 = .6350640664u;
y18 = .5565481581u;
draw fullcircle scaled D shifted z18 withcolor RoyalBlue withpen pencircle scaled 1.0;

x19 = .6033180725u;
y19 = .6484350862u;
draw fullcircle scaled D shifted z19 withcolor RoyalBlue withpen pencircle scaled 1.0;

x20 = .6044574590u;
y20 = .7101075743u;
draw fullcircle scaled D shifted z20 withcolor RoyalBlue withpen pencircle scaled 1.0;

x21 = .5401502540u;
y21 = .9089901265u;
draw fullcircle scaled D shifted z21 withcolor RoyalBlue withpen pencircle scaled 1.0;

x22 = .5164420824u;
y22 = .9619371268u;
draw fullcircle scaled D shifted z22 withcolor RoyalBlue withpen pencircle scaled 1.0;

x23 = .5067942180u;
y23 = 1.0215338625u;
draw fullcircle scaled D shifted z23 withcolor RoyalBlue withpen pencircle scaled 1.0;

x24 = .4933104308u;
y24 = 1.0778144122u;
draw fullcircle scaled D shifted z24 withcolor RoyalBlue withpen pencircle scaled 1.0;

x25 = .4905454155u;
y25 = 1.1391785883u;
draw fullcircle scaled D shifted z25 withcolor RoyalBlue withpen pencircle scaled 1.0;

x26 = .6321508180u;
y26 = .3800880765u;
draw fullcircle scaled D shifted z26 withcolor RoyalBlue withpen pencircle scaled 1.0;

x27 = .6465400776u;
y27 = .4309216301u;
draw fullcircle scaled D shifted z27 withcolor RoyalBlue withpen pencircle scaled 1.0;

x28 = .6323933737u;
y28 = .4942784759u;
draw fullcircle scaled D shifted z28 withcolor RoyalBlue withpen pencircle scaled 1.0;

x29 = .5918860159u;
y29 = .5872413189u;
draw fullcircle scaled D shifted z29 withcolor RoyalBlue withpen pencircle scaled 1.0;

x30 = .5496759221u;
y30 = .6350349555u;
draw fullcircle scaled D shifted z30 withcolor RoyalBlue withpen pencircle scaled 1.0;

x31 = .5574106762u;
y31 = .6948694045u;
draw fullcircle scaled D shifted z31 withcolor RoyalBlue withpen pencircle scaled 1.0;

x32 = .5357167397u;
y32 = .7889970822u;
draw fullcircle scaled D shifted z32 withcolor RoyalBlue withpen pencircle scaled 1.0;

x33 = .5097125214u;
y33 = .8584841078u;
draw fullcircle scaled D shifted z33 withcolor RoyalBlue withpen pencircle scaled 1.0;

x34 = .4819875584u;
y34 = .9158806770u;
draw fullcircle scaled D shifted z34 withcolor RoyalBlue withpen pencircle scaled 1.0;

x35 = .4610187100u;
y35 = .9778029001u;
draw fullcircle scaled D shifted z35 withcolor RoyalBlue withpen pencircle scaled 1.0;

x36 = .4486059873u;
y36 = 1.0348515887u;
draw fullcircle scaled D shifted z36 withcolor RoyalBlue withpen pencircle scaled 1.0;

x37 = .4416442904u;
y37 = 1.0979364542u;
draw fullcircle scaled D shifted z37 withcolor RoyalBlue withpen pencircle scaled 1.0;

x38 = .4419397385u;
y38 = 1.1609441766u;
draw fullcircle scaled D shifted z38 withcolor RoyalBlue withpen pencircle scaled 1.0;

x39 = .6025643263u;
y39 = .2919635975u;
draw fullcircle scaled D shifted z39 withcolor RoyalBlue withpen pencircle scaled 1.0;

x40 = .5842105156u;
y40 = .3450269674u;
draw fullcircle scaled D shifted z40 withcolor RoyalBlue withpen pencircle scaled 1.0;

x41 = .5899954506u;
y41 = .4059130280u;
draw fullcircle scaled D shifted z41 withcolor RoyalBlue withpen pencircle scaled 1.0;

x42 = .5915574047u;
y42 = .4623706461u;
draw fullcircle scaled D shifted z42 withcolor RoyalBlue withpen pencircle scaled 1.0;

x43 = .5827686159u;
y43 = .5227128852u;
draw fullcircle scaled D shifted z43 withcolor RoyalBlue withpen pencircle scaled 1.0;

x44 = .5466475412u;
y44 = .5780862272u;
draw fullcircle scaled D shifted z44 withcolor RoyalBlue withpen pencircle scaled 1.0;

x45 = .4988577904u;
y45 = .7318268150u;
draw fullcircle scaled D shifted z45 withcolor RoyalBlue withpen pencircle scaled 1.0;

x46 = .4657322242u;
y46 = .7808551920u;
draw fullcircle scaled D shifted z46 withcolor RoyalBlue withpen pencircle scaled 1.0;

x47 = .4481096366u;
y47 = .8772080396u;
draw fullcircle scaled D shifted z47 withcolor RoyalBlue withpen pencircle scaled 1.0;

x48 = .4258451950u;
y48 = .9398604168u;
draw fullcircle scaled D shifted z48 withcolor RoyalBlue withpen pencircle scaled 1.0;

x49 = .4011433295u;
y49 = .9992589086u;
draw fullcircle scaled D shifted z49 withcolor RoyalBlue withpen pencircle scaled 1.0;

x50 = .3982298944u;
y50 = 1.0513327682u;
draw fullcircle scaled D shifted z50 withcolor RoyalBlue withpen pencircle scaled 1.0;

x51 = .3963365383u;
y51 = 1.1169454208u;
draw fullcircle scaled D shifted z51 withcolor RoyalBlue withpen pencircle scaled 1.0;

x52 = .3791771888u;
y52 = 1.1681030405u;
draw fullcircle scaled D shifted z52 withcolor RoyalBlue withpen pencircle scaled 1.0;

x53 = .4079085752u;
y53 = .7876480491u;
draw fullcircle scaled D shifted z53 withcolor RoyalBlue withpen pencircle scaled 1.0;

x54 = .4057978622u;
y54 = .8352155875u;
draw fullcircle scaled D shifted z54 withcolor RoyalBlue withpen pencircle scaled 1.0;

x55 = .3927016165u;
y55 = .8968874432u;
draw fullcircle scaled D shifted z55 withcolor RoyalBlue withpen pencircle scaled 1.0;

x56 = .3667678989u;
y56 = .9562747320u;
draw fullcircle scaled D shifted z56 withcolor RoyalBlue withpen pencircle scaled 1.0;

x57 = .3423779917u;
y57 = 1.0159862313u;
draw fullcircle scaled D shifted z57 withcolor RoyalBlue withpen pencircle scaled 1.0;

x58 = .3627054802u;
y58 = 1.0683119482u;
draw fullcircle scaled D shifted z58 withcolor RoyalBlue withpen pencircle scaled 1.0;

x59 = .3358434634u;
y59 = 1.1189053492u;
draw fullcircle scaled D shifted z59 withcolor RoyalBlue withpen pencircle scaled 1.0;

x60 = .3183502650u;
y60 = 1.1695209439u;
draw fullcircle scaled D shifted z60 withcolor RoyalBlue withpen pencircle scaled 1.0;

x61 = .3564777602u;
y61 = .8497964628u;
draw fullcircle scaled D shifted z61 withcolor RoyalBlue withpen pencircle scaled 1.0;

x62 = .3349268039u;
y62 = .9090214254u;
draw fullcircle scaled D shifted z62 withcolor RoyalBlue withpen pencircle scaled 1.0;

x63 = .3155564997u;
y63 = .9581776778u;
draw fullcircle scaled D shifted z63 withcolor RoyalBlue withpen pencircle scaled 1.0;

x64 = .3089037610u;
y64 = 1.0619846965u;
draw fullcircle scaled D shifted z64 withcolor RoyalBlue withpen pencircle scaled 1.0;

x65 = .2811342264u;
y65 = 1.1185846021u;
draw fullcircle scaled D shifted z65 withcolor RoyalBlue withpen pencircle scaled 1.0;

x66 = .2591941030u;
y66 = 1.1734409720u;
draw fullcircle scaled D shifted z66 withcolor RoyalBlue withpen pencircle scaled 1.0;

x67 = .2231554997u;
y67 = 1.1204689231u;
draw fullcircle scaled D shifted z67 withcolor RoyalBlue withpen pencircle scaled 1.0;

x68 = .2051715531u;
y68 = 1.1799067422u;
draw fullcircle scaled D shifted z68 withcolor RoyalBlue withpen pencircle scaled 1.0;

x69 = .1625938234u;
y69 = 1.1076270503u;
draw fullcircle scaled D shifted z69 withcolor RoyalBlue withpen pencircle scaled 1.0;

x70 = .1642643332u;
y70 = 1.1513019084u;
draw fullcircle scaled D shifted z70 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.7477047925u,        0.5493615527u)--(        0.7700529876u,        0.5763793380u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7502440084u,        0.6053104094u)--(        0.7700529876u,        0.5763793380u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7477047925u,        0.5493615527u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7477047925u,        0.5493615527u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7502440084u,        0.6053104094u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7502440084u,        0.6053104094u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7700529876u,        0.5763793380u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7700529876u,        0.5763793380u) withcolor Green;
draw (        0.7477047925u,        0.5493615527u)--(        0.7557057364u,        0.5152237811u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7206992929u,        0.5172116103u)--(        0.7557057364u,        0.5152237811u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7477047925u,        0.5493615527u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7477047925u,        0.5493615527u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7206992929u,        0.5172116103u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7206992929u,        0.5172116103u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7557057364u,        0.5152237811u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7557057364u,        0.5152237811u) withcolor Green;
draw (        0.7477047925u,        0.5493615527u)--(        0.7126983491u,        0.5513493818u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7206992929u,        0.5172116103u)--(        0.7126983491u,        0.5513493818u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7477047925u,        0.5493615527u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7477047925u,        0.5493615527u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7206992929u,        0.5172116103u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7206992929u,        0.5172116103u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7126983491u,        0.5513493818u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7126983491u,        0.5513493818u) withcolor MidnightBlue;
draw (        0.7477047925u,        0.5493615527u)--(        0.7126985501u,        0.5513529189u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6943759482u,        0.5812474804u)--(        0.7126985501u,        0.5513529189u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7477047925u,        0.5493615527u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7477047925u,        0.5493615527u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6943759482u,        0.5812474804u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6943759482u,        0.5812474804u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7126985501u,        0.5513529189u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7126985501u,        0.5513529189u) withcolor MidnightBlue;
draw (        0.7502440084u,        0.6053104094u)--(        0.7850000000u,        0.6100000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7502440084u,        0.6053104094u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7502440084u,        0.6053104094u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7850000000u,        0.6100000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7850000000u,        0.6100000000u) withcolor Yellow;
draw (        0.7502440084u,        0.6053104094u)--(        0.7154116416u,        0.6093239823u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7005450825u,        0.6410791074u)--(        0.7154116416u,        0.6093239823u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7502440084u,        0.6053104094u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7502440084u,        0.6053104094u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7005450825u,        0.6410791074u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7005450825u,        0.6410791074u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7154116416u,        0.6093239823u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7154116416u,        0.6093239823u) withcolor MidnightBlue;
draw (        0.7502440084u,        0.6053104094u)--(        0.7353774494u,        0.6370655345u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7005450825u,        0.6410791074u)--(        0.7353774494u,        0.6370655345u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7502440084u,        0.6053104094u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7502440084u,        0.6053104094u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7005450825u,        0.6410791074u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7005450825u,        0.6410791074u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7353774494u,        0.6370655345u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7353774494u,        0.6370655345u) withcolor Green;
draw (        0.7502440084u,        0.6053104094u)--(        0.7154090261u,        0.6093012178u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6943759482u,        0.5812474804u)--(        0.7154090261u,        0.6093012178u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7502440084u,        0.6053104094u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7502440084u,        0.6053104094u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6943759482u,        0.5812474804u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6943759482u,        0.5812474804u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7154090261u,        0.6093012178u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7154090261u,        0.6093012178u) withcolor MidnightBlue;
draw (        0.7005450825u,        0.6410791074u)--(        0.7153913263u,        0.6093144793u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6943759482u,        0.5812474804u)--(        0.7153913263u,        0.6093144793u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7005450825u,        0.6410791074u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7005450825u,        0.6410791074u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6943759482u,        0.5812474804u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6943759482u,        0.5812474804u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7153913263u,        0.6093144793u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7153913263u,        0.6093144793u) withcolor MidnightBlue;
draw (        0.7005450825u,        0.6410791074u)--(        0.6657963429u,        0.6457618036u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6485346317u,        0.6152423693u)--(        0.6657963429u,        0.6457618036u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7005450825u,        0.6410791074u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7005450825u,        0.6410791074u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6485346317u,        0.6152423693u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6485346317u,        0.6152423693u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6657963429u,        0.6457618036u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6657963429u,        0.6457618036u) withcolor MidnightBlue;
draw (        0.7005450825u,        0.6410791074u)--(        0.6657969826u,        0.6457665478u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6556387253u,        0.6793256353u)--(        0.6657969826u,        0.6457665478u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7005450825u,        0.6410791074u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7005450825u,        0.6410791074u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6556387253u,        0.6793256353u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6556387253u,        0.6793256353u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6657969826u,        0.6457665478u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6657969826u,        0.6457665478u) withcolor MidnightBlue;
draw (        0.6266584787u,        0.9954769032u)--(        0.6550000000u,        0.9750000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6266584787u,        0.9954769032u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6266584787u,        0.9954769032u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6550000000u,        0.9750000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6550000000u,        0.9750000000u) withcolor Yellow;
draw (        0.6591878639u,        1.0465328070u)--(        0.6900000000u,        1.0300000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6591878639u,        1.0465328070u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6591878639u,        1.0465328070u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6900000000u,        1.0300000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6900000000u,        1.0300000000u) withcolor Yellow;
draw (        0.6076157576u,        1.0897917912u)--(        0.5900000000u,        1.1200000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6076157576u,        1.0897917912u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6076157576u,        1.0897917912u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5900000000u,        1.1200000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5900000000u,        1.1200000000u) withcolor Yellow;
draw (        0.7206992929u,        0.5172116103u)--(        0.7023859191u,        0.5471118258u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6943759482u,        0.5812474804u)--(        0.7023859191u,        0.5471118258u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7206992929u,        0.5172116103u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7206992929u,        0.5172116103u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6943759482u,        0.5812474804u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6943759482u,        0.5812474804u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7023859191u,        0.5471118258u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7023859191u,        0.5471118258u) withcolor MidnightBlue;
draw (        0.7206992929u,        0.5172116103u)--(        0.7126893220u,        0.5513472649u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6943759482u,        0.5812474804u)--(        0.7126893220u,        0.5513472649u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7206992929u,        0.5172116103u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7206992929u,        0.5172116103u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6943759482u,        0.5812474804u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6943759482u,        0.5812474804u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7126893220u,        0.5513472649u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7126893220u,        0.5513472649u) withcolor MidnightBlue;
draw (        0.7206992929u,        0.5172116103u)--(        0.7095332968u,        0.4839742355u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6852954712u,        0.4586378910u)--(        0.7095332968u,        0.4839742355u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7206992929u,        0.5172116103u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7206992929u,        0.5172116103u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6852954712u,        0.4586378910u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6852954712u,        0.4586378910u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7095332968u,        0.4839742355u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7095332968u,        0.4839742355u) withcolor Green;
draw (        0.7206992929u,        0.5172116103u)--(        0.6964614672u,        0.4918752657u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6852954712u,        0.4586378910u)--(        0.6964614672u,        0.4918752657u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7206992929u,        0.5172116103u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7206992929u,        0.5172116103u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6852954712u,        0.4586378910u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6852954712u,        0.4586378910u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6964614672u,        0.4918752657u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6964614672u,        0.4918752657u) withcolor MidnightBlue;
draw (        0.7206992929u,        0.5172116103u)--(        0.6964533759u,        0.4918830087u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6781399816u,        0.5217832117u)--(        0.6964533759u,        0.4918830087u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7206992929u,        0.5172116103u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7206992929u,        0.5172116103u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6781399816u,        0.5217832117u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6781399816u,        0.5217832117u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6964533759u,        0.4918830087u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6964533759u,        0.4918830087u) withcolor MidnightBlue;
draw (        0.7206992929u,        0.5172116103u)--(        0.7023858986u,        0.5471118133u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6781399816u,        0.5217832117u)--(        0.7023858986u,        0.5471118133u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7206992929u,        0.5172116103u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7206992929u,        0.5172116103u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6781399816u,        0.5217832117u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6781399816u,        0.5217832117u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7023858986u,        0.5471118133u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7023858986u,        0.5471118133u) withcolor MidnightBlue;
draw (        0.6943759482u,        0.5812474804u)--(        0.6593188014u,        0.5818791525u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6485346317u,        0.6152423693u)--(        0.6593188014u,        0.5818791525u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6943759482u,        0.5812474804u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6943759482u,        0.5812474804u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6485346317u,        0.6152423693u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6485346317u,        0.6152423693u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6593188014u,        0.5818791525u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6593188014u,        0.5818791525u) withcolor MidnightBlue;
draw (        0.6943759482u,        0.5812474804u)--(        0.7023858921u,        0.5471118195u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6781399816u,        0.5217832117u)--(        0.7023858921u,        0.5471118195u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6943759482u,        0.5812474804u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6943759482u,        0.5812474804u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6781399816u,        0.5217832117u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6781399816u,        0.5217832117u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7023858921u,        0.5471118195u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7023858921u,        0.5471118195u) withcolor MidnightBlue;
draw (        0.6943759482u,        0.5812474804u)--(        0.6701300376u,        0.5559188726u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6781399816u,        0.5217832117u)--(        0.6701300376u,        0.5559188726u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6943759482u,        0.5812474804u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6943759482u,        0.5812474804u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6781399816u,        0.5217832117u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6781399816u,        0.5217832117u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6701300376u,        0.5559188726u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6701300376u,        0.5559188726u) withcolor MidnightBlue;
draw (        0.6943759482u,        0.5812474804u)--(        0.6701214035u,        0.5559271405u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6350640664u,        0.5565481581u)--(        0.6701214035u,        0.5559271405u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6943759482u,        0.5812474804u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6943759482u,        0.5812474804u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6350640664u,        0.5565481581u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6350640664u,        0.5565481581u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6701214035u,        0.5559271405u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6701214035u,        0.5559271405u) withcolor MidnightBlue;
draw (        0.6943759482u,        0.5812474804u)--(        0.6593186111u,        0.5818684980u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6350640664u,        0.5565481581u)--(        0.6593186111u,        0.5818684980u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6943759482u,        0.5812474804u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6943759482u,        0.5812474804u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6350640664u,        0.5565481581u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6350640664u,        0.5565481581u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6593186111u,        0.5818684980u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6593186111u,        0.5818684980u) withcolor MidnightBlue;
draw (        0.6485346317u,        0.6152423693u)--(        0.6383822561u,        0.6488032367u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6556387253u,        0.6793256353u)--(        0.6383822561u,        0.6488032367u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6485346317u,        0.6152423693u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6485346317u,        0.6152423693u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6556387253u,        0.6793256353u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6556387253u,        0.6793256353u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6383822561u,        0.6488032367u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6383822561u,        0.6488032367u) withcolor MidnightBlue;
draw (        0.6485346317u,        0.6152423693u)--(        0.6657911008u,        0.6457647680u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6556387253u,        0.6793256353u)--(        0.6657911008u,        0.6457647680u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6485346317u,        0.6152423693u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6485346317u,        0.6152423693u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6556387253u,        0.6793256353u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6556387253u,        0.6793256353u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6657911008u,        0.6457647680u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6657911008u,        0.6457647680u) withcolor MidnightBlue;
draw (        0.6485346317u,        0.6152423693u)--(        0.6593103397u,        0.5818764186u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6350640664u,        0.5565481581u)--(        0.6593103397u,        0.5818764186u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6485346317u,        0.6152423693u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6485346317u,        0.6152423693u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6350640664u,        0.5565481581u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6350640664u,        0.5565481581u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6593103397u,        0.5818764186u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6593103397u,        0.5818764186u) withcolor MidnightBlue;
draw (        0.6485346317u,        0.6152423693u)--(        0.6134737169u,        0.6148752079u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6033180725u,        0.6484350862u)--(        0.6134737169u,        0.6148752079u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6485346317u,        0.6152423693u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6485346317u,        0.6152423693u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6033180725u,        0.6484350862u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6033180725u,        0.6484350862u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6134737169u,        0.6148752079u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6134737169u,        0.6148752079u) withcolor MidnightBlue;
draw (        0.6485346317u,        0.6152423693u)--(        0.6383789872u,        0.6488022476u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6033180725u,        0.6484350862u)--(        0.6383789872u,        0.6488022476u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6485346317u,        0.6152423693u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6485346317u,        0.6152423693u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6033180725u,        0.6484350862u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6033180725u,        0.6484350862u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6383789872u,        0.6488022476u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6383789872u,        0.6488022476u) withcolor MidnightBlue;
draw (        0.6485346317u,        0.6152423693u)--(        0.6134737666u,        0.6148704986u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5918860159u,        0.5872413189u)--(        0.6134737666u,        0.6148704986u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6485346317u,        0.6152423693u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6485346317u,        0.6152423693u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5918860159u,        0.5872413189u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5918860159u,        0.5872413189u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6134737666u,        0.6148704986u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6134737666u,        0.6148704986u) withcolor MidnightBlue;
draw (        0.6556387253u,        0.6793256353u)--(        0.6383789572u,        0.6488051020u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6033180725u,        0.6484350862u)--(        0.6383789572u,        0.6488051020u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6556387253u,        0.6793256353u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6556387253u,        0.6793256353u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6033180725u,        0.6484350862u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6033180725u,        0.6484350862u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6383789572u,        0.6488051020u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6383789572u,        0.6488051020u) withcolor MidnightBlue;
draw (        0.6556387253u,        0.6793256353u)--(        0.6205778406u,        0.6789556195u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6033180725u,        0.6484350862u)--(        0.6205778406u,        0.6789556195u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6556387253u,        0.6793256353u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6556387253u,        0.6793256353u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6033180725u,        0.6484350862u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6033180725u,        0.6484350862u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6205778406u,        0.6789556195u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6205778406u,        0.6789556195u) withcolor MidnightBlue;
draw (        0.6556387253u,        0.6793256353u)--(        0.6205776905u,        0.6789701261u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6044574590u,        0.7101075743u)--(        0.6205776905u,        0.6789701261u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6556387253u,        0.6793256353u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6556387253u,        0.6793256353u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6044574590u,        0.7101075743u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6044574590u,        0.7101075743u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6205776905u,        0.6789701261u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6205776905u,        0.6789701261u) withcolor MidnightBlue;
draw (        0.5400277005u,        1.1497207711u)--(        0.5750000000u,        1.1500000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5400277005u,        1.1497207711u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5400277005u,        1.1497207711u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5750000000u,        1.1500000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5750000000u,        1.1500000000u) withcolor Yellow;
draw (        0.5454565025u,        1.1953242446u)--(        0.5500000000u,        1.2300000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5454565025u,        1.1953242446u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5454565025u,        1.1953242446u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5500000000u,        1.2300000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5500000000u,        1.2300000000u) withcolor Yellow;
draw (        0.6852954712u,        0.4586378910u)--(        0.6669873787u,        0.4885413406u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6781399816u,        0.5217832117u)--(        0.6669873787u,        0.4885413406u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6852954712u,        0.4586378910u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6852954712u,        0.4586378910u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6781399816u,        0.5217832117u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6781399816u,        0.5217832117u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6669873787u,        0.4885413406u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6669873787u,        0.4885413406u) withcolor MidnightBlue;
draw (        0.6852954712u,        0.4586378910u)--(        0.6964480741u,        0.4918797621u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6781399816u,        0.5217832117u)--(        0.6964480741u,        0.4918797621u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6852954712u,        0.4586378910u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6852954712u,        0.4586378910u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6781399816u,        0.5217832117u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6781399816u,        0.5217832117u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6964480741u,        0.4918797621u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6964480741u,        0.4918797621u) withcolor MidnightBlue;
draw (        0.6852954712u,        0.4586378910u)--(        0.6669832382u,        0.4885388052u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6323933737u,        0.4942784759u)--(        0.6669832382u,        0.4885388052u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6852954712u,        0.4586378910u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6852954712u,        0.4586378910u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6323933737u,        0.4942784759u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6323933737u,        0.4942784759u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6669832382u,        0.4885388052u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6669832382u,        0.4885388052u) withcolor MidnightBlue;
draw (        0.6781399816u,        0.5217832117u)--(        0.6701212189u,        0.5559168021u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6350640664u,        0.5565481581u)--(        0.6701212189u,        0.5559168021u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6781399816u,        0.5217832117u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6781399816u,        0.5217832117u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6350640664u,        0.5565481581u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6350640664u,        0.5565481581u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6701212189u,        0.5559168021u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6701212189u,        0.5559168021u) withcolor MidnightBlue;
draw (        0.6781399816u,        0.5217832117u)--(        0.6669838549u,        0.4885425230u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6323933737u,        0.4942784759u)--(        0.6669838549u,        0.4885425230u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6781399816u,        0.5217832117u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6781399816u,        0.5217832117u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6323933737u,        0.4942784759u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6323933737u,        0.4942784759u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6669838549u,        0.4885425230u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6669838549u,        0.4885425230u) withcolor MidnightBlue;
draw (        0.6350640664u,        0.5565481581u)--(        0.6176735990u,        0.5261019072u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6323933737u,        0.4942784759u)--(        0.6176735990u,        0.5261019072u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6350640664u,        0.5565481581u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6350640664u,        0.5565481581u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6323933737u,        0.4942784759u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6323933737u,        0.4942784759u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6176735990u,        0.5261019072u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6176735990u,        0.5261019072u) withcolor MidnightBlue;
draw (        0.6350640664u,        0.5565481581u)--(        0.6001643681u,        0.5531697554u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5918860159u,        0.5872413189u)--(        0.6001643681u,        0.5531697554u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6350640664u,        0.5565481581u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6350640664u,        0.5565481581u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5918860159u,        0.5872413189u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5918860159u,        0.5872413189u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6001643681u,        0.5531697554u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6001643681u,        0.5531697554u) withcolor MidnightBlue;
draw (        0.6350640664u,        0.5565481581u)--(        0.6176669135u,        0.5261057269u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5827686159u,        0.5227128852u)--(        0.6176669135u,        0.5261057269u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6350640664u,        0.5565481581u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6350640664u,        0.5565481581u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5827686159u,        0.5227128852u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5827686159u,        0.5227128852u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6176669135u,        0.5261057269u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6176669135u,        0.5261057269u) withcolor MidnightBlue;
draw (        0.6350640664u,        0.5565481581u)--(        0.6001657688u,        0.5531553164u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5827686159u,        0.5227128852u)--(        0.6001657688u,        0.5531553164u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6350640664u,        0.5565481581u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6350640664u,        0.5565481581u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5827686159u,        0.5227128852u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5827686159u,        0.5227128852u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6001657688u,        0.5531553164u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6001657688u,        0.5531553164u) withcolor MidnightBlue;
draw (        0.6033180725u,        0.6484350862u)--(        0.6205643699u,        0.6789632335u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6044574590u,        0.7101075743u)--(        0.6205643699u,        0.6789632335u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6033180725u,        0.6484350862u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6033180725u,        0.6484350862u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6044574590u,        0.7101075743u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6044574590u,        0.7101075743u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6205643699u,        0.6789632335u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6205643699u,        0.6789632335u) withcolor MidnightBlue;
draw (        0.6033180725u,        0.6484350862u)--(        0.6134694074u,        0.6148739040u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5918860159u,        0.5872413189u)--(        0.6134694074u,        0.6148739040u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6033180725u,        0.6484350862u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6033180725u,        0.6484350862u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5918860159u,        0.5872413189u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5918860159u,        0.5872413189u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6134694074u,        0.6148739040u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6134694074u,        0.6148739040u) withcolor MidnightBlue;
draw (        0.6033180725u,        0.6484350862u)--(        0.5817346810u,        0.6208025010u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5918860159u,        0.5872413189u)--(        0.5817346810u,        0.6208025010u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6033180725u,        0.6484350862u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6033180725u,        0.6484350862u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5918860159u,        0.5872413189u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5918860159u,        0.5872413189u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5817346810u,        0.6208025010u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5817346810u,        0.6208025010u) withcolor MidnightBlue;
draw (        0.6033180725u,        0.6484350862u)--(        0.5817239768u,        0.6208108653u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5496759221u,        0.6350349555u)--(        0.5817239768u,        0.6208108653u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6033180725u,        0.6484350862u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6033180725u,        0.6484350862u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5496759221u,        0.6350349555u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5496759221u,        0.6350349555u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5817239768u,        0.6208108653u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5817239768u,        0.6208108653u) withcolor MidnightBlue;
draw (        0.6033180725u,        0.6484350862u)--(        0.5712700178u,        0.6626591764u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5496759221u,        0.6350349555u)--(        0.5712700178u,        0.6626591764u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6033180725u,        0.6484350862u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6033180725u,        0.6484350862u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5496759221u,        0.6350349555u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5496759221u,        0.6350349555u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5712700178u,        0.6626591764u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5712700178u,        0.6626591764u) withcolor MidnightBlue;
draw (        0.6033180725u,        0.6484350862u)--(        0.5712715306u,        0.6626625843u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5574106762u,        0.6948694045u)--(        0.5712715306u,        0.6626625843u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6033180725u,        0.6484350862u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6033180725u,        0.6484350862u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5574106762u,        0.6948694045u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5574106762u,        0.6948694045u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5712715306u,        0.6626625843u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5712715306u,        0.6626625843u) withcolor MidnightBlue;
draw (        0.6465400776u,        0.4309216301u)--(        0.6265246654u,        0.4597102706u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6323933737u,        0.4942784759u)--(        0.6265246654u,        0.4597102706u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6465400776u,        0.4309216301u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6465400776u,        0.4309216301u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6323933737u,        0.4942784759u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6323933737u,        0.4942784759u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6265246654u,        0.4597102706u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6265246654u,        0.4597102706u) withcolor MidnightBlue;
draw (        0.6465400776u,        0.4309216301u)--(        0.6115757449u,        0.4335480321u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5899954506u,        0.4059130280u)--(        0.6115757449u,        0.4335480321u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6465400776u,        0.4309216301u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6465400776u,        0.4309216301u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5899954506u,        0.4059130280u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5899954506u,        0.4059130280u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6115757449u,        0.4335480321u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6115757449u,        0.4335480321u) withcolor MidnightBlue;
draw (        0.6465400776u,        0.4309216301u)--(        0.6115786208u,        0.4335860416u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5915574047u,        0.4623706461u)--(        0.6115786208u,        0.4335860416u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6465400776u,        0.4309216301u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6465400776u,        0.4309216301u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5915574047u,        0.4623706461u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5915574047u,        0.4623706461u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6115786208u,        0.4335860416u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6115786208u,        0.4335860416u) withcolor MidnightBlue;
draw (        0.6465400776u,        0.4309216301u)--(        0.6265188614u,        0.4597062345u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5915574047u,        0.4623706461u)--(        0.6265188614u,        0.4597062345u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6465400776u,        0.4309216301u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6465400776u,        0.4309216301u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5915574047u,        0.4623706461u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5915574047u,        0.4623706461u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6265188614u,        0.4597062345u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6265188614u,        0.4597062345u) withcolor MidnightBlue;
draw (        0.6323933737u,        0.4942784759u)--(        0.6265192389u,        0.4597111923u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5915574047u,        0.4623706461u)--(        0.6265192389u,        0.4597111923u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6323933737u,        0.4942784759u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6323933737u,        0.4942784759u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5915574047u,        0.4623706461u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5915574047u,        0.4623706461u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6265192389u,        0.4597111923u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6265192389u,        0.4597111923u) withcolor MidnightBlue;
draw (        0.6323933737u,        0.4942784759u)--(        0.6176675560u,        0.5260991113u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5827686159u,        0.5227128852u)--(        0.6176675560u,        0.5260991113u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6323933737u,        0.4942784759u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6323933737u,        0.4942784759u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5827686159u,        0.5227128852u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5827686159u,        0.5227128852u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6176675560u,        0.5260991113u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6176675560u,        0.5260991113u) withcolor MidnightBlue;
draw (        0.5918860159u,        0.5872413189u)--(        0.5817180227u,        0.6207974579u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5496759221u,        0.6350349555u)--(        0.5817180227u,        0.6207974579u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5918860159u,        0.5872413189u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5918860159u,        0.5872413189u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5496759221u,        0.6350349555u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5496759221u,        0.6350349555u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5817180227u,        0.6207974579u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5817180227u,        0.6207974579u) withcolor MidnightBlue;
draw (        0.5918860159u,        0.5872413189u)--(        0.6001476081u,        0.5531656877u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5827686159u,        0.5227128852u)--(        0.6001476081u,        0.5531656877u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5918860159u,        0.5872413189u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5918860159u,        0.5872413189u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5827686159u,        0.5227128852u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5827686159u,        0.5227128852u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6001476081u,        0.5531656877u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6001476081u,        0.5531656877u) withcolor MidnightBlue;
draw (        0.5918860159u,        0.5872413189u)--(        0.5745070237u,        0.5567885164u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5827686159u,        0.5227128852u)--(        0.5745070237u,        0.5567885164u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5918860159u,        0.5872413189u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5918860159u,        0.5872413189u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5827686159u,        0.5227128852u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5827686159u,        0.5227128852u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5745070237u,        0.5567885164u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5745070237u,        0.5567885164u) withcolor MidnightBlue;
draw (        0.5918860159u,        0.5872413189u)--(        0.5745027663u,        0.5567909464u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5466475412u,        0.5780862272u)--(        0.5745027663u,        0.5567909464u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5918860159u,        0.5872413189u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5918860159u,        0.5872413189u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5466475412u,        0.5780862272u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5466475412u,        0.5780862272u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5745027663u,        0.5567909464u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5745027663u,        0.5567909464u) withcolor MidnightBlue;
draw (        0.5496759221u,        0.6350349555u)--(        0.5712677430u,        0.6626609545u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5574106762u,        0.6948694045u)--(        0.5712677430u,        0.6626609545u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5496759221u,        0.6350349555u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5496759221u,        0.6350349555u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5574106762u,        0.6948694045u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5574106762u,        0.6948694045u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5712677430u,        0.6626609545u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5712677430u,        0.6626609545u) withcolor MidnightBlue;
draw (        0.5097125214u,        0.8584841078u)--(        0.5300000000u,        0.8300000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5097125214u,        0.8584841078u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5097125214u,        0.8584841078u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5300000000u,        0.8300000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5300000000u,        0.8300000000u) withcolor Yellow;
draw (        0.6025643263u,        0.2919635975u)--(        0.6250000000u,        0.2650000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6025643263u,        0.2919635975u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6025643263u,        0.2919635975u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6250000000u,        0.2650000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6250000000u,        0.2650000000u) withcolor Yellow;
draw (        0.5899954506u,        0.4059130280u)--(        0.5700015241u,        0.4347165947u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5915574047u,        0.4623706461u)--(        0.5700015241u,        0.4347165947u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5899954506u,        0.4059130280u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5899954506u,        0.4059130280u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5915574047u,        0.4623706461u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5915574047u,        0.4623706461u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5700015241u,        0.4347165947u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5700015241u,        0.4347165947u) withcolor Green;
draw (        0.5899954506u,        0.4059130280u)--(        0.6115513312u,        0.4335670795u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5915574047u,        0.4623706461u)--(        0.6115513312u,        0.4335670795u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5899954506u,        0.4059130280u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5899954506u,        0.4059130280u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5915574047u,        0.4623706461u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5915574047u,        0.4623706461u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6115513312u,        0.4335670795u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6115513312u,        0.4335670795u) withcolor MidnightBlue;
draw (        0.5827686159u,        0.5227128852u)--(        0.5549166413u,        0.5440124171u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5466475412u,        0.5780862272u)--(        0.5549166413u,        0.5440124171u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5827686159u,        0.5227128852u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5827686159u,        0.5227128852u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5466475412u,        0.5780862272u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5466475412u,        0.5780862272u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5549166413u,        0.5440124171u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5549166413u,        0.5440124171u) withcolor Green;
draw (        0.5827686159u,        0.5227128852u)--(        0.5744995158u,        0.5567866953u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5466475412u,        0.5780862272u)--(        0.5744995158u,        0.5567866953u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5827686159u,        0.5227128852u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5827686159u,        0.5227128852u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5466475412u,        0.5780862272u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5466475412u,        0.5780862272u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5744995158u,        0.5567866953u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5744995158u,        0.5567866953u) withcolor MidnightBlue;
draw (        0.4988577904u,        0.7318268150u)--(        0.5100000000u,        0.7650000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4988577904u,        0.7318268150u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4988577904u,        0.7318268150u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5100000000u,        0.7650000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5100000000u,        0.7650000000u) withcolor Yellow;
draw (        0.4657322242u,        0.7808551920u)--(        0.4950000000u,        0.8000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4657322242u,        0.7808551920u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4657322242u,        0.7808551920u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4950000000u,        0.8000000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4950000000u,        0.8000000000u) withcolor Yellow;
draw (        0.3183502650u,        1.1695209439u)--(        0.2900000000u,        1.1900000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3183502650u,        1.1695209439u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3183502650u,        1.1695209439u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2900000000u,        1.1900000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2900000000u,        1.1900000000u) withcolor Yellow;
draw (        0.3564777602u,        0.8497964628u)--(        0.3600000000u,        0.8150000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3564777602u,        0.8497964628u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3564777602u,        0.8497964628u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3600000000u,        0.8150000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3600000000u,        0.8150000000u) withcolor Yellow;
draw (        0.2591941030u,        1.1734409720u)--(        0.2900000000u,        1.1900000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2591941030u,        1.1734409720u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2591941030u,        1.1734409720u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2900000000u,        1.1900000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2900000000u,        1.1900000000u) withcolor Yellow;
draw (        0.1625938234u,        1.1076270503u)--(        0.1500000000u,        1.0750000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.1625938234u,        1.1076270503u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1625938234u,        1.1076270503u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1500000000u,        1.0750000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1500000000u,        1.0750000000u) withcolor Yellow;
% r =         0.0350628371
endfig;
end;
