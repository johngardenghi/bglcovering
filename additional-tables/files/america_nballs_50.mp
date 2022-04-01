prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(50);
u := 20cm;
path pols[];
pair p[];
R := .0429520825u;
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

x1 = .6035046527u;
y1 = .3021865157u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .5862484713u;
y2 = .3724651797u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .6511916823u;
y3 = .4165541982u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .5926703281u;
y4 = .4487076027u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .6699440384u;
y5 = .4764967695u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .5916752348u;
y6 = .5222716261u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .6390131818u;
y7 = .5157180109u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .7212166072u;
y8 = .5276398376u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .5663095687u;
y9 = .5855409858u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = .6258828090u;
y10 = .5969294055u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;

x11 = .6772533160u;
y11 = .5754527037u;
draw fullcircle scaled D shifted z11 withcolor RoyalBlue withpen pencircle scaled 1.0;

x12 = .7635334415u;
y12 = .5786840888u;
draw fullcircle scaled D shifted z12 withcolor RoyalBlue withpen pencircle scaled 1.0;

x13 = .5434424078u;
y13 = .6454104164u;
draw fullcircle scaled D shifted z13 withcolor RoyalBlue withpen pencircle scaled 1.0;

x14 = .5957248452u;
y14 = .6637303889u;
draw fullcircle scaled D shifted z14 withcolor RoyalBlue withpen pencircle scaled 1.0;

x15 = .6695245620u;
y15 = .6589871315u;
draw fullcircle scaled D shifted z15 withcolor RoyalBlue withpen pencircle scaled 1.0;

x16 = .7214935675u;
y16 = .6134960948u;
draw fullcircle scaled D shifted z16 withcolor RoyalBlue withpen pencircle scaled 1.0;

x17 = .4221353683u;
y17 = .8078313427u;
draw fullcircle scaled D shifted z17 withcolor RoyalBlue withpen pencircle scaled 1.0;

x18 = .4682932178u;
y18 = .7663599950u;
draw fullcircle scaled D shifted z18 withcolor RoyalBlue withpen pencircle scaled 1.0;

x19 = .5321482373u;
y19 = .7192906623u;
draw fullcircle scaled D shifted z19 withcolor RoyalBlue withpen pencircle scaled 1.0;

x20 = .6114321210u;
y20 = .7050001093u;
draw fullcircle scaled D shifted z20 withcolor RoyalBlue withpen pencircle scaled 1.0;

x21 = .3759941786u;
y21 = .8524724606u;
draw fullcircle scaled D shifted z21 withcolor RoyalBlue withpen pencircle scaled 1.0;

x22 = .4437018377u;
y22 = .8857425436u;
draw fullcircle scaled D shifted z22 withcolor RoyalBlue withpen pencircle scaled 1.0;

x23 = .5420285379u;
y23 = .8322101708u;
draw fullcircle scaled D shifted z23 withcolor RoyalBlue withpen pencircle scaled 1.0;

x24 = .3330102677u;
y24 = .9057342245u;
draw fullcircle scaled D shifted z24 withcolor RoyalBlue withpen pencircle scaled 1.0;

x25 = .3994953157u;
y25 = .9304791740u;
draw fullcircle scaled D shifted z25 withcolor RoyalBlue withpen pencircle scaled 1.0;

x26 = .4981429690u;
y26 = .8979876984u;
draw fullcircle scaled D shifted z26 withcolor RoyalBlue withpen pencircle scaled 1.0;

x27 = .5534262774u;
y27 = .9127564013u;
draw fullcircle scaled D shifted z27 withcolor RoyalBlue withpen pencircle scaled 1.0;

x28 = .3252780259u;
y28 = .9688640491u;
draw fullcircle scaled D shifted z28 withcolor RoyalBlue withpen pencircle scaled 1.0;

x29 = .3855642387u;
y29 = .9885731336u;
draw fullcircle scaled D shifted z29 withcolor RoyalBlue withpen pencircle scaled 1.0;

x30 = .4624686007u;
y30 = .9666603525u;
draw fullcircle scaled D shifted z30 withcolor RoyalBlue withpen pencircle scaled 1.0;

x31 = .5663260054u;
y31 = 1.0165327428u;
draw fullcircle scaled D shifted z31 withcolor RoyalBlue withpen pencircle scaled 1.0;

x32 = .5322462144u;
y32 = .9652748060u;
draw fullcircle scaled D shifted z32 withcolor RoyalBlue withpen pencircle scaled 1.0;

x33 = .1724703160u;
y33 = 1.1116137488u;
draw fullcircle scaled D shifted z33 withcolor RoyalBlue withpen pencircle scaled 1.0;

x34 = .2462503031u;
y34 = 1.1281665331u;
draw fullcircle scaled D shifted z34 withcolor RoyalBlue withpen pencircle scaled 1.0;

x35 = .3371139229u;
y35 = 1.0479769349u;
draw fullcircle scaled D shifted z35 withcolor RoyalBlue withpen pencircle scaled 1.0;

x36 = .4191140307u;
y36 = 1.0346686360u;
draw fullcircle scaled D shifted z36 withcolor RoyalBlue withpen pencircle scaled 1.0;

x37 = .4891872607u;
y37 = 1.0321904626u;
draw fullcircle scaled D shifted z37 withcolor RoyalBlue withpen pencircle scaled 1.0;

x38 = .6128542018u;
y38 = .9665997098u;
draw fullcircle scaled D shifted z38 withcolor RoyalBlue withpen pencircle scaled 1.0;

x39 = .6470638850u;
y39 = 1.0286151029u;
draw fullcircle scaled D shifted z39 withcolor RoyalBlue withpen pencircle scaled 1.0;

x40 = .1665533653u;
y40 = 1.1643144095u;
draw fullcircle scaled D shifted z40 withcolor RoyalBlue withpen pencircle scaled 1.0;

x41 = .2336278643u;
y41 = 1.1717113471u;
draw fullcircle scaled D shifted z41 withcolor RoyalBlue withpen pencircle scaled 1.0;

x42 = .3187286688u;
y42 = 1.0950635824u;
draw fullcircle scaled D shifted z42 withcolor RoyalBlue withpen pencircle scaled 1.0;

x43 = .3932999430u;
y43 = 1.0929389529u;
draw fullcircle scaled D shifted z43 withcolor RoyalBlue withpen pencircle scaled 1.0;

x44 = .4673750396u;
y44 = 1.0993719135u;
draw fullcircle scaled D shifted z44 withcolor RoyalBlue withpen pencircle scaled 1.0;

x45 = .5242873119u;
y45 = 1.1508760564u;
draw fullcircle scaled D shifted z45 withcolor RoyalBlue withpen pencircle scaled 1.0;

x46 = .6015930033u;
y46 = 1.0786391807u;
draw fullcircle scaled D shifted z46 withcolor RoyalBlue withpen pencircle scaled 1.0;

x47 = .3046716111u;
y47 = 1.1662707609u;
draw fullcircle scaled D shifted z47 withcolor RoyalBlue withpen pencircle scaled 1.0;

x48 = .3721427304u;
y48 = 1.1553464658u;
draw fullcircle scaled D shifted z48 withcolor RoyalBlue withpen pencircle scaled 1.0;

x49 = .4439876753u;
y49 = 1.1576474332u;
draw fullcircle scaled D shifted z49 withcolor RoyalBlue withpen pencircle scaled 1.0;

x50 = .5584735889u;
y50 = 1.1878280279u;
draw fullcircle scaled D shifted z50 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.6035046527u,        0.3021865157u)--(        0.6250000000u,        0.2650000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6035046527u,        0.3021865157u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6035046527u,        0.3021865157u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6250000000u,        0.2650000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6250000000u,        0.2650000000u) withcolor Yellow;
draw (        0.6035046527u,        0.3021865157u)--(        0.5723997715u,        0.3318069099u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5862484713u,        0.3724651797u)--(        0.5723997715u,        0.3318069099u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6035046527u,        0.3021865157u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6035046527u,        0.3021865157u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5862484713u,        0.3724651797u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5862484713u,        0.3724651797u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5723997715u,        0.3318069099u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5723997715u,        0.3318069099u) withcolor Green;
draw (        0.6035046527u,        0.3021865157u)--(        0.6173533525u,        0.3428447855u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5862484713u,        0.3724651797u)--(        0.6173533525u,        0.3428447855u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6035046527u,        0.3021865157u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6035046527u,        0.3021865157u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5862484713u,        0.3724651797u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5862484713u,        0.3724651797u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6173533525u,        0.3428447855u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6173533525u,        0.3428447855u) withcolor Green;
draw (        0.5862484713u,        0.3724651797u)--(        0.6089185557u,        0.4089473488u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6511916823u,        0.4165541982u)--(        0.6089185557u,        0.4089473488u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5862484713u,        0.3724651797u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5862484713u,        0.3724651797u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6511916823u,        0.4165541982u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6511916823u,        0.4165541982u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6089185557u,        0.4089473488u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6089185557u,        0.4089473488u) withcolor MidnightBlue;
draw (        0.5862484713u,        0.3724651797u)--(        0.6285215979u,        0.3800720291u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6511916823u,        0.4165541982u)--(        0.6285215979u,        0.3800720291u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5862484713u,        0.3724651797u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5862484713u,        0.3724651797u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6511916823u,        0.4165541982u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6511916823u,        0.4165541982u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6285215979u,        0.3800720291u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6285215979u,        0.3800720291u) withcolor Green;
draw (        0.5862484713u,        0.3724651797u)--(        0.5700002596u,        0.4122254237u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5926703281u,        0.4487076027u)--(        0.5700002596u,        0.4122254237u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5862484713u,        0.3724651797u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5862484713u,        0.3724651797u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5926703281u,        0.4487076027u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5926703281u,        0.4487076027u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5700002596u,        0.4122254237u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5700002596u,        0.4122254237u) withcolor Green;
draw (        0.5862484713u,        0.3724651797u)--(        0.6089185399u,        0.4089473586u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5926703281u,        0.4487076027u)--(        0.6089185399u,        0.4089473586u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5862484713u,        0.3724651797u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5862484713u,        0.3724651797u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5926703281u,        0.4487076027u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5926703281u,        0.4487076027u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6089185399u,        0.4089473586u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6089185399u,        0.4089473586u) withcolor MidnightBlue;
draw (        0.6511916823u,        0.4165541982u)--(        0.6089185530u,        0.4089473640u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5926703281u,        0.4487076027u)--(        0.6089185530u,        0.4089473640u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6511916823u,        0.4165541982u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6511916823u,        0.4165541982u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5926703281u,        0.4487076027u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5926703281u,        0.4487076027u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6089185530u,        0.4089473640u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6089185530u,        0.4089473640u) withcolor MidnightBlue;
draw (        0.6511916823u,        0.4165541982u)--(        0.6885346946u,        0.4377763757u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6699440384u,        0.4764967695u)--(        0.6885346946u,        0.4377763757u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6511916823u,        0.4165541982u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6511916823u,        0.4165541982u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6699440384u,        0.4764967695u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6699440384u,        0.4764967695u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6885346946u,        0.4377763757u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6885346946u,        0.4377763757u) withcolor Green;
draw (        0.5926703281u,        0.4487076027u)--(        0.6270398796u,        0.4744683348u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6699440384u,        0.4764967695u)--(        0.6270398796u,        0.4744683348u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5926703281u,        0.4487076027u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5926703281u,        0.4487076027u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6699440384u,        0.4764967695u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6699440384u,        0.4764967695u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6270398796u,        0.4744683348u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6270398796u,        0.4744683348u) withcolor MidnightBlue;
draw (        0.5926703281u,        0.4487076027u)--(        0.5700001070u,        0.4851896868u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5916752348u,        0.5222716261u)--(        0.5700001070u,        0.4851896868u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5926703281u,        0.4487076027u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5926703281u,        0.4487076027u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5916752348u,        0.5222716261u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5916752348u,        0.5222716261u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5700001070u,        0.4851896868u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5700001070u,        0.4851896868u) withcolor Green;
draw (        0.5926703281u,        0.4487076027u)--(        0.6270397130u,        0.4744685571u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6390131818u,        0.5157180109u)--(        0.6270397130u,        0.4744685571u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5926703281u,        0.4487076027u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5926703281u,        0.4487076027u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6390131818u,        0.5157180109u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6390131818u,        0.5157180109u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6270397130u,        0.4744685571u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6270397130u,        0.4744685571u) withcolor MidnightBlue;
draw (        0.6699440384u,        0.4764967695u)--(        0.6270398713u,        0.4744685112u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6390131818u,        0.5157180109u)--(        0.6270398713u,        0.4744685112u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6699440384u,        0.4764967695u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6699440384u,        0.4764967695u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6390131818u,        0.5157180109u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6390131818u,        0.5157180109u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6270398713u,        0.4744685112u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6270398713u,        0.4744685112u) withcolor MidnightBlue;
draw (        0.6699440384u,        0.4764967695u)--(        0.7118960993u,        0.4857112135u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7212166072u,        0.5276398376u)--(        0.7118960993u,        0.4857112135u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6699440384u,        0.4764967695u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6699440384u,        0.4764967695u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7212166072u,        0.5276398376u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7212166072u,        0.5276398376u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7118960993u,        0.4857112135u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7118960993u,        0.4857112135u) withcolor Green;
draw (        0.5916752348u,        0.5222716261u)--(        0.6202388470u,        0.5543496785u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6390131818u,        0.5157180109u)--(        0.6202388470u,        0.5543496785u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5916752348u,        0.5222716261u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5916752348u,        0.5222716261u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6390131818u,        0.5157180109u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6390131818u,        0.5157180109u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6202388470u,        0.5543496785u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6202388470u,        0.5543496785u) withcolor MidnightBlue;
draw (        0.5916752348u,        0.5222716261u)--(        0.5547301395u,        0.5441791893u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5663095687u,        0.5855409858u)--(        0.5547301395u,        0.5441791893u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5916752348u,        0.5222716261u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5916752348u,        0.5222716261u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5663095687u,        0.5855409858u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5663095687u,        0.5855409858u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5547301395u,        0.5441791893u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5547301395u,        0.5441791893u) withcolor Green;
draw (        0.5916752348u,        0.5222716261u)--(        0.6202387530u,        0.5543497622u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6258828090u,        0.5969294055u)--(        0.6202387530u,        0.5543497622u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5916752348u,        0.5222716261u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5916752348u,        0.5222716261u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6258828090u,        0.5969294055u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6258828090u,        0.5969294055u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6202387530u,        0.5543497622u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6202387530u,        0.5543497622u) withcolor MidnightBlue;
draw (        0.6390131818u,        0.5157180109u)--(        0.6785426497u,        0.5325198724u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7212166072u,        0.5276398376u)--(        0.6785426497u,        0.5325198724u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6390131818u,        0.5157180109u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6390131818u,        0.5157180109u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7212166072u,        0.5276398376u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7212166072u,        0.5276398376u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6785426497u,        0.5325198724u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6785426497u,        0.5325198724u) withcolor MidnightBlue;
draw (        0.6390131818u,        0.5157180109u)--(        0.6202389622u,        0.5543497345u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6258828090u,        0.5969294055u)--(        0.6202389622u,        0.5543497345u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6390131818u,        0.5157180109u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6390131818u,        0.5157180109u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6258828090u,        0.5969294055u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6258828090u,        0.5969294055u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6202389622u,        0.5543497345u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6202389622u,        0.5543497345u) withcolor MidnightBlue;
draw (        0.6390131818u,        0.5157180109u)--(        0.6785426057u,        0.5325199759u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6772533160u,        0.5754527037u)--(        0.6785426057u,        0.5325199759u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6390131818u,        0.5157180109u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6390131818u,        0.5157180109u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6772533160u,        0.5754527037u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6772533160u,        0.5754527037u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6785426057u,        0.5325199759u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6785426057u,        0.5325199759u) withcolor MidnightBlue;
draw (        0.7212166072u,        0.5276398376u)--(        0.6785426617u,        0.5325199776u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6772533160u,        0.5754527037u)--(        0.6785426617u,        0.5325199776u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7212166072u,        0.5276398376u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7212166072u,        0.5276398376u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6772533160u,        0.5754527037u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6772533160u,        0.5754527037u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6785426617u,        0.5325199776u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6785426617u,        0.5325199776u) withcolor MidnightBlue;
draw (        0.7212166072u,        0.5276398376u)--(        0.7199272615u,        0.5705725638u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6772533160u,        0.5754527037u)--(        0.7199272615u,        0.5705725638u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7212166072u,        0.5276398376u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7212166072u,        0.5276398376u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6772533160u,        0.5754527037u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6772533160u,        0.5754527037u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7199272615u,        0.5705725638u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7199272615u,        0.5705725638u) withcolor MidnightBlue;
draw (        0.7212166072u,        0.5276398376u)--(        0.7213505680u,        0.5705917112u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7635334415u,        0.5786840888u)--(        0.7213505680u,        0.5705917112u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7212166072u,        0.5276398376u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7212166072u,        0.5276398376u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7635334415u,        0.5786840888u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7635334415u,        0.5786840888u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7213505680u,        0.5705917112u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7213505680u,        0.5705917112u) withcolor MidnightBlue;
draw (        0.7212166072u,        0.5276398376u)--(        0.7633994807u,        0.5357322152u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7635334415u,        0.5786840888u)--(        0.7633994807u,        0.5357322152u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7212166072u,        0.5276398376u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7212166072u,        0.5276398376u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7635334415u,        0.5786840888u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7635334415u,        0.5786840888u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7633994807u,        0.5357322152u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7633994807u,        0.5357322152u) withcolor Green;
draw (        0.7212166072u,        0.5276398376u)--(        0.7199275137u,        0.5705725714u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7214935675u,        0.6134960948u)--(        0.7199275137u,        0.5705725714u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7212166072u,        0.5276398376u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7212166072u,        0.5276398376u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7214935675u,        0.6134960948u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7214935675u,        0.6134960948u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7199275137u,        0.5705725714u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7199275137u,        0.5705725714u) withcolor MidnightBlue;
draw (        0.5663095687u,        0.5855409858u)--(        0.5903848298u,        0.6211115240u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6258828090u,        0.5969294055u)--(        0.5903848298u,        0.6211115240u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5663095687u,        0.5855409858u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5663095687u,        0.5855409858u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6258828090u,        0.5969294055u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6258828090u,        0.5969294055u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5903848298u,        0.6211115240u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5903848298u,        0.6211115240u) withcolor MidnightBlue;
draw (        0.5663095687u,        0.5855409858u)--(        0.5281567661u,        0.6052702800u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5434424078u,        0.6454104164u)--(        0.5281567661u,        0.6052702800u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5663095687u,        0.5855409858u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5663095687u,        0.5855409858u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5434424078u,        0.6454104164u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5434424078u,        0.6454104164u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5281567661u,        0.6052702800u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5281567661u,        0.6052702800u) withcolor Green;
draw (        0.5663095687u,        0.5855409858u)--(        0.5903847869u,        0.6211115531u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5957248452u,        0.6637303889u)--(        0.5903847869u,        0.6211115531u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5663095687u,        0.5855409858u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5663095687u,        0.5855409858u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5957248452u,        0.6637303889u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5957248452u,        0.6637303889u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5903847869u,        0.6211115531u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5903847869u,        0.6211115531u) withcolor MidnightBlue;
draw (        0.6258828090u,        0.5969294055u)--(        0.6641843940u,        0.6163682835u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6772533160u,        0.5754527037u)--(        0.6641843940u,        0.6163682835u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6258828090u,        0.5969294055u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6258828090u,        0.5969294055u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6772533160u,        0.5754527037u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6772533160u,        0.5754527037u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6641843940u,        0.6163682835u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6641843940u,        0.6163682835u) withcolor MidnightBlue;
draw (        0.6258828090u,        0.5969294055u)--(        0.5903848447u,        0.6211115459u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5957248452u,        0.6637303889u)--(        0.5903848447u,        0.6211115459u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6258828090u,        0.5969294055u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6258828090u,        0.5969294055u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5957248452u,        0.6637303889u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5957248452u,        0.6637303889u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5903848447u,        0.6211115459u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5903848447u,        0.6211115459u) withcolor MidnightBlue;
draw (        0.6258828090u,        0.5969294055u)--(        0.6312228096u,        0.6395482485u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5957248452u,        0.6637303889u)--(        0.6312228096u,        0.6395482485u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6258828090u,        0.5969294055u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6258828090u,        0.5969294055u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5957248452u,        0.6637303889u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5957248452u,        0.6637303889u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6312228096u,        0.6395482485u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6312228096u,        0.6395482485u) withcolor MidnightBlue;
draw (        0.6258828090u,        0.5969294055u)--(        0.6312229911u,        0.6395482258u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6695245620u,        0.6589871315u)--(        0.6312229911u,        0.6395482258u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6258828090u,        0.5969294055u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6258828090u,        0.5969294055u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6695245620u,        0.6589871315u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6695245620u,        0.6589871315u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6312229911u,        0.6395482258u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6312229911u,        0.6395482258u) withcolor MidnightBlue;
draw (        0.6258828090u,        0.5969294055u)--(        0.6641843799u,        0.6163683112u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6695245620u,        0.6589871315u)--(        0.6641843799u,        0.6163683112u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6258828090u,        0.5969294055u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6258828090u,        0.5969294055u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6695245620u,        0.6589871315u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6695245620u,        0.6589871315u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6641843799u,        0.6163683112u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6641843799u,        0.6163683112u) withcolor MidnightBlue;
draw (        0.6772533160u,        0.5754527037u)--(        0.6641844524u,        0.6163683021u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6695245620u,        0.6589871315u)--(        0.6641844524u,        0.6163683021u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6772533160u,        0.5754527037u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6772533160u,        0.5754527037u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6695245620u,        0.6589871315u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6695245620u,        0.6589871315u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6641844524u,        0.6163683021u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6641844524u,        0.6163683021u) withcolor MidnightBlue;
draw (        0.6772533160u,        0.5754527037u)--(        0.7199272634u,        0.5705725805u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7214935675u,        0.6134960948u)--(        0.7199272634u,        0.5705725805u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6772533160u,        0.5754527037u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6772533160u,        0.5754527037u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7214935675u,        0.6134960948u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7214935675u,        0.6134960948u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7199272634u,        0.5705725805u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7199272634u,        0.5705725805u) withcolor MidnightBlue;
draw (        0.7635334415u,        0.5786840888u)--(        0.7213597063u,        0.5705442209u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7214935675u,        0.6134960948u)--(        0.7213597063u,        0.5705442209u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7635334415u,        0.5786840888u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7635334415u,        0.5786840888u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7214935675u,        0.6134960948u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7214935675u,        0.6134960948u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7213597063u,        0.5705442209u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7213597063u,        0.5705442209u) withcolor MidnightBlue;
draw (        0.7635334415u,        0.5786840888u)--(        0.7636673027u,        0.6216359627u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7214935675u,        0.6134960948u)--(        0.7636673027u,        0.6216359627u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7635334415u,        0.5786840888u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7635334415u,        0.5786840888u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7214935675u,        0.6134960948u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7214935675u,        0.6134960948u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7636673027u,        0.6216359627u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7636673027u,        0.6216359627u) withcolor Green;
draw (        0.5434424078u,        0.6454104164u)--(        0.5587280611u,        0.6855505483u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5957248452u,        0.6637303889u)--(        0.5587280611u,        0.6855505483u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5434424078u,        0.6454104164u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5434424078u,        0.6454104164u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5957248452u,        0.6637303889u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5957248452u,        0.6637303889u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5587280611u,        0.6855505483u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5587280611u,        0.6855505483u) withcolor MidnightBlue;
draw (        0.5434424078u,        0.6454104164u)--(        0.5587280492u,        0.6855505528u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5321482373u,        0.7192906623u)--(        0.5587280492u,        0.6855505528u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5434424078u,        0.6454104164u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5434424078u,        0.6454104164u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5321482373u,        0.7192906623u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5321482373u,        0.7192906623u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5587280492u,        0.6855505528u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5587280492u,        0.6855505528u) withcolor MidnightBlue;
draw (        0.5957248452u,        0.6637303889u)--(        0.6312229070u,        0.6395483915u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6695245620u,        0.6589871315u)--(        0.6312229070u,        0.6395483915u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5957248452u,        0.6637303889u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5957248452u,        0.6637303889u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6695245620u,        0.6589871315u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6695245620u,        0.6589871315u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6312229070u,        0.6395483915u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6312229070u,        0.6395483915u) withcolor MidnightBlue;
draw (        0.5957248452u,        0.6637303889u)--(        0.5587280765u,        0.6855505744u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5321482373u,        0.7192906623u)--(        0.5587280765u,        0.6855505744u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5957248452u,        0.6637303889u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5957248452u,        0.6637303889u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5321482373u,        0.7192906623u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5321482373u,        0.7192906623u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5587280765u,        0.6855505744u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5587280765u,        0.6855505744u) withcolor MidnightBlue;
draw (        0.5957248452u,        0.6637303889u)--(        0.5691450060u,        0.6974704769u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5321482373u,        0.7192906623u)--(        0.5691450060u,        0.6974704769u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5957248452u,        0.6637303889u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5957248452u,        0.6637303889u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5321482373u,        0.7192906623u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5321482373u,        0.7192906623u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5691450060u,        0.6974704769u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5691450060u,        0.6974704769u) withcolor MidnightBlue;
draw (        0.5957248452u,        0.6637303889u)--(        0.5691451519u,        0.6974705918u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6114321210u,        0.7050001093u)--(        0.5691451519u,        0.6974705918u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5957248452u,        0.6637303889u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5957248452u,        0.6637303889u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6114321210u,        0.7050001093u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6114321210u,        0.7050001093u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5691451519u,        0.6974705918u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5691451519u,        0.6974705918u) withcolor MidnightBlue;
draw (        0.6695245620u,        0.6589871315u)--(        0.7123315493u,        0.6554596363u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7214935675u,        0.6134960948u)--(        0.7123315493u,        0.6554596363u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6695245620u,        0.6589871315u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6695245620u,        0.6589871315u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7214935675u,        0.6134960948u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7214935675u,        0.6134960948u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7123315493u,        0.6554596363u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7123315493u,        0.6554596363u) withcolor Green;
draw (        0.6695245620u,        0.6589871315u)--(        0.6539661989u,        0.6990223501u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6114321210u,        0.7050001093u)--(        0.6539661989u,        0.6990223501u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6695245620u,        0.6589871315u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6695245620u,        0.6589871315u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6114321210u,        0.7050001093u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6114321210u,        0.7050001093u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6539661989u,        0.6990223501u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6539661989u,        0.6990223501u) withcolor Green;
draw (        0.4221353683u,        0.8078313427u)--(        0.4250000000u,        0.7650000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4221353683u,        0.8078313427u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4221353683u,        0.8078313427u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4250000000u,        0.7650000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4250000000u,        0.7650000000u) withcolor Yellow;
draw (        0.4221353683u,        0.8078313427u)--(        0.4253626501u,        0.7650006756u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4682932178u,        0.7663599950u)--(        0.4253626501u,        0.7650006756u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4221353683u,        0.8078313427u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4221353683u,        0.8078313427u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4682932178u,        0.7663599950u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4682932178u,        0.7663599950u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4253626501u,        0.7650006756u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4253626501u,        0.7650006756u) withcolor Green;
draw (        0.4221353683u,        0.8078313427u)--(        0.4189080885u,        0.8506620099u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3759941786u,        0.8524724606u)--(        0.4189080885u,        0.8506620099u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4221353683u,        0.8078313427u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4221353683u,        0.8078313427u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3759941786u,        0.8524724606u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3759941786u,        0.8524724606u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4189080885u,        0.8506620099u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4189080885u,        0.8506620099u) withcolor MidnightBlue;
draw (        0.4221353683u,        0.8078313427u)--(        0.4189170654u,        0.8506626854u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4437018377u,        0.8857425436u)--(        0.4189170654u,        0.8506626854u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4221353683u,        0.8078313427u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4221353683u,        0.8078313427u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4437018377u,        0.8857425436u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4437018377u,        0.8857425436u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4189170654u,        0.8506626854u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4189170654u,        0.8506626854u) withcolor MidnightBlue;
draw (        0.4682932178u,        0.7663599950u)--(        0.4950000000u,        0.8000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4682932178u,        0.7663599950u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4682932178u,        0.7663599950u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4950000000u,        0.8000000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4950000000u,        0.8000000000u) withcolor Yellow;
draw (        0.4682932178u,        0.7663599950u)--(        0.5099999095u,        0.7560919254u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5321482373u,        0.7192906623u)--(        0.5099999095u,        0.7560919254u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4682932178u,        0.7663599950u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4682932178u,        0.7663599950u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5321482373u,        0.7192906623u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5321482373u,        0.7192906623u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5099999095u,        0.7560919254u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5099999095u,        0.7560919254u) withcolor Green;
draw (        0.4682932178u,        0.7663599950u)--(        0.4904415456u,        0.7295587319u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5321482373u,        0.7192906623u)--(        0.4904415456u,        0.7295587319u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4682932178u,        0.7663599950u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4682932178u,        0.7663599950u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5321482373u,        0.7192906623u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5321482373u,        0.7192906623u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4904415456u,        0.7295587319u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4904415456u,        0.7295587319u) withcolor Green;
draw (        0.5321482373u,        0.7192906623u)--(        0.5744352245u,        0.7268200781u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6114321210u,        0.7050001093u)--(        0.5744352245u,        0.7268200781u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5321482373u,        0.7192906623u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5321482373u,        0.7192906623u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6114321210u,        0.7050001093u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6114321210u,        0.7050001093u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5744352245u,        0.7268200781u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5744352245u,        0.7268200781u) withcolor Green;
draw (        0.5321482373u,        0.7192906623u)--(        0.5691451338u,        0.6974706935u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6114321210u,        0.7050001093u)--(        0.5691451338u,        0.6974706935u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5321482373u,        0.7192906623u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5321482373u,        0.7192906623u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6114321210u,        0.7050001093u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6114321210u,        0.7050001093u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5691451338u,        0.6974706935u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5691451338u,        0.6974706935u) withcolor MidnightBlue;
draw (        0.3759941786u,        0.8524724606u)--(        0.3550000000u,        0.8150000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3759941786u,        0.8524724606u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3759941786u,        0.8524724606u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3550000000u,        0.8150000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3550000000u,        0.8150000000u) withcolor Yellow;
draw (        0.3759941786u,        0.8524724606u)--(        0.4007876408u,        0.8875461775u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4437018377u,        0.8857425436u)--(        0.4007876408u,        0.8875461775u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3759941786u,        0.8524724606u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3759941786u,        0.8524724606u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4437018377u,        0.8857425436u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4437018377u,        0.8857425436u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4007876408u,        0.8875461775u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4007876408u,        0.8875461775u) withcolor MidnightBlue;
draw (        0.3759941786u,        0.8524724606u)--(        0.4189083755u,        0.8506688266u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4437018377u,        0.8857425436u)--(        0.4189083755u,        0.8506688266u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3759941786u,        0.8524724606u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3759941786u,        0.8524724606u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4437018377u,        0.8857425436u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4437018377u,        0.8857425436u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4189083755u,        0.8506688266u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4189083755u,        0.8506688266u) withcolor MidnightBlue;
draw (        0.3759941786u,        0.8524724606u)--(        0.3343025663u,        0.8628015871u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3330102677u,        0.9057342245u)--(        0.3343025663u,        0.8628015871u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3759941786u,        0.8524724606u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3759941786u,        0.8524724606u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3330102677u,        0.9057342245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3330102677u,        0.9057342245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3343025663u,        0.8628015871u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3343025663u,        0.8628015871u) withcolor Green;
draw (        0.3759941786u,        0.8524724606u)--(        0.3747018800u,        0.8954050980u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3330102677u,        0.9057342245u)--(        0.3747018800u,        0.8954050980u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3759941786u,        0.8524724606u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3759941786u,        0.8524724606u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3330102677u,        0.9057342245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3330102677u,        0.9057342245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3747018800u,        0.8954050980u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3747018800u,        0.8954050980u) withcolor MidnightBlue;
draw (        0.3759941786u,        0.8524724606u)--(        0.3747023419u,        0.8954051119u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3994953157u,        0.9304791740u)--(        0.3747023419u,        0.8954051119u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3759941786u,        0.8524724606u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3759941786u,        0.8524724606u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3994953157u,        0.9304791740u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3994953157u,        0.9304791740u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3747023419u,        0.8954051119u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3747023419u,        0.8954051119u) withcolor MidnightBlue;
draw (        0.3759941786u,        0.8524724606u)--(        0.4007871525u,        0.8875465227u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3994953157u,        0.9304791740u)--(        0.4007871525u,        0.8875465227u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3759941786u,        0.8524724606u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3759941786u,        0.8524724606u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3994953157u,        0.9304791740u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3994953157u,        0.9304791740u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4007871525u,        0.8875465227u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4007871525u,        0.8875465227u) withcolor MidnightBlue;
draw (        0.4437018377u,        0.8857425436u)--(        0.4007876560u,        0.8875465379u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3994953157u,        0.9304791740u)--(        0.4007876560u,        0.8875465379u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4437018377u,        0.8857425436u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4437018377u,        0.8857425436u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3994953157u,        0.9304791740u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3994953157u,        0.9304791740u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4007876560u,        0.8875465379u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4007876560u,        0.8875465379u) withcolor MidnightBlue;
draw (        0.4437018377u,        0.8857425436u)--(        0.4424094974u,        0.9286751797u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3994953157u,        0.9304791740u)--(        0.4424094974u,        0.9286751797u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4437018377u,        0.8857425436u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4437018377u,        0.8857425436u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3994953157u,        0.9304791740u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3994953157u,        0.9304791740u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4424094974u,        0.9286751797u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4424094974u,        0.9286751797u) withcolor MidnightBlue;
draw (        0.4437018377u,        0.8857425436u)--(        0.4637562232u,        0.9237254741u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4981429690u,        0.8979876984u)--(        0.4637562232u,        0.9237254741u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4437018377u,        0.8857425436u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4437018377u,        0.8857425436u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4981429690u,        0.8979876984u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4981429690u,        0.8979876984u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4637562232u,        0.9237254741u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4637562232u,        0.9237254741u) withcolor MidnightBlue;
draw (        0.4437018377u,        0.8857425436u)--(        0.4780885835u,        0.8600047678u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4981429690u,        0.8979876984u)--(        0.4780885835u,        0.8600047678u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4437018377u,        0.8857425436u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4437018377u,        0.8857425436u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4981429690u,        0.8979876984u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4981429690u,        0.8979876984u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4780885835u,        0.8600047678u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4780885835u,        0.8600047678u) withcolor Green;
draw (        0.4437018377u,        0.8857425436u)--(        0.4424179804u,        0.9286754342u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4624686007u,        0.9666603525u)--(        0.4424179804u,        0.9286754342u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4437018377u,        0.8857425436u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4437018377u,        0.8857425436u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4624686007u,        0.9666603525u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4624686007u,        0.9666603525u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4424179804u,        0.9286754342u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4424179804u,        0.9286754342u) withcolor MidnightBlue;
draw (        0.4437018377u,        0.8857425436u)--(        0.4637524580u,        0.9237274618u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4624686007u,        0.9666603525u)--(        0.4637524580u,        0.9237274618u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4437018377u,        0.8857425436u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4437018377u,        0.8857425436u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4624686007u,        0.9666603525u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4624686007u,        0.9666603525u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4637524580u,        0.9237274618u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4637524580u,        0.9237274618u) withcolor MidnightBlue;
draw (        0.5420285379u,        0.8322101708u)--(        0.5500000000u,        0.7900000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5420285379u,        0.8322101708u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5420285379u,        0.8322101708u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5500000000u,        0.7900000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5500000000u,        0.7900000000u) withcolor Yellow;
draw (        0.5420285379u,        0.8322101708u)--(        0.5340481144u,        0.8744143677u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4981429690u,        0.8979876984u)--(        0.5340481144u,        0.8744143677u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5420285379u,        0.8322101708u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5420285379u,        0.8322101708u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4981429690u,        0.8979876984u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4981429690u,        0.8979876984u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5340481144u,        0.8744143677u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5340481144u,        0.8744143677u) withcolor MidnightBlue;
draw (        0.5420285379u,        0.8322101708u)--(        0.5340620796u,        0.8744170059u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5534262774u,        0.9127564013u)--(        0.5340620796u,        0.8744170059u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5420285379u,        0.8322101708u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5420285379u,        0.8322101708u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5534262774u,        0.9127564013u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5534262774u,        0.9127564013u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5340620796u,        0.8744170059u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5340620796u,        0.8744170059u) withcolor MidnightBlue;
draw (        0.3330102677u,        0.9057342245u)--(        0.3578036317u,        0.9408080109u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3994953157u,        0.9304791740u)--(        0.3578036317u,        0.9408080109u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3330102677u,        0.9057342245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3330102677u,        0.9057342245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3994953157u,        0.9304791740u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3994953157u,        0.9304791740u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3578036317u,        0.9408080109u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3578036317u,        0.9408080109u) withcolor MidnightBlue;
draw (        0.3330102677u,        0.9057342245u)--(        0.3747019518u,        0.8954053876u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3994953157u,        0.9304791740u)--(        0.3747019518u,        0.8954053876u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3330102677u,        0.9057342245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3330102677u,        0.9057342245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3994953157u,        0.9304791740u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3994953157u,        0.9304791740u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3747019518u,        0.8954053876u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3747019518u,        0.8954053876u) withcolor MidnightBlue;
draw (        0.3330102677u,        0.9057342245u)--(        0.3004863325u,        0.9337890819u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3252780259u,        0.9688640491u)--(        0.3004863325u,        0.9337890819u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3330102677u,        0.9057342245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3330102677u,        0.9057342245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3252780259u,        0.9688640491u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3252780259u,        0.9688640491u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3004863325u,        0.9337890819u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3004863325u,        0.9337890819u) withcolor Green;
draw (        0.3330102677u,        0.9057342245u)--(        0.3578019611u,        0.9408091917u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3252780259u,        0.9688640491u)--(        0.3578019611u,        0.9408091917u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3330102677u,        0.9057342245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3330102677u,        0.9057342245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3252780259u,        0.9688640491u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3252780259u,        0.9688640491u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3578019611u,        0.9408091917u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3578019611u,        0.9408091917u) withcolor MidnightBlue;
draw (        0.3994953157u,        0.9304791740u)--(        0.3578047157u,        0.9408123853u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3252780259u,        0.9688640491u)--(        0.3578047157u,        0.9408123853u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3994953157u,        0.9304791740u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3994953157u,        0.9304791740u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3252780259u,        0.9688640491u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3252780259u,        0.9688640491u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3578047157u,        0.9408123853u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3578047157u,        0.9408123853u) withcolor MidnightBlue;
draw (        0.3994953157u,        0.9304791740u)--(        0.4424096919u,        0.9286798106u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4624686007u,        0.9666603525u)--(        0.4424096919u,        0.9286798106u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3994953157u,        0.9304791740u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3994953157u,        0.9304791740u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4624686007u,        0.9666603525u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4624686007u,        0.9666603525u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4424096919u,        0.9286798106u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4424096919u,        0.9286798106u) withcolor MidnightBlue;
draw (        0.4981429690u,        0.8979876984u)--(        0.5340528928u,        0.8744216473u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5534262774u,        0.9127564013u)--(        0.5340528928u,        0.8744216473u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4981429690u,        0.8979876984u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4981429690u,        0.8979876984u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5534262774u,        0.9127564013u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5534262774u,        0.9127564013u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5340528928u,        0.8744216473u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5340528928u,        0.8744216473u) withcolor MidnightBlue;
draw (        0.4981429690u,        0.8979876984u)--(        0.4637578316u,        0.9237276229u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4624686007u,        0.9666603525u)--(        0.4637578316u,        0.9237276229u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4981429690u,        0.8979876984u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4981429690u,        0.8979876984u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4624686007u,        0.9666603525u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4624686007u,        0.9666603525u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4637578316u,        0.9237276229u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4637578316u,        0.9237276229u) withcolor MidnightBlue;
draw (        0.4981429690u,        0.8979876984u)--(        0.4968537381u,        0.9409204280u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4624686007u,        0.9666603525u)--(        0.4968537381u,        0.9409204280u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4981429690u,        0.8979876984u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4981429690u,        0.8979876984u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4624686007u,        0.9666603525u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4624686007u,        0.9666603525u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4968537381u,        0.9409204280u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4968537381u,        0.9409204280u) withcolor MidnightBlue;
draw (        0.4981429690u,        0.8979876984u)--(        0.4968659455u,        0.9409207928u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5322462144u,        0.9652748060u)--(        0.4968659455u,        0.9409207928u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4981429690u,        0.8979876984u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4981429690u,        0.8979876984u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5322462144u,        0.9652748060u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5322462144u,        0.9652748060u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4968659455u,        0.9409207928u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4968659455u,        0.9409207928u) withcolor MidnightBlue;
draw (        0.5534262774u,        0.9127564013u)--(        0.5727904751u,        0.9510957966u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5322462144u,        0.9652748060u)--(        0.5727904751u,        0.9510957966u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5534262774u,        0.9127564013u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5534262774u,        0.9127564013u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5322462144u,        0.9652748060u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5322462144u,        0.9652748060u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5727904751u,        0.9510957966u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5727904751u,        0.9510957966u) withcolor MidnightBlue;
draw (        0.5534262774u,        0.9127564013u)--(        0.5727996626u,        0.9510911549u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6128542018u,        0.9665997098u)--(        0.5727996626u,        0.9510911549u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5534262774u,        0.9127564013u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5534262774u,        0.9127564013u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6128542018u,        0.9665997098u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6128542018u,        0.9665997098u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5727996626u,        0.9510911549u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5727996626u,        0.9510911549u) withcolor MidnightBlue;
draw (        0.3855642387u,        0.9885731336u)--(        0.4283167365u,        0.9927089909u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4624686007u,        0.9666603525u)--(        0.4283167365u,        0.9927089909u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3855642387u,        0.9885731336u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3855642387u,        0.9885731336u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4624686007u,        0.9666603525u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4624686007u,        0.9666603525u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4283167365u,        0.9927089909u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4283167365u,        0.9927089909u) withcolor MidnightBlue;
draw (        0.3855642387u,        0.9885731336u)--(        0.3763620013u,        1.0305278737u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3371139229u,        1.0479769349u)--(        0.3763620013u,        1.0305278737u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3855642387u,        0.9885731336u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3855642387u,        0.9885731336u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3371139229u,        1.0479769349u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3371139229u,        1.0479769349u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3763620013u,        1.0305278737u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3763620013u,        1.0305278737u) withcolor MidnightBlue;
draw (        0.3855642387u,        0.9885731336u)--(        0.3763620075u,        1.0305278751u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4191140307u,        1.0346686360u)--(        0.3763620075u,        1.0305278751u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3855642387u,        0.9885731336u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3855642387u,        0.9885731336u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4191140307u,        1.0346686360u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4191140307u,        1.0346686360u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3763620075u,        1.0305278751u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3763620075u,        1.0305278751u) withcolor MidnightBlue;
draw (        0.3855642387u,        0.9885731336u)--(        0.4283162619u,        0.9927138946u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4191140307u,        1.0346686360u)--(        0.4283162619u,        0.9927138946u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3855642387u,        0.9885731336u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3855642387u,        0.9885731336u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4191140307u,        1.0346686360u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4191140307u,        1.0346686360u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4283162619u,        0.9927138946u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4283162619u,        0.9927138946u) withcolor MidnightBlue;
draw (        0.4624686007u,        0.9666603525u)--(        0.4968602284u,        0.9409291004u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5322462144u,        0.9652748060u)--(        0.4968602284u,        0.9409291004u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4624686007u,        0.9666603525u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4624686007u,        0.9666603525u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5322462144u,        0.9652748060u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5322462144u,        0.9652748060u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4968602284u,        0.9409291004u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4968602284u,        0.9409291004u) withcolor MidnightBlue;
draw (        0.4624686007u,        0.9666603525u)--(        0.4283213248u,        0.9927150054u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4191140307u,        1.0346686360u)--(        0.4283213248u,        0.9927150054u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4624686007u,        0.9666603525u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4624686007u,        0.9666603525u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4191140307u,        1.0346686360u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4191140307u,        1.0346686360u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4283213248u,        0.9927150054u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4283213248u,        0.9927150054u) withcolor MidnightBlue;
draw (        0.4624686007u,        0.9666603525u)--(        0.4532613066u,        1.0086139832u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4191140307u,        1.0346686360u)--(        0.4532613066u,        1.0086139832u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4624686007u,        0.9666603525u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4624686007u,        0.9666603525u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4191140307u,        1.0346686360u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4191140307u,        1.0346686360u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4532613066u,        1.0086139832u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4532613066u,        1.0086139832u) withcolor MidnightBlue;
draw (        0.4624686007u,        0.9666603525u)--(        0.4532813045u,        1.0086183670u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4891872607u,        1.0321904626u)--(        0.4532813045u,        1.0086183670u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4624686007u,        0.9666603525u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4624686007u,        0.9666603525u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4891872607u,        1.0321904626u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4891872607u,        1.0321904626u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4532813045u,        1.0086183670u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4532813045u,        1.0086183670u) withcolor MidnightBlue;
draw (        0.5663260054u,        1.0165327428u)--(        0.5243361173u,        1.0074922403u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5322462144u,        0.9652748060u)--(        0.5243361173u,        1.0074922403u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5663260054u,        1.0165327428u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5663260054u,        1.0165327428u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5322462144u,        0.9652748060u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5322462144u,        0.9652748060u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5243361173u,        1.0074922403u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5243361173u,        1.0074922403u) withcolor MidnightBlue;
draw (        0.5663260054u,        1.0165327428u)--(        0.5243342397u,        1.0075009656u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4891872607u,        1.0321904626u)--(        0.5243342397u,        1.0075009656u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5663260054u,        1.0165327428u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5663260054u,        1.0165327428u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4891872607u,        1.0321904626u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4891872607u,        1.0321904626u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5243342397u,        1.0075009656u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5243342397u,        1.0075009656u) withcolor MidnightBlue;
draw (        0.5663260054u,        1.0165327428u)--(        0.5311790265u,        1.0412222398u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4891872607u,        1.0321904626u)--(        0.5311790265u,        1.0412222398u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5663260054u,        1.0165327428u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5663260054u,        1.0165327428u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4891872607u,        1.0321904626u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4891872607u,        1.0321904626u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5311790265u,        1.0412222398u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5311790265u,        1.0412222398u) withcolor Green;
draw (        0.5663260054u,        1.0165327428u)--(        0.6086728798u,        1.0093477850u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6128542018u,        0.9665997098u)--(        0.6086728798u,        1.0093477850u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5663260054u,        1.0165327428u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5663260054u,        1.0165327428u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6128542018u,        0.9665997098u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6128542018u,        0.9665997098u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6086728798u,        1.0093477850u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6086728798u,        1.0093477850u) withcolor MidnightBlue;
draw (        0.5663260054u,        1.0165327428u)--(        0.6047162952u,        1.0357958431u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6470638850u,        1.0286151029u)--(        0.6047162952u,        1.0357958431u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5663260054u,        1.0165327428u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5663260054u,        1.0165327428u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6470638850u,        1.0286151029u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6470638850u,        1.0286151029u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6047162952u,        1.0357958431u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6047162952u,        1.0357958431u) withcolor MidnightBlue;
draw (        0.5663260054u,        1.0165327428u)--(        0.6086735952u,        1.0093520026u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6470638850u,        1.0286151029u)--(        0.6086735952u,        1.0093520026u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5663260054u,        1.0165327428u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5663260054u,        1.0165327428u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6470638850u,        1.0286151029u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6470638850u,        1.0286151029u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6086735952u,        1.0093520026u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6086735952u,        1.0093520026u) withcolor MidnightBlue;
draw (        0.5663260054u,        1.0165327428u)--(        0.6047138930u,        1.0358006299u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6015930033u,        1.0786391807u)--(        0.6047138930u,        1.0358006299u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5663260054u,        1.0165327428u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5663260054u,        1.0165327428u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6015930033u,        1.0786391807u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6015930033u,        1.0786391807u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6047138930u,        1.0358006299u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6047138930u,        1.0358006299u) withcolor MidnightBlue;
draw (        0.5322462144u,        0.9652748060u)--(        0.5243268933u,        1.0074905110u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4891872607u,        1.0321904626u)--(        0.5243268933u,        1.0074905110u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5322462144u,        0.9652748060u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5322462144u,        0.9652748060u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4891872607u,        1.0321904626u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4891872607u,        1.0321904626u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5243268933u,        1.0074905110u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5243268933u,        1.0074905110u) withcolor MidnightBlue;
draw (        0.5322462144u,        0.9652748060u)--(        0.5727939831u,        0.9511058314u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6128542018u,        0.9665997098u)--(        0.5727939831u,        0.9511058314u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5322462144u,        0.9652748060u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5322462144u,        0.9652748060u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6128542018u,        0.9665997098u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6128542018u,        0.9665997098u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5727939831u,        0.9511058314u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5727939831u,        0.9511058314u) withcolor MidnightBlue;
draw (        0.1724703160u,        1.1116137488u)--(        0.1500000000u,        1.0750000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.1724703160u,        1.1116137488u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1724703160u,        1.1116137488u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1500000000u,        1.0750000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1500000000u,        1.0750000000u) withcolor Yellow;
draw (        0.1724703160u,        1.1116137488u)--(        0.2048979969u,        1.1397798083u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2462503031u,        1.1281665331u)--(        0.2048979969u,        1.1397798083u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.1724703160u,        1.1116137488u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1724703160u,        1.1116137488u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2462503031u,        1.1281665331u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2462503031u,        1.1281665331u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2048979969u,        1.1397798083u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2048979969u,        1.1397798083u) withcolor MidnightBlue;
draw (        0.1724703160u,        1.1116137488u)--(        0.2138226222u,        1.1000004737u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2462503031u,        1.1281665331u)--(        0.2138226222u,        1.1000004737u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.1724703160u,        1.1116137488u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1724703160u,        1.1116137488u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2462503031u,        1.1281665331u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2462503031u,        1.1281665331u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2138226222u,        1.1000004737u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2138226222u,        1.1000004737u) withcolor Green;
draw (        0.1724703160u,        1.1116137488u)--(        0.1359324551u,        1.1341939638u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1665533653u,        1.1643144095u)--(        0.1359324551u,        1.1341939638u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.1724703160u,        1.1116137488u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1724703160u,        1.1116137488u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1665533653u,        1.1643144095u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1665533653u,        1.1643144095u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1359324551u,        1.1341939638u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1359324551u,        1.1341939638u) withcolor Green;
draw (        0.1724703160u,        1.1116137488u)--(        0.2048888170u,        1.1397903736u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2336278643u,        1.1717113471u)--(        0.2048888170u,        1.1397903736u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.1724703160u,        1.1116137488u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1724703160u,        1.1116137488u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2336278643u,        1.1717113471u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2336278643u,        1.1717113471u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2048888170u,        1.1397903736u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2048888170u,        1.1397903736u) withcolor MidnightBlue;
draw (        0.2462503031u,        1.1281665331u)--(        0.2048985132u,        1.1397816466u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2336278643u,        1.1717113471u)--(        0.2048985132u,        1.1397816466u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2462503031u,        1.1281665331u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2462503031u,        1.1281665331u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2336278643u,        1.1717113471u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2336278643u,        1.1717113471u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2048985132u,        1.1397816466u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2048985132u,        1.1397816466u) withcolor MidnightBlue;
draw (        0.2462503031u,        1.1281665331u)--(        0.2891580650u,        1.1262158020u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3187286688u,        1.0950635824u)--(        0.2891580650u,        1.1262158020u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2462503031u,        1.1281665331u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2462503031u,        1.1281665331u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3187286688u,        1.0950635824u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3187286688u,        1.0950635824u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2891580650u,        1.1262158020u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2891580650u,        1.1262158020u) withcolor MidnightBlue;
draw (        0.2462503031u,        1.1281665331u)--(        0.2891581701u,        1.1262181138u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3046716111u,        1.1662707609u)--(        0.2891581701u,        1.1262181138u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2462503031u,        1.1281665331u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2462503031u,        1.1281665331u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3046716111u,        1.1662707609u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3046716111u,        1.1662707609u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2891581701u,        1.1262181138u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2891581701u,        1.1262181138u) withcolor MidnightBlue;
draw (        0.3371139229u,        1.0479769349u)--(        0.3798659471u,        1.0521176854u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4191140307u,        1.0346686360u)--(        0.3798659471u,        1.0521176854u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3371139229u,        1.0479769349u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3371139229u,        1.0479769349u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4191140307u,        1.0346686360u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4191140307u,        1.0346686360u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3798659471u,        1.0521176854u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3798659471u,        1.0521176854u) withcolor MidnightBlue;
draw (        0.3371139229u,        1.0479769349u)--(        0.3763620065u,        1.0305278855u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4191140307u,        1.0346686360u)--(        0.3763620065u,        1.0305278855u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3371139229u,        1.0479769349u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3371139229u,        1.0479769349u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4191140307u,        1.0346686360u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4191140307u,        1.0346686360u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3763620065u,        1.0305278855u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3763620065u,        1.0305278855u) withcolor MidnightBlue;
draw (        0.3371139229u,        1.0479769349u)--(        0.3798635288u,        1.0521425786u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3932999430u,        1.0929389529u)--(        0.3798635288u,        1.0521425786u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3371139229u,        1.0479769349u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3371139229u,        1.0479769349u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3932999430u,        1.0929389529u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3932999430u,        1.0929389529u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3798635288u,        1.0521425786u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3798635288u,        1.0521425786u) withcolor MidnightBlue;
draw (        0.4191140307u,        1.0346686360u)--(        0.4550276900u,        1.0582289940u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4891872607u,        1.0321904626u)--(        0.4550276900u,        1.0582289940u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4191140307u,        1.0346686360u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4191140307u,        1.0346686360u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4891872607u,        1.0321904626u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4891872607u,        1.0321904626u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4550276900u,        1.0582289940u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4550276900u,        1.0582289940u) withcolor MidnightBlue;
draw (        0.4191140307u,        1.0346686360u)--(        0.4532736015u,        1.0086301047u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4891872607u,        1.0321904626u)--(        0.4532736015u,        1.0086301047u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4191140307u,        1.0346686360u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4191140307u,        1.0346686360u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4891872607u,        1.0321904626u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4891872607u,        1.0321904626u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4532736015u,        1.0086301047u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4532736015u,        1.0086301047u) withcolor MidnightBlue;
draw (        0.4191140307u,        1.0346686360u)--(        0.3798752984u,        1.0521387041u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3932999430u,        1.0929389529u)--(        0.3798752984u,        1.0521387041u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4191140307u,        1.0346686360u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4191140307u,        1.0346686360u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3932999430u,        1.0929389529u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3932999430u,        1.0929389529u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3798752984u,        1.0521387041u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3798752984u,        1.0521387041u) withcolor MidnightBlue;
draw (        0.4191140307u,        1.0346686360u)--(        0.4550245574u,        1.0582337685u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4673750396u,        1.0993719135u)--(        0.4550245574u,        1.0582337685u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4191140307u,        1.0346686360u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4191140307u,        1.0346686360u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4673750396u,        1.0993719135u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4673750396u,        1.0993719135u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4550245574u,        1.0582337685u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4550245574u,        1.0582337685u) withcolor MidnightBlue;
draw (        0.4891872607u,        1.0321904626u)--(        0.4550300687u,        1.0582321143u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4673750396u,        1.0993719135u)--(        0.4550300687u,        1.0582321143u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4891872607u,        1.0321904626u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4891872607u,        1.0321904626u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4673750396u,        1.0993719135u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4673750396u,        1.0993719135u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4550300687u,        1.0582321143u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4550300687u,        1.0582321143u) withcolor MidnightBlue;
draw (        0.4891872607u,        1.0321904626u)--(        0.5015322316u,        1.0733302619u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4673750396u,        1.0993719135u)--(        0.5015322316u,        1.0733302619u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4891872607u,        1.0321904626u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4891872607u,        1.0321904626u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4673750396u,        1.0993719135u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4673750396u,        1.0993719135u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5015322316u,        1.0733302619u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5015322316u,        1.0733302619u) withcolor Green;
draw (        0.6128542018u,        0.9665997098u)--(        0.6550000000u,        0.9750000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6128542018u,        0.9665997098u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6128542018u,        0.9665997098u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6550000000u,        0.9750000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6550000000u,        0.9750000000u) withcolor Yellow;
draw (        0.6128542018u,        0.9665997098u)--(        0.6086755793u,        1.0093480489u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6470638850u,        1.0286151029u)--(        0.6086755793u,        1.0093480489u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6128542018u,        0.9665997098u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6128542018u,        0.9665997098u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6470638850u,        1.0286151029u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6470638850u,        1.0286151029u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6086755793u,        1.0093480489u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6086755793u,        1.0093480489u) withcolor MidnightBlue;
draw (        0.6470638850u,        1.0286151029u)--(        0.6900000000u,        1.0300000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6470638850u,        1.0286151029u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6470638850u,        1.0286151029u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6900000000u,        1.0300000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6900000000u,        1.0300000000u) withcolor Yellow;
draw (        0.6470638850u,        1.0286151029u)--(        0.6047171474u,        1.0358008671u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6015930033u,        1.0786391807u)--(        0.6047171474u,        1.0358008671u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6470638850u,        1.0286151029u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6470638850u,        1.0286151029u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6015930033u,        1.0786391807u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6015930033u,        1.0786391807u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6047171474u,        1.0358008671u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6047171474u,        1.0358008671u) withcolor MidnightBlue;
draw (        0.6470638850u,        1.0286151029u)--(        0.6439397408u,        1.0714534165u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6015930033u,        1.0786391807u)--(        0.6439397408u,        1.0714534165u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6470638850u,        1.0286151029u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6470638850u,        1.0286151029u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6015930033u,        1.0786391807u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6015930033u,        1.0786391807u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6439397408u,        1.0714534165u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6439397408u,        1.0714534165u) withcolor Green;
draw (        0.1665533653u,        1.1643144095u)--(        0.1971771349u,        1.1944319481u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2336278643u,        1.1717113471u)--(        0.1971771349u,        1.1944319481u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.1665533653u,        1.1643144095u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1665533653u,        1.1643144095u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2336278643u,        1.1717113471u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2336278643u,        1.1717113471u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1971771349u,        1.1944319481u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1971771349u,        1.1944319481u) withcolor Green;
draw (        0.3187286688u,        1.0950635824u)--(        0.3566208193u,        1.1152889703u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3932999430u,        1.0929389529u)--(        0.3566208193u,        1.1152889703u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3187286688u,        1.0950635824u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3187286688u,        1.0950635824u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3932999430u,        1.0929389529u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3932999430u,        1.0929389529u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3566208193u,        1.1152889703u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3566208193u,        1.1152889703u) withcolor MidnightBlue;
draw (        0.3187286688u,        1.0950635824u)--(        0.2891598253u,        1.1262174728u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3046716111u,        1.1662707609u)--(        0.2891598253u,        1.1262174728u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3187286688u,        1.0950635824u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3187286688u,        1.0950635824u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3046716111u,        1.1662707609u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3046716111u,        1.1662707609u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2891598253u,        1.1262174728u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2891598253u,        1.1262174728u) withcolor MidnightBlue;
draw (        0.3187286688u,        1.0950635824u)--(        0.3342404546u,        1.1351168706u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3046716111u,        1.1662707609u)--(        0.3342404546u,        1.1351168706u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3187286688u,        1.0950635824u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3187286688u,        1.0950635824u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3046716111u,        1.1662707609u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3046716111u,        1.1662707609u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3342404546u,        1.1351168706u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3342404546u,        1.1351168706u) withcolor MidnightBlue;
draw (        0.3187286688u,        1.0950635824u)--(        0.3342560538u,        1.1351108258u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3721427304u,        1.1553464658u)--(        0.3342560538u,        1.1351108258u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3187286688u,        1.0950635824u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3187286688u,        1.0950635824u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3721427304u,        1.1553464658u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3721427304u,        1.1553464658u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3342560538u,        1.1351108258u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3342560538u,        1.1351108258u) withcolor MidnightBlue;
draw (        0.3187286688u,        1.0950635824u)--(        0.3566153454u,        1.1152992224u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3721427304u,        1.1553464658u)--(        0.3566153454u,        1.1152992224u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3187286688u,        1.0950635824u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3187286688u,        1.0950635824u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3721427304u,        1.1553464658u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3721427304u,        1.1553464658u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3566153454u,        1.1152992224u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3566153454u,        1.1152992224u) withcolor MidnightBlue;
draw (        0.3932999430u,        1.0929389529u)--(        0.4284763344u,        1.1175865264u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4673750396u,        1.0993719135u)--(        0.4284763344u,        1.1175865264u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3932999430u,        1.0929389529u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3932999430u,        1.0929389529u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4673750396u,        1.0993719135u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4673750396u,        1.0993719135u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4284763344u,        1.1175865264u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4284763344u,        1.1175865264u) withcolor MidnightBlue;
draw (        0.3932999430u,        1.0929389529u)--(        0.3566248277u,        1.1152955472u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3721427304u,        1.1553464658u)--(        0.3566248277u,        1.1152955472u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3932999430u,        1.0929389529u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3932999430u,        1.0929389529u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3721427304u,        1.1553464658u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3721427304u,        1.1553464658u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3566248277u,        1.1152955472u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3566248277u,        1.1152955472u) withcolor MidnightBlue;
draw (        0.3932999430u,        1.0929389529u)--(        0.4088178458u,        1.1329898716u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3721427304u,        1.1553464658u)--(        0.4088178458u,        1.1329898716u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3932999430u,        1.0929389529u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3932999430u,        1.0929389529u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3721427304u,        1.1553464658u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3721427304u,        1.1553464658u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4088178458u,        1.1329898716u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4088178458u,        1.1329898716u) withcolor MidnightBlue;
draw (        0.3932999430u,        1.0929389529u)--(        0.4088184483u,        1.1329896381u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4439876753u,        1.1576474332u)--(        0.4088184483u,        1.1329896381u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3932999430u,        1.0929389529u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3932999430u,        1.0929389529u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4439876753u,        1.1576474332u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4439876753u,        1.1576474332u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4088184483u,        1.1329896381u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4088184483u,        1.1329896381u) withcolor MidnightBlue;
draw (        0.3932999430u,        1.0929389529u)--(        0.4284691700u,        1.1175967480u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4439876753u,        1.1576474332u)--(        0.4284691700u,        1.1175967480u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3932999430u,        1.0929389529u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3932999430u,        1.0929389529u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4439876753u,        1.1576474332u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4439876753u,        1.1576474332u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4284691700u,        1.1175967480u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4284691700u,        1.1175967480u) withcolor MidnightBlue;
draw (        0.4673750396u,        1.0993719135u)--(        0.4828900849u,        1.1394239392u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5242873119u,        1.1508760564u)--(        0.4828900849u,        1.1394239392u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4673750396u,        1.0993719135u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4673750396u,        1.0993719135u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5242873119u,        1.1508760564u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5242873119u,        1.1508760564u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4828900849u,        1.1394239392u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4828900849u,        1.1394239392u) withcolor MidnightBlue;
draw (        0.4673750396u,        1.0993719135u)--(        0.4284792865u,        1.1175928296u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4439876753u,        1.1576474332u)--(        0.4284792865u,        1.1175928296u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4673750396u,        1.0993719135u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4673750396u,        1.0993719135u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4439876753u,        1.1576474332u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4439876753u,        1.1576474332u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4284792865u,        1.1175928296u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4284792865u,        1.1175928296u) withcolor MidnightBlue;
draw (        0.4673750396u,        1.0993719135u)--(        0.4828834284u,        1.1394265171u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4439876753u,        1.1576474332u)--(        0.4828834284u,        1.1394265171u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4673750396u,        1.0993719135u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4673750396u,        1.0993719135u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4439876753u,        1.1576474332u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4439876753u,        1.1576474332u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4828834284u,        1.1394265171u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4828834284u,        1.1394265171u) withcolor MidnightBlue;
draw (        0.5242873119u,        1.1508760564u)--(        0.4828871656u,        1.1394344971u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4439876753u,        1.1576474332u)--(        0.4828871656u,        1.1394344971u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5242873119u,        1.1508760564u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5242873119u,        1.1508760564u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4439876753u,        1.1576474332u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4439876753u,        1.1576474332u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4828871656u,        1.1394344971u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4828871656u,        1.1394344971u) withcolor MidnightBlue;
draw (        0.6015930033u,        1.0786391807u)--(        0.5900000000u,        1.1200000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6015930033u,        1.0786391807u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6015930033u,        1.0786391807u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5900000000u,        1.1200000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5900000000u,        1.1200000000u) withcolor Yellow;
draw (        0.3046716111u,        1.1662707609u)--(        0.3425656764u,        1.1864925612u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3721427304u,        1.1553464658u)--(        0.3425656764u,        1.1864925612u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3046716111u,        1.1662707609u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3046716111u,        1.1662707609u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3721427304u,        1.1553464658u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3721427304u,        1.1553464658u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3425656764u,        1.1864925612u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3425656764u,        1.1864925612u) withcolor Green;
draw (        0.3046716111u,        1.1662707609u)--(        0.3342486652u,        1.1351246655u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3721427304u,        1.1553464658u)--(        0.3342486652u,        1.1351246655u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3046716111u,        1.1662707609u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3046716111u,        1.1662707609u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3721427304u,        1.1553464658u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3721427304u,        1.1553464658u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3342486652u,        1.1351246655u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3342486652u,        1.1351246655u) withcolor MidnightBlue;
draw (        0.3721427304u,        1.1553464658u)--(        0.4088180499u,        1.1329902064u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4439876753u,        1.1576474332u)--(        0.4088180499u,        1.1329902064u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3721427304u,        1.1553464658u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3721427304u,        1.1553464658u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4439876753u,        1.1576474332u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4439876753u,        1.1576474332u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4088180499u,        1.1329902064u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4088180499u,        1.1329902064u) withcolor MidnightBlue;
draw (        0.5584735889u,        1.1878280279u)--(        0.5500000000u,        1.2300000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5584735889u,        1.1878280279u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5584735889u,        1.1878280279u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5500000000u,        1.2300000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5500000000u,        1.2300000000u) withcolor Yellow;
% r =         0.0429520825
endfig;
end;
