prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(20);
u := 20cm;
path pols[];
pair p[];
R := .0705661937u;
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

x1 = .5953220926u;
y1 = .3290231517u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .7136798384u;
y2 = .5112637737u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .7144302956u;
y3 = .6089314979u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .6209127272u;
y4 = .4437675487u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .5925540902u;
y5 = .5595393461u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .6192980787u;
y6 = .6741574197u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .5331108925u;
y7 = .6630871697u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .4985911919u;
y8 = .7939870128u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .5050400342u;
y9 = .9202157046u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = .5848025304u;
y10 = .9825732635u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;

x11 = .6316205802u;
y11 = 1.0656702809u;
draw fullcircle scaled D shifted z11 withcolor RoyalBlue withpen pencircle scaled 1.0;

x12 = .4050817355u;
y12 = .8292637272u;
draw fullcircle scaled D shifted z12 withcolor RoyalBlue withpen pencircle scaled 1.0;

x13 = .4720816091u;
y13 = 1.0464318653u;
draw fullcircle scaled D shifted z13 withcolor RoyalBlue withpen pencircle scaled 1.0;

x14 = .3284911833u;
y14 = .9197602356u;
draw fullcircle scaled D shifted z14 withcolor RoyalBlue withpen pencircle scaled 1.0;

x15 = .4198610712u;
y15 = .9570008549u;
draw fullcircle scaled D shifted z15 withcolor RoyalBlue withpen pencircle scaled 1.0;

x16 = .5083040772u;
y16 = 1.1730490854u;
draw fullcircle scaled D shifted z16 withcolor RoyalBlue withpen pencircle scaled 1.0;

x17 = .3398354990u;
y17 = 1.0539931867u;
draw fullcircle scaled D shifted z17 withcolor RoyalBlue withpen pencircle scaled 1.0;

x18 = .2780509917u;
y18 = 1.1613381810u;
draw fullcircle scaled D shifted z18 withcolor RoyalBlue withpen pencircle scaled 1.0;

x19 = .4010484353u;
y19 = 1.1448158944u;
draw fullcircle scaled D shifted z19 withcolor RoyalBlue withpen pencircle scaled 1.0;

x20 = .1840626172u;
y20 = 1.1368144925u;
draw fullcircle scaled D shifted z20 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.5953220926u,        0.3290231517u)--(        0.6100000000u,        0.2600000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5953220926u,        0.3290231517u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5953220926u,        0.3290231517u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6100000000u,        0.2600000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6100000000u,        0.2600000000u) withcolor Yellow;
draw (        0.5953220926u,        0.3290231517u)--(        0.6250000000u,        0.2650000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5953220926u,        0.3290231517u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5953220926u,        0.3290231517u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6250000000u,        0.2650000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6250000000u,        0.2650000000u) withcolor Yellow;
draw (        0.5953220926u,        0.3290231517u)--(        0.5700113294u,        0.3948938809u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6209127272u,        0.4437675487u)--(        0.5700113294u,        0.3948938809u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5953220926u,        0.3290231517u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5953220926u,        0.3290231517u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6209127272u,        0.4437675487u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6209127272u,        0.4437675487u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5700113294u,        0.3948938809u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5700113294u,        0.3948938809u) withcolor Green;
draw (        0.7136798384u,        0.5112637737u)--(        0.6631183648u,        0.5604890221u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7144302956u,        0.6089314979u)--(        0.6631183648u,        0.5604890221u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7136798384u,        0.5112637737u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7136798384u,        0.5112637737u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7144302956u,        0.6089314979u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7144302956u,        0.6089314979u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6631183648u,        0.5604890221u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6631183648u,        0.5604890221u) withcolor MidnightBlue;
draw (        0.7136798384u,        0.5112637737u)--(        0.6631139541u,        0.5604844913u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5925540902u,        0.5595393461u)--(        0.6631139541u,        0.5604844913u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7136798384u,        0.5112637737u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7136798384u,        0.5112637737u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5925540902u,        0.5595393461u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5925540902u,        0.5595393461u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6631139541u,        0.5604844913u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6631139541u,        0.5604844913u) withcolor MidnightBlue;
draw (        0.7144302956u,        0.6089314979u)--(        0.7850000000u,        0.6100000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7144302956u,        0.6089314979u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7144302956u,        0.6089314979u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7850000000u,        0.6100000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7850000000u,        0.6100000000u) withcolor Yellow;
draw (        0.7144302956u,        0.6089314979u)--(        0.6631138284u,        0.5604938277u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5925540902u,        0.5595393461u)--(        0.6631138284u,        0.5604938277u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7144302956u,        0.6089314979u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7144302956u,        0.6089314979u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5925540902u,        0.5595393461u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5925540902u,        0.5595393461u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6631138284u,        0.5604938277u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6631138284u,        0.5604938277u) withcolor MidnightBlue;
draw (        0.7144302956u,        0.6089314979u)--(        0.6438705573u,        0.6079770164u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5925540902u,        0.5595393461u)--(        0.6438705573u,        0.6079770164u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7144302956u,        0.6089314979u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7144302956u,        0.6089314979u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5925540902u,        0.5595393461u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5925540902u,        0.5595393461u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6438705573u,        0.6079770164u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6438705573u,        0.6079770164u) withcolor MidnightBlue;
draw (        0.7144302956u,        0.6089314979u)--(        0.6438701505u,        0.6080075823u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6192980787u,        0.6741574197u)--(        0.6438701505u,        0.6080075823u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7144302956u,        0.6089314979u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7144302956u,        0.6089314979u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6192980787u,        0.6741574197u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6192980787u,        0.6741574197u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6438701505u,        0.6080075823u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6438701505u,        0.6080075823u) withcolor MidnightBlue;
draw (        0.6209127272u,        0.4437675487u)--(        0.6600000000u,        0.3850000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6209127272u,        0.4437675487u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6209127272u,        0.4437675487u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6600000000u,        0.3850000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6600000000u,        0.3850000000u) withcolor Yellow;
draw (        0.5925540902u,        0.5595393461u)--(        0.6438490805u,        0.6079997594u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6192980787u,        0.6741574197u)--(        0.6438490805u,        0.6079997594u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5925540902u,        0.5595393461u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5925540902u,        0.5595393461u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6192980787u,        0.6741574197u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6192980787u,        0.6741574197u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6438490805u,        0.6079997594u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6438490805u,        0.6079997594u) withcolor MidnightBlue;
draw (        0.5050400342u,        0.9202157046u)--(        0.5146416470u,        0.9901256233u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5848025304u,        0.9825732635u)--(        0.5146416470u,        0.9901256233u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5050400342u,        0.9202157046u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5050400342u,        0.9202157046u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5848025304u,        0.9825732635u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5848025304u,        0.9825732635u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5146416470u,        0.9901256233u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5146416470u,        0.9901256233u) withcolor MidnightBlue;
draw (        0.5050400342u,        0.9202157046u)--(        0.4413738639u,        0.8897820642u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4050817355u,        0.8292637272u)--(        0.4413738639u,        0.8897820642u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5050400342u,        0.9202157046u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5050400342u,        0.9202157046u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4050817355u,        0.8292637272u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4050817355u,        0.8292637272u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4413738639u,        0.8897820642u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4413738639u,        0.8897820642u) withcolor MidnightBlue;
draw (        0.5050400342u,        0.9202157046u)--(        0.5146201444u,        0.9901285731u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4720816091u,        1.0464318653u)--(        0.5146201444u,        0.9901285731u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5050400342u,        0.9202157046u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5050400342u,        0.9202157046u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4720816091u,        1.0464318653u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4720816091u,        1.0464318653u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5146201444u,        0.9901285731u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5146201444u,        0.9901285731u) withcolor MidnightBlue;
draw (        0.5050400342u,        0.9202157046u)--(        0.4413689969u,        0.8897922481u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4198610712u,        0.9570008549u)--(        0.4413689969u,        0.8897922481u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5050400342u,        0.9202157046u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5050400342u,        0.9202157046u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4198610712u,        0.9570008549u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4198610712u,        0.9570008549u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4413689969u,        0.8897922481u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4413689969u,        0.8897922481u) withcolor MidnightBlue;
draw (        0.5848025304u,        0.9825732635u)--(        0.6550000000u,        0.9750000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5848025304u,        0.9825732635u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5848025304u,        0.9825732635u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6550000000u,        0.9750000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6550000000u,        0.9750000000u) withcolor Yellow;
draw (        0.5848025304u,        0.9825732635u)--(        0.5146439005u,        0.9901465293u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4720816091u,        1.0464318653u)--(        0.5146439005u,        0.9901465293u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5848025304u,        0.9825732635u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5848025304u,        0.9825732635u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4720816091u,        1.0464318653u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4720816091u,        1.0464318653u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5146439005u,        0.9901465293u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5146439005u,        0.9901465293u) withcolor MidnightBlue;
draw (        0.4050817355u,        0.8292637272u)--(        0.4413610950u,        0.8897897198u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4198610712u,        0.9570008549u)--(        0.4413610950u,        0.8897897198u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4050817355u,        0.8292637272u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4050817355u,        0.8292637272u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4198610712u,        0.9570008549u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4198610712u,        0.9570008549u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4413610950u,        0.8897897198u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4413610950u,        0.8897897198u) withcolor MidnightBlue;
draw (        0.4720816091u,        1.0464318653u)--(        0.4045712560u,        1.0258906891u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4198610712u,        0.9570008549u)--(        0.4045712560u,        1.0258906891u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4720816091u,        1.0464318653u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4720816091u,        1.0464318653u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4198610712u,        0.9570008549u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4198610712u,        0.9570008549u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4045712560u,        1.0258906891u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4045712560u,        1.0258906891u) withcolor MidnightBlue;
draw (        0.4720816091u,        1.0464318653u)--(        0.4658022024u,        1.1167181140u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5083040772u,        1.1730490854u)--(        0.4658022024u,        1.1167181140u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4720816091u,        1.0464318653u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4720816091u,        1.0464318653u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5083040772u,        1.1730490854u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5083040772u,        1.1730490854u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4658022024u,        1.1167181140u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4658022024u,        1.1167181140u) withcolor MidnightBlue;
draw (        0.4720816091u,        1.0464318653u)--(        0.4045684524u,        1.0258999056u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3398354990u,        1.0539931867u)--(        0.4045684524u,        1.0258999056u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4720816091u,        1.0464318653u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4720816091u,        1.0464318653u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3398354990u,        1.0539931867u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3398354990u,        1.0539931867u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4045684524u,        1.0258999056u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4045684524u,        1.0258999056u) withcolor MidnightBlue;
draw (        0.4720816091u,        1.0464318653u)--(        0.4073486557u,        1.0745251464u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3398354990u,        1.0539931867u)--(        0.4073486557u,        1.0745251464u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4720816091u,        1.0464318653u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4720816091u,        1.0464318653u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3398354990u,        1.0539931867u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3398354990u,        1.0539931867u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4073486557u,        1.0745251464u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4073486557u,        1.0745251464u) withcolor MidnightBlue;
draw (        0.4720816091u,        1.0464318653u)--(        0.4073515289u,        1.0745317660u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4010484353u,        1.1448158944u)--(        0.4073515289u,        1.0745317660u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4720816091u,        1.0464318653u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4720816091u,        1.0464318653u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4010484353u,        1.1448158944u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4010484353u,        1.1448158944u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4073515289u,        1.0745317660u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4073515289u,        1.0745317660u) withcolor MidnightBlue;
draw (        0.4720816091u,        1.0464318653u)--(        0.4657785155u,        1.1167159938u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4010484353u,        1.1448158944u)--(        0.4657785155u,        1.1167159938u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4720816091u,        1.0464318653u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4720816091u,        1.0464318653u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4010484353u,        1.1448158944u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4010484353u,        1.1448158944u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4657785155u,        1.1167159938u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4657785155u,        1.1167159938u) withcolor MidnightBlue;
draw (        0.3284911833u,        0.9197602356u)--(        0.3551324722u,        0.9851041674u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4198610712u,        0.9570008549u)--(        0.3551324722u,        0.9851041674u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3284911833u,        0.9197602356u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3284911833u,        0.9197602356u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4198610712u,        0.9570008549u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4198610712u,        0.9570008549u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3551324722u,        0.9851041674u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3551324722u,        0.9851041674u) withcolor MidnightBlue;
draw (        0.3284911833u,        0.9197602356u)--(        0.3551312428u,        0.9851046686u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3398354990u,        1.0539931867u)--(        0.3551312428u,        0.9851046686u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3284911833u,        0.9197602356u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3284911833u,        0.9197602356u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3398354990u,        1.0539931867u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3398354990u,        1.0539931867u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3551312428u,        0.9851046686u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3551312428u,        0.9851046686u) withcolor MidnightBlue;
draw (        0.4198610712u,        0.9570008549u)--(        0.3551328442u,        0.9851050242u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3398354990u,        1.0539931867u)--(        0.3551328442u,        0.9851050242u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4198610712u,        0.9570008549u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4198610712u,        0.9570008549u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3398354990u,        1.0539931867u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3398354990u,        1.0539931867u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3551328442u,        0.9851050242u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3551328442u,        0.9851050242u) withcolor MidnightBlue;
draw (        0.4198610712u,        0.9570008549u)--(        0.4045637259u,        1.0258890174u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3398354990u,        1.0539931867u)--(        0.4045637259u,        1.0258890174u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4198610712u,        0.9570008549u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4198610712u,        0.9570008549u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3398354990u,        1.0539931867u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3398354990u,        1.0539931867u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4045637259u,        1.0258890174u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4045637259u,        1.0258890174u) withcolor MidnightBlue;
draw (        0.5083040772u,        1.1730490854u)--(        0.5750000000u,        1.1500000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5083040772u,        1.1730490854u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5083040772u,        1.1730490854u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5750000000u,        1.1500000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5750000000u,        1.1500000000u) withcolor Yellow;
draw (        0.5083040772u,        1.1730490854u)--(        0.5500000000u,        1.2300000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5083040772u,        1.1730490854u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5083040772u,        1.1730490854u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5500000000u,        1.2300000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5500000000u,        1.2300000000u) withcolor Yellow;
draw (        0.5083040772u,        1.1730490854u)--(        0.4657850528u,        1.1167310575u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4010484353u,        1.1448158944u)--(        0.4657850528u,        1.1167310575u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5083040772u,        1.1730490854u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5083040772u,        1.1730490854u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4010484353u,        1.1448158944u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4010484353u,        1.1448158944u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4657850528u,        1.1167310575u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4657850528u,        1.1167310575u) withcolor MidnightBlue;
draw (        0.3398354990u,        1.0539931867u)--(        0.4073467719u,        1.0745313395u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4010484353u,        1.1448158944u)--(        0.4073467719u,        1.0745313395u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3398354990u,        1.0539931867u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3398354990u,        1.0539931867u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4010484353u,        1.1448158944u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4010484353u,        1.1448158944u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4073467719u,        1.0745313395u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4073467719u,        1.0745313395u) withcolor MidnightBlue;
draw (        0.1840626172u,        1.1368144925u)--(        0.1500000000u,        1.0750000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.1840626172u,        1.1368144925u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1840626172u,        1.1368144925u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1500000000u,        1.0750000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1500000000u,        1.0750000000u) withcolor Yellow;
% r =         0.0705661937
endfig;
end;
