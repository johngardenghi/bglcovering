prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(60);
u := 20cm;
path pols[];
pair p[];
R := .0383522944u;
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

x1 = .5990193452u;
y1 = .2926772424u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .5928958498u;
y2 = .3527062786u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .6030319231u;
y3 = .4017013042u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .6561398587u;
y4 = .4177973706u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .6948856339u;
y5 = .4750087115u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .7332615210u;
y6 = .5297479033u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .5891581949u;
y7 = .4531241522u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .6302584835u;
y8 = .4753030747u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .6642295878u;
y9 = .5269368446u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = .6996016511u;
y10 = .5729877217u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;

x11 = .7605091497u;
y11 = .5881588980u;
draw fullcircle scaled D shifted z11 withcolor RoyalBlue withpen pencircle scaled 1.0;

x12 = .5910948730u;
y12 = .5174021761u;
draw fullcircle scaled D shifted z12 withcolor RoyalBlue withpen pencircle scaled 1.0;

x13 = .6151407502u;
y13 = .5672830237u;
draw fullcircle scaled D shifted z13 withcolor RoyalBlue withpen pencircle scaled 1.0;

x14 = .6596684004u;
y14 = .6004999699u;
draw fullcircle scaled D shifted z14 withcolor RoyalBlue withpen pencircle scaled 1.0;

x15 = .7196950410u;
y15 = .6369005855u;
draw fullcircle scaled D shifted z15 withcolor RoyalBlue withpen pencircle scaled 1.0;

x16 = .5568559730u;
y16 = .5770353555u;
draw fullcircle scaled D shifted z16 withcolor RoyalBlue withpen pencircle scaled 1.0;

x17 = .5999921296u;
y17 = .6322418872u;
draw fullcircle scaled D shifted z17 withcolor RoyalBlue withpen pencircle scaled 1.0;

x18 = .6604720565u;
y18 = .6621836518u;
draw fullcircle scaled D shifted z18 withcolor RoyalBlue withpen pencircle scaled 1.0;

x19 = .5407001816u;
y19 = .6383317133u;
draw fullcircle scaled D shifted z19 withcolor RoyalBlue withpen pencircle scaled 1.0;

x20 = .5688105745u;
y20 = .6964183853u;
draw fullcircle scaled D shifted z20 withcolor RoyalBlue withpen pencircle scaled 1.0;

x21 = .6153691022u;
y21 = .7000647377u;
draw fullcircle scaled D shifted z21 withcolor RoyalBlue withpen pencircle scaled 1.0;

x22 = .5010691918u;
y22 = .7281734295u;
draw fullcircle scaled D shifted z22 withcolor RoyalBlue withpen pencircle scaled 1.0;

x23 = .4608619437u;
y23 = .7824999999u;
draw fullcircle scaled D shifted z23 withcolor RoyalBlue withpen pencircle scaled 1.0;

x24 = .5374484676u;
y24 = .8262497379u;
draw fullcircle scaled D shifted z24 withcolor RoyalBlue withpen pencircle scaled 1.0;

x25 = .4139944147u;
y25 = .8011871977u;
draw fullcircle scaled D shifted z25 withcolor RoyalBlue withpen pencircle scaled 1.0;

x26 = .4260522514u;
y26 = .8672275259u;
draw fullcircle scaled D shifted z26 withcolor RoyalBlue withpen pencircle scaled 1.0;

x27 = .4928968169u;
y27 = .8836642620u;
draw fullcircle scaled D shifted z27 withcolor RoyalBlue withpen pencircle scaled 1.0;

x28 = .5420500412u;
y28 = .8968005333u;
draw fullcircle scaled D shifted z28 withcolor RoyalBlue withpen pencircle scaled 1.0;

x29 = .5615879005u;
y29 = .9585227374u;
draw fullcircle scaled D shifted z29 withcolor RoyalBlue withpen pencircle scaled 1.0;

x30 = .6172117132u;
y30 = .9799937933u;
draw fullcircle scaled D shifted z30 withcolor RoyalBlue withpen pencircle scaled 1.0;

x31 = .3661072462u;
y31 = .8517117216u;
draw fullcircle scaled D shifted z31 withcolor RoyalBlue withpen pencircle scaled 1.0;

x32 = .4455719853u;
y32 = .9289563292u;
draw fullcircle scaled D shifted z32 withcolor RoyalBlue withpen pencircle scaled 1.0;

x33 = .5002120193u;
y33 = .9550258571u;
draw fullcircle scaled D shifted z33 withcolor RoyalBlue withpen pencircle scaled 1.0;

x34 = .5255661946u;
y34 = 1.0168977635u;
draw fullcircle scaled D shifted z34 withcolor RoyalBlue withpen pencircle scaled 1.0;

x35 = .5849102882u;
y35 = 1.0242725299u;
draw fullcircle scaled D shifted z35 withcolor RoyalBlue withpen pencircle scaled 1.0;

x36 = .6529612200u;
y36 = 1.0400540135u;
draw fullcircle scaled D shifted z36 withcolor RoyalBlue withpen pencircle scaled 1.0;

x37 = .3221193845u;
y37 = .9070774036u;
draw fullcircle scaled D shifted z37 withcolor RoyalBlue withpen pencircle scaled 1.0;

x38 = .3803892025u;
y38 = .9180417992u;
draw fullcircle scaled D shifted z38 withcolor RoyalBlue withpen pencircle scaled 1.0;

x39 = .3980602034u;
y39 = .9773427301u;
draw fullcircle scaled D shifted z39 withcolor RoyalBlue withpen pencircle scaled 1.0;

x40 = .4572953307u;
y40 = .9966515659u;
draw fullcircle scaled D shifted z40 withcolor RoyalBlue withpen pencircle scaled 1.0;

x41 = .4804555594u;
y41 = 1.0526291751u;
draw fullcircle scaled D shifted z41 withcolor RoyalBlue withpen pencircle scaled 1.0;

x42 = .6014974179u;
y42 = 1.0833863447u;
draw fullcircle scaled D shifted z42 withcolor RoyalBlue withpen pencircle scaled 1.0;

x43 = .3301783185u;
y43 = .9630220499u;
draw fullcircle scaled D shifted z43 withcolor RoyalBlue withpen pencircle scaled 1.0;

x44 = .3505692260u;
y44 = 1.0231789213u;
draw fullcircle scaled D shifted z44 withcolor RoyalBlue withpen pencircle scaled 1.0;

x45 = .4105703335u;
y45 = 1.0452503476u;
draw fullcircle scaled D shifted z45 withcolor RoyalBlue withpen pencircle scaled 1.0;

x46 = .4279770410u;
y46 = 1.0991399966u;
draw fullcircle scaled D shifted z46 withcolor RoyalBlue withpen pencircle scaled 1.0;

x47 = .4867275029u;
y47 = 1.1177078606u;
draw fullcircle scaled D shifted z47 withcolor RoyalBlue withpen pencircle scaled 1.0;

x48 = .5375535653u;
y48 = 1.1583213228u;
draw fullcircle scaled D shifted z48 withcolor RoyalBlue withpen pencircle scaled 1.0;

x49 = .3099710325u;
y49 = 1.0549867056u;
draw fullcircle scaled D shifted z49 withcolor RoyalBlue withpen pencircle scaled 1.0;

x50 = .3646722590u;
y50 = 1.0915197181u;
draw fullcircle scaled D shifted z50 withcolor RoyalBlue withpen pencircle scaled 1.0;

x51 = .4623000928u;
y51 = 1.1602107096u;
draw fullcircle scaled D shifted z51 withcolor RoyalBlue withpen pencircle scaled 1.0;

x52 = .5589057283u;
y52 = 1.1937674399u;
draw fullcircle scaled D shifted z52 withcolor RoyalBlue withpen pencircle scaled 1.0;

x53 = .3021129054u;
y53 = 1.1149840488u;
draw fullcircle scaled D shifted z53 withcolor RoyalBlue withpen pencircle scaled 1.0;

x54 = .3435696844u;
y54 = 1.1591009039u;
draw fullcircle scaled D shifted z54 withcolor RoyalBlue withpen pencircle scaled 1.0;

x55 = .3973133608u;
y55 = 1.1546950559u;
draw fullcircle scaled D shifted z55 withcolor RoyalBlue withpen pencircle scaled 1.0;

x56 = .1657595975u;
y56 = 1.1099734207u;
draw fullcircle scaled D shifted z56 withcolor RoyalBlue withpen pencircle scaled 1.0;

x57 = .2347964888u;
y57 = 1.1211528506u;
draw fullcircle scaled D shifted z57 withcolor RoyalBlue withpen pencircle scaled 1.0;

x58 = .2836255647u;
y58 = 1.1720169523u;
draw fullcircle scaled D shifted z58 withcolor RoyalBlue withpen pencircle scaled 1.0;

x59 = .1679437295u;
y59 = 1.1552482269u;
draw fullcircle scaled D shifted z59 withcolor RoyalBlue withpen pencircle scaled 1.0;

x60 = .2180866253u;
y60 = 1.1806641696u;
draw fullcircle scaled D shifted z60 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.4608619437u,        0.7824999999u)--(        0.4950000000u,        0.8000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4608619437u,        0.7824999999u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4608619437u,        0.7824999999u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4950000000u,        0.8000000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4950000000u,        0.8000000000u) withcolor Yellow;
draw (        0.5374484676u,        0.8262497379u)--(        0.5500000000u,        0.7900000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5374484676u,        0.8262497379u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5374484676u,        0.8262497379u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5500000000u,        0.7900000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5500000000u,        0.7900000000u) withcolor Yellow;
draw (        0.5374484676u,        0.8262497379u)--(        0.5248673512u,        0.8624797530u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4928968169u,        0.8836642620u)--(        0.5248673512u,        0.8624797530u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5374484676u,        0.8262497379u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5374484676u,        0.8262497379u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4928968169u,        0.8836642620u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4928968169u,        0.8836642620u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5248673512u,        0.8624797530u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5248673512u,        0.8624797530u) withcolor MidnightBlue;
draw (        0.5374484676u,        0.8262497379u)--(        0.5249062782u,        0.8624932473u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5420500412u,        0.8968005333u)--(        0.5249062782u,        0.8624932473u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5374484676u,        0.8262497379u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5374484676u,        0.8262497379u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5420500412u,        0.8968005333u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5420500412u,        0.8968005333u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5249062782u,        0.8624932473u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5249062782u,        0.8624932473u) withcolor MidnightBlue;
draw (        0.4139944147u,        0.8011871977u)--(        0.4017719081u,        0.8375397604u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4260522514u,        0.8672275259u)--(        0.4017719081u,        0.8375397604u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4139944147u,        0.8011871977u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4139944147u,        0.8011871977u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4260522514u,        0.8672275259u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4260522514u,        0.8672275259u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4017719081u,        0.8375397604u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4017719081u,        0.8375397604u) withcolor MidnightBlue;
draw (        0.4139944147u,        0.8011871977u)--(        0.3783607496u,        0.8153695953u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3661072462u,        0.8517117216u)--(        0.3783607496u,        0.8153695953u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4139944147u,        0.8011871977u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4139944147u,        0.8011871977u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3661072462u,        0.8517117216u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3661072462u,        0.8517117216u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3783607496u,        0.8153695953u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3783607496u,        0.8153695953u) withcolor Green;
draw (        0.4139944147u,        0.8011871977u)--(        0.4017409112u,        0.8375293239u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3661072462u,        0.8517117216u)--(        0.4017409112u,        0.8375293239u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4139944147u,        0.8011871977u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4139944147u,        0.8011871977u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3661072462u,        0.8517117216u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3661072462u,        0.8517117216u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4017409112u,        0.8375293239u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4017409112u,        0.8375293239u) withcolor MidnightBlue;
draw (        0.4260522514u,        0.8672275259u)--(        0.4554341920u,        0.8918770706u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4928968169u,        0.8836642620u)--(        0.4554341920u,        0.8918770706u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4260522514u,        0.8672275259u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4260522514u,        0.8672275259u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4928968169u,        0.8836642620u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4928968169u,        0.8836642620u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4554341920u,        0.8918770706u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4554341920u,        0.8918770706u) withcolor MidnightBlue;
draw (        0.4260522514u,        0.8672275259u)--(        0.4635148763u,        0.8590147173u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4928968169u,        0.8836642620u)--(        0.4635148763u,        0.8590147173u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4260522514u,        0.8672275259u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4260522514u,        0.8672275259u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4928968169u,        0.8836642620u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4928968169u,        0.8836642620u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4635148763u,        0.8590147173u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4635148763u,        0.8590147173u) withcolor Green;
draw (        0.4260522514u,        0.8672275259u)--(        0.4017516515u,        0.8375563390u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3661072462u,        0.8517117216u)--(        0.4017516515u,        0.8375563390u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4260522514u,        0.8672275259u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4260522514u,        0.8672275259u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3661072462u,        0.8517117216u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3661072462u,        0.8517117216u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4017516515u,        0.8375563390u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4017516515u,        0.8375563390u) withcolor MidnightBlue;
draw (        0.4260522514u,        0.8672275259u)--(        0.4162014842u,        0.9042931551u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4455719853u,        0.9289563292u)--(        0.4162014842u,        0.9042931551u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4260522514u,        0.8672275259u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4260522514u,        0.8672275259u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4455719853u,        0.9289563292u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4455719853u,        0.9289563292u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4162014842u,        0.9042931551u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4162014842u,        0.9042931551u) withcolor MidnightBlue;
draw (        0.4260522514u,        0.8672275259u)--(        0.4554227525u,        0.8918906999u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4455719853u,        0.9289563292u)--(        0.4554227525u,        0.8918906999u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4260522514u,        0.8672275259u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4260522514u,        0.8672275259u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4455719853u,        0.9289563292u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4455719853u,        0.9289563292u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4554227525u,        0.8918906999u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4554227525u,        0.8918906999u) withcolor MidnightBlue;
draw (        0.4260522514u,        0.8672275259u)--(        0.4161914340u,        0.9042904827u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3803892025u,        0.9180417992u)--(        0.4161914340u,        0.9042904827u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4260522514u,        0.8672275259u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4260522514u,        0.8672275259u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3803892025u,        0.9180417992u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3803892025u,        0.9180417992u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4161914340u,        0.9042904827u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4161914340u,        0.9042904827u) withcolor MidnightBlue;
draw (        0.4928968169u,        0.8836642620u)--(        0.5100631120u,        0.9179602793u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5420500412u,        0.8968005333u)--(        0.5100631120u,        0.9179602793u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4928968169u,        0.8836642620u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4928968169u,        0.8836642620u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5420500412u,        0.8968005333u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5420500412u,        0.8968005333u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5100631120u,        0.9179602793u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5100631120u,        0.9179602793u) withcolor MidnightBlue;
draw (        0.4928968169u,        0.8836642620u)--(        0.5248837461u,        0.8625045161u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5420500412u,        0.8968005333u)--(        0.5248837461u,        0.8625045161u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4928968169u,        0.8836642620u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4928968169u,        0.8836642620u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5420500412u,        0.8968005333u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5420500412u,        0.8968005333u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5248837461u,        0.8625045161u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5248837461u,        0.8625045161u) withcolor MidnightBlue;
draw (        0.4928968169u,        0.8836642620u)--(        0.4554380779u,        0.8918947763u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4455719853u,        0.9289563292u)--(        0.4554380779u,        0.8918947763u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4928968169u,        0.8836642620u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4928968169u,        0.8836642620u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4455719853u,        0.9289563292u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4455719853u,        0.9289563292u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4554380779u,        0.8918947763u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4554380779u,        0.8918947763u) withcolor MidnightBlue;
draw (        0.4928968169u,        0.8836642620u)--(        0.4830307242u,        0.9207258150u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4455719853u,        0.9289563292u)--(        0.4830307242u,        0.9207258150u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4928968169u,        0.8836642620u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4928968169u,        0.8836642620u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4455719853u,        0.9289563292u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4455719853u,        0.9289563292u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4830307242u,        0.9207258150u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4830307242u,        0.9207258150u) withcolor MidnightBlue;
draw (        0.4928968169u,        0.8836642620u)--(        0.4830457699u,        0.9207298170u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5002120193u,        0.9550258571u)--(        0.4830457699u,        0.9207298170u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4928968169u,        0.8836642620u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4928968169u,        0.8836642620u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5002120193u,        0.9550258571u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5002120193u,        0.9550258571u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4830457699u,        0.9207298170u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4830457699u,        0.9207298170u) withcolor MidnightBlue;
draw (        0.4928968169u,        0.8836642620u)--(        0.5100630663u,        0.9179603022u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5002120193u,        0.9550258571u)--(        0.5100630663u,        0.9179603022u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4928968169u,        0.8836642620u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4928968169u,        0.8836642620u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5002120193u,        0.9550258571u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5002120193u,        0.9550258571u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5100630663u,        0.9179603022u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5100630663u,        0.9179603022u) withcolor MidnightBlue;
draw (        0.5420500412u,        0.8968005333u)--(        0.5322095650u,        0.9338688961u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5615879005u,        0.9585227374u)--(        0.5322095650u,        0.9338688961u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5420500412u,        0.8968005333u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5420500412u,        0.8968005333u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5615879005u,        0.9585227374u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5615879005u,        0.9585227374u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5322095650u,        0.9338688961u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5322095650u,        0.9338688961u) withcolor MidnightBlue;
draw (        0.5420500412u,        0.8968005333u)--(        0.5100631401u,        0.9179603218u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5002120193u,        0.9550258571u)--(        0.5100631401u,        0.9179603218u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5420500412u,        0.8968005333u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5420500412u,        0.8968005333u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5002120193u,        0.9550258571u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5002120193u,        0.9550258571u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5100631401u,        0.9179603218u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5100631401u,        0.9179603218u) withcolor MidnightBlue;
draw (        0.5420500412u,        0.8968005333u)--(        0.5321989204u,        0.9338660686u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5002120193u,        0.9550258571u)--(        0.5321989204u,        0.9338660686u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5420500412u,        0.8968005333u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5420500412u,        0.8968005333u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5002120193u,        0.9550258571u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5002120193u,        0.9550258571u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5321989204u,        0.9338660686u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5321989204u,        0.9338660686u) withcolor MidnightBlue;
draw (        0.5615879005u,        0.9585227374u)--(        0.5322046580u,        0.9338747447u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5002120193u,        0.9550258571u)--(        0.5322046580u,        0.9338747447u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5615879005u,        0.9585227374u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5615879005u,        0.9585227374u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5002120193u,        0.9550258571u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5002120193u,        0.9550258571u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5322046580u,        0.9338747447u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5322046580u,        0.9338747447u) withcolor MidnightBlue;
draw (        0.5615879005u,        0.9585227374u)--(        0.5581837337u,        0.9967236554u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5255661946u,        1.0168977635u)--(        0.5581837337u,        0.9967236554u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5615879005u,        0.9585227374u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5615879005u,        0.9585227374u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5255661946u,        1.0168977635u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5255661946u,        1.0168977635u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5581837337u,        0.9967236554u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5581837337u,        0.9967236554u) withcolor MidnightBlue;
draw (        0.5615879005u,        0.9585227374u)--(        0.5582240595u,        0.9967272275u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5849102882u,        1.0242725299u)--(        0.5582240595u,        0.9967272275u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5615879005u,        0.9585227374u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5615879005u,        0.9585227374u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5849102882u,        1.0242725299u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5849102882u,        1.0242725299u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5582240595u,        0.9967272275u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5582240595u,        0.9967272275u) withcolor MidnightBlue;
draw (        0.3661072462u,        0.8517117216u)--(        0.3550000000u,        0.8150000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3661072462u,        0.8517117216u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3661072462u,        0.8517117216u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3550000000u,        0.8150000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3550000000u,        0.8150000000u) withcolor Yellow;
draw (        0.3661072462u,        0.8517117216u)--(        0.3557483569u,        0.8886385669u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3221193845u,        0.9070774036u)--(        0.3557483569u,        0.8886385669u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3661072462u,        0.8517117216u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3661072462u,        0.8517117216u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3221193845u,        0.9070774036u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3221193845u,        0.9070774036u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3557483569u,        0.8886385669u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3557483569u,        0.8886385669u) withcolor MidnightBlue;
draw (        0.3661072462u,        0.8517117216u)--(        0.3557608710u,        0.8886420751u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3803892025u,        0.9180417992u)--(        0.3557608710u,        0.8886420751u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3661072462u,        0.8517117216u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3661072462u,        0.8517117216u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3803892025u,        0.9180417992u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3803892025u,        0.9180417992u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3557608710u,        0.8886420751u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3557608710u,        0.8886420751u) withcolor MidnightBlue;
draw (        0.4455719853u,        0.9289563292u)--(        0.4830330064u,        0.9207362085u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5002120193u,        0.9550258571u)--(        0.4830330064u,        0.9207362085u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4455719853u,        0.9289563292u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4455719853u,        0.9289563292u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5002120193u,        0.9550258571u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5002120193u,        0.9550258571u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4830330064u,        0.9207362085u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4830330064u,        0.9207362085u) withcolor MidnightBlue;
draw (        0.4455719853u,        0.9289563292u)--(        0.4161952919u,        0.9043005312u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3803892025u,        0.9180417992u)--(        0.4161952919u,        0.9043005312u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4455719853u,        0.9289563292u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4455719853u,        0.9289563292u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3803892025u,        0.9180417992u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3803892025u,        0.9180417992u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4161952919u,        0.9043005312u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4161952919u,        0.9043005312u) withcolor MidnightBlue;
draw (        0.4455719853u,        0.9289563292u)--(        0.4346048471u,        0.9657071165u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3980602034u,        0.9773427301u)--(        0.4346048471u,        0.9657071165u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4455719853u,        0.9289563292u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4455719853u,        0.9289563292u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3980602034u,        0.9773427301u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3980602034u,        0.9773427301u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4346048471u,        0.9657071165u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4346048471u,        0.9657071165u) withcolor MidnightBlue;
draw (        0.4455719853u,        0.9289563292u)--(        0.4346286473u,        0.9657142105u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4572953307u,        0.9966515659u)--(        0.4346286473u,        0.9657142105u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4455719853u,        0.9289563292u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4455719853u,        0.9289563292u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4572953307u,        0.9966515659u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4572953307u,        0.9966515659u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4346286473u,        0.9657142105u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4346286473u,        0.9657142105u) withcolor MidnightBlue;
draw (        0.5002120193u,        0.9550258571u)--(        0.4955000307u,        0.9930875918u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5255661946u,        1.0168977635u)--(        0.4955000307u,        0.9930875918u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5002120193u,        0.9550258571u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5002120193u,        0.9550258571u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5255661946u,        1.0168977635u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5255661946u,        1.0168977635u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4955000307u,        0.9930875918u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4955000307u,        0.9930875918u) withcolor MidnightBlue;
draw (        0.5002120193u,        0.9550258571u)--(        0.4954814557u,        0.9930852877u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4572953307u,        0.9966515659u)--(        0.4954814557u,        0.9930852877u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5002120193u,        0.9550258571u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5002120193u,        0.9550258571u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4572953307u,        0.9966515659u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4572953307u,        0.9966515659u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4954814557u,        0.9930852877u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4954814557u,        0.9930852877u) withcolor MidnightBlue;
draw (        0.5255661946u,        1.0168977635u)--(        0.5582002233u,        0.9967503403u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5849102882u,        1.0242725299u)--(        0.5582002233u,        0.9967503403u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5255661946u,        1.0168977635u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5255661946u,        1.0168977635u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5849102882u,        1.0242725299u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5849102882u,        1.0242725299u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5582002233u,        0.9967503403u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5582002233u,        0.9967503403u) withcolor MidnightBlue;
draw (        0.5255661946u,        1.0168977635u)--(        0.4954836026u,        0.9931083510u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4572953307u,        0.9966515659u)--(        0.4954836026u,        0.9931083510u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5255661946u,        1.0168977635u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5255661946u,        1.0168977635u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4572953307u,        0.9966515659u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4572953307u,        0.9966515659u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4954836026u,        0.9931083510u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4954836026u,        0.9931083510u) withcolor MidnightBlue;
draw (        0.5849102882u,        1.0242725299u)--(        0.6153573369u,        1.0475936730u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6529612200u,        1.0400540135u)--(        0.6153573369u,        1.0475936730u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5849102882u,        1.0242725299u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5849102882u,        1.0242725299u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6529612200u,        1.0400540135u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6529612200u,        1.0400540135u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6153573369u,        1.0475936730u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6153573369u,        1.0475936730u) withcolor MidnightBlue;
draw (        0.5849102882u,        1.0242725299u)--(        0.5710694874u,        1.0600402629u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6014974179u,        1.0833863447u)--(        0.5710694874u,        1.0600402629u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5849102882u,        1.0242725299u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5849102882u,        1.0242725299u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6014974179u,        1.0833863447u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6014974179u,        1.0833863447u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5710694874u,        1.0600402629u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5710694874u,        1.0600402629u) withcolor Green;
draw (        0.5849102882u,        1.0242725299u)--(        0.6153382187u,        1.0476186117u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6014974179u,        1.0833863447u)--(        0.6153382187u,        1.0476186117u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5849102882u,        1.0242725299u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5849102882u,        1.0242725299u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6014974179u,        1.0833863447u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6014974179u,        1.0833863447u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6153382187u,        1.0476186117u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6153382187u,        1.0476186117u) withcolor MidnightBlue;
draw (        0.6529612200u,        1.0400540135u)--(        0.6900000000u,        1.0300000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6529612200u,        1.0400540135u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6529612200u,        1.0400540135u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6900000000u,        1.0300000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6900000000u,        1.0300000000u) withcolor Yellow;
draw (        0.6529612200u,        1.0400540135u)--(        0.6153643867u,        1.0476287488u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6014974179u,        1.0833863447u)--(        0.6153643867u,        1.0476287488u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6529612200u,        1.0400540135u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6529612200u,        1.0400540135u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6014974179u,        1.0833863447u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6014974179u,        1.0833863447u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6153643867u,        1.0476287488u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6153643867u,        1.0476287488u) withcolor MidnightBlue;
draw (        0.6529612200u,        1.0400540135u)--(        0.6390942512u,        1.0758116094u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6014974179u,        1.0833863447u)--(        0.6390942512u,        1.0758116094u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6529612200u,        1.0400540135u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6529612200u,        1.0400540135u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6014974179u,        1.0833863447u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6014974179u,        1.0833863447u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6390942512u,        1.0758116094u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6390942512u,        1.0758116094u) withcolor Green;
draw (        0.3221193845u,        0.9070774036u)--(        0.3557536132u,        0.8886481565u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3803892025u,        0.9180417992u)--(        0.3557536132u,        0.8886481565u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3221193845u,        0.9070774036u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3221193845u,        0.9070774036u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3803892025u,        0.9180417992u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3803892025u,        0.9180417992u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3557536132u,        0.8886481565u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3557536132u,        0.8886481565u) withcolor MidnightBlue;
draw (        0.3803892025u,        0.9180417992u)--(        0.3675040033u,        0.9541647853u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3980602034u,        0.9773427301u)--(        0.3675040033u,        0.9541647853u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3803892025u,        0.9180417992u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3803892025u,        0.9180417992u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3980602034u,        0.9773427301u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3980602034u,        0.9773427301u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3675040033u,        0.9541647853u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3675040033u,        0.9541647853u) withcolor MidnightBlue;
draw (        0.3803892025u,        0.9180417992u)--(        0.3674928876u,        0.9541608184u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3301783185u,        0.9630220499u)--(        0.3674928876u,        0.9541608184u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3803892025u,        0.9180417992u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3803892025u,        0.9180417992u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3301783185u,        0.9630220499u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3301783185u,        0.9630220499u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3674928876u,        0.9541608184u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3674928876u,        0.9541608184u) withcolor MidnightBlue;
draw (        0.3980602034u,        0.9773427301u)--(        0.4207443565u,        1.0082672783u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4572953307u,        0.9966515659u)--(        0.4207443565u,        1.0082672783u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3980602034u,        0.9773427301u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3980602034u,        0.9773427301u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4572953307u,        0.9966515659u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4572953307u,        0.9966515659u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4207443565u,        1.0082672783u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4207443565u,        1.0082672783u) withcolor MidnightBlue;
draw (        0.3980602034u,        0.9773427301u)--(        0.4346111776u,        0.9657270176u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4572953307u,        0.9966515659u)--(        0.4346111776u,        0.9657270176u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3980602034u,        0.9773427301u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3980602034u,        0.9773427301u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4572953307u,        0.9966515659u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4572953307u,        0.9966515659u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4346111776u,        0.9657270176u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4346111776u,        0.9657270176u) withcolor MidnightBlue;
draw (        0.3980602034u,        0.9773427301u)--(        0.3674962538u,        0.9541750053u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3301783185u,        0.9630220499u)--(        0.3674962538u,        0.9541750053u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3980602034u,        0.9773427301u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3980602034u,        0.9773427301u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3301783185u,        0.9630220499u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3301783185u,        0.9630220499u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3674962538u,        0.9541750053u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3674962538u,        0.9541750053u) withcolor MidnightBlue;
draw (        0.3980602034u,        0.9773427301u)--(        0.3607422682u,        0.9861897746u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3301783185u,        0.9630220499u)--(        0.3607422682u,        0.9861897746u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3980602034u,        0.9773427301u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3980602034u,        0.9773427301u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3301783185u,        0.9630220499u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3301783185u,        0.9630220499u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3607422682u,        0.9861897746u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3607422682u,        0.9861897746u) withcolor MidnightBlue;
draw (        0.3980602034u,        0.9773427301u)--(        0.3607449763u,        0.9862011903u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3505692260u,        1.0231789213u)--(        0.3607449763u,        0.9862011903u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3980602034u,        0.9773427301u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3980602034u,        0.9773427301u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3505692260u,        1.0231789213u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3505692260u,        1.0231789213u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3607449763u,        0.9862011903u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3607449763u,        0.9862011903u) withcolor MidnightBlue;
draw (        0.3980602034u,        0.9773427301u)--(        0.3878844531u,        1.0143204611u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3505692260u,        1.0231789213u)--(        0.3878844531u,        1.0143204611u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3980602034u,        0.9773427301u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3980602034u,        0.9773427301u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3505692260u,        1.0231789213u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3505692260u,        1.0231789213u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3878844531u,        1.0143204611u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3878844531u,        1.0143204611u) withcolor MidnightBlue;
draw (        0.3980602034u,        0.9773427301u)--(        0.3878910024u,        1.0143222628u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4105703335u,        1.0452503476u)--(        0.3878910024u,        1.0143222628u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3980602034u,        0.9773427301u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3980602034u,        0.9773427301u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4105703335u,        1.0452503476u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4105703335u,        1.0452503476u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3878910024u,        1.0143222628u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3878910024u,        1.0143222628u) withcolor MidnightBlue;
draw (        0.3980602034u,        0.9773427301u)--(        0.4207395344u,        1.0082708149u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4105703335u,        1.0452503476u)--(        0.4207395344u,        1.0082708149u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3980602034u,        0.9773427301u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3980602034u,        0.9773427301u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4105703335u,        1.0452503476u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4105703335u,        1.0452503476u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4207395344u,        1.0082708149u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4207395344u,        1.0082708149u) withcolor MidnightBlue;
draw (        0.4572953307u,        0.9966515659u)--(        0.4471375360u,        1.0336342334u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4804555594u,        1.0526291751u)--(        0.4471375360u,        1.0336342334u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4572953307u,        0.9966515659u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4572953307u,        0.9966515659u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4804555594u,        1.0526291751u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4804555594u,        1.0526291751u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4471375360u,        1.0336342334u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4471375360u,        1.0336342334u) withcolor MidnightBlue;
draw (        0.4572953307u,        0.9966515659u)--(        0.4207460505u,        1.0082726074u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4105703335u,        1.0452503476u)--(        0.4207460505u,        1.0082726074u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4572953307u,        0.9966515659u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4572953307u,        0.9966515659u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4105703335u,        1.0452503476u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4105703335u,        1.0452503476u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4207460505u,        1.0082726074u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4207460505u,        1.0082726074u) withcolor MidnightBlue;
draw (        0.4572953307u,        0.9966515659u)--(        0.4471196137u,        1.0336293061u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4105703335u,        1.0452503476u)--(        0.4471196137u,        1.0336293061u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4572953307u,        0.9966515659u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4572953307u,        0.9966515659u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4105703335u,        1.0452503476u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4105703335u,        1.0452503476u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4471196137u,        1.0336293061u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4471196137u,        1.0336293061u) withcolor MidnightBlue;
draw (        0.4804555594u,        1.0526291751u)--(        0.4471270333u,        1.0336526675u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4105703335u,        1.0452503476u)--(        0.4471270333u,        1.0336526675u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4804555594u,        1.0526291751u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4804555594u,        1.0526291751u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4105703335u,        1.0452503476u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4105703335u,        1.0452503476u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4471270333u,        1.0336526675u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4471270333u,        1.0336526675u) withcolor MidnightBlue;
draw (        0.4804555594u,        1.0526291751u)--(        0.4438988595u,        1.0642268552u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4105703335u,        1.0452503476u)--(        0.4438988595u,        1.0642268552u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4804555594u,        1.0526291751u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4804555594u,        1.0526291751u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4105703335u,        1.0452503476u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4105703335u,        1.0452503476u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4438988595u,        1.0642268552u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4438988595u,        1.0642268552u) withcolor MidnightBlue;
draw (        0.4804555594u,        1.0526291751u)--(        0.4439070196u,        1.0642525452u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4279770410u,        1.0991399966u)--(        0.4439070196u,        1.0642525452u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4804555594u,        1.0526291751u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4804555594u,        1.0526291751u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4279770410u,        1.0991399966u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4279770410u,        1.0991399966u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4439070196u,        1.0642525452u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4439070196u,        1.0642525452u) withcolor MidnightBlue;
draw (        0.6014974179u,        1.0833863447u)--(        0.5900000000u,        1.1200000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6014974179u,        1.0833863447u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6014974179u,        1.0833863447u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5900000000u,        1.1200000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5900000000u,        1.1200000000u) withcolor Yellow;
draw (        0.3301783185u,        0.9630220499u)--(        0.3200000000u,        1.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3301783185u,        0.9630220499u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3301783185u,        0.9630220499u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3200000000u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3200000000u,        1.0000000000u) withcolor Yellow;
draw (        0.3301783185u,        0.9630220499u)--(        0.3607355734u,        0.9861986040u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3505692260u,        1.0231789213u)--(        0.3607355734u,        0.9861986040u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3301783185u,        0.9630220499u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3301783185u,        0.9630220499u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3505692260u,        1.0231789213u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3505692260u,        1.0231789213u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3607355734u,        0.9861986040u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3607355734u,        0.9861986040u) withcolor MidnightBlue;
draw (        0.3505692260u,        1.0231789213u)--(        0.3200000000u,        1.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3505692260u,        1.0231789213u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3505692260u,        1.0231789213u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3200000000u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3200000000u,        1.0000000000u) withcolor Yellow;
draw (        0.3505692260u,        1.0231789213u)--(        0.3732537713u,        1.0541031819u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4105703335u,        1.0452503476u)--(        0.3732537713u,        1.0541031819u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3505692260u,        1.0231789213u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3505692260u,        1.0231789213u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4105703335u,        1.0452503476u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4105703335u,        1.0452503476u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3732537713u,        1.0541031819u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3732537713u,        1.0541031819u) withcolor MidnightBlue;
draw (        0.3505692260u,        1.0231789213u)--(        0.3878857882u,        1.0143260870u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4105703335u,        1.0452503476u)--(        0.3878857882u,        1.0143260870u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3505692260u,        1.0231789213u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3505692260u,        1.0231789213u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4105703335u,        1.0452503476u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4105703335u,        1.0452503476u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3878857882u,        1.0143260870u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3878857882u,        1.0143260870u) withcolor MidnightBlue;
draw (        0.3505692260u,        1.0231789213u)--(        0.3732156629u,        1.0541311003u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3646722590u,        1.0915197181u)--(        0.3732156629u,        1.0541311003u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3505692260u,        1.0231789213u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3505692260u,        1.0231789213u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3646722590u,        1.0915197181u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3646722590u,        1.0915197181u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3732156629u,        1.0541311003u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3732156629u,        1.0541311003u) withcolor MidnightBlue;
draw (        0.4105703335u,        1.0452503476u)--(        0.4438889268u,        1.0642442894u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4279770410u,        1.0991399966u)--(        0.4438889268u,        1.0642442894u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4105703335u,        1.0452503476u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4105703335u,        1.0452503476u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4279770410u,        1.0991399966u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4279770410u,        1.0991399966u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4438889268u,        1.0642442894u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4438889268u,        1.0642442894u) withcolor MidnightBlue;
draw (        0.4105703335u,        1.0452503476u)--(        0.3732629888u,        1.0541419459u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3646722590u,        1.0915197181u)--(        0.3732629888u,        1.0541419459u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4105703335u,        1.0452503476u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4105703335u,        1.0452503476u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3646722590u,        1.0915197181u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3646722590u,        1.0915197181u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3732629888u,        1.0541419459u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3732629888u,        1.0541419459u) withcolor MidnightBlue;
draw (        0.4279770410u,        1.0991399966u)--(        0.3937767583u,        1.1164962391u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3646722590u,        1.0915197181u)--(        0.3937767583u,        1.1164962391u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4279770410u,        1.0991399966u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4279770410u,        1.0991399966u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3646722590u,        1.0915197181u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3646722590u,        1.0915197181u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3937767583u,        1.1164962391u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3937767583u,        1.1164962391u) withcolor MidnightBlue;
draw (        0.4279770410u,        1.0991399966u)--(        0.4315221052u,        1.1373280969u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4623000928u,        1.1602107096u)--(        0.4315221052u,        1.1373280969u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4279770410u,        1.0991399966u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4279770410u,        1.0991399966u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4623000928u,        1.1602107096u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4623000928u,        1.1602107096u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4315221052u,        1.1373280969u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4315221052u,        1.1373280969u) withcolor MidnightBlue;
draw (        0.4279770410u,        1.0991399966u)--(        0.3937815741u,        1.1165057253u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3973133608u,        1.1546950559u)--(        0.3937815741u,        1.1165057253u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4279770410u,        1.0991399966u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4279770410u,        1.0991399966u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3973133608u,        1.1546950559u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3973133608u,        1.1546950559u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3937815741u,        1.1165057253u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3937815741u,        1.1165057253u) withcolor MidnightBlue;
draw (        0.4279770410u,        1.0991399966u)--(        0.4315088277u,        1.1373293272u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3973133608u,        1.1546950559u)--(        0.4315088277u,        1.1373293272u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4279770410u,        1.0991399966u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4279770410u,        1.0991399966u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3973133608u,        1.1546950559u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3973133608u,        1.1546950559u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4315088277u,        1.1373293272u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4315088277u,        1.1373293272u) withcolor MidnightBlue;
draw (        0.5375535653u,        1.1583213228u)--(        0.5750000000u,        1.1500000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5375535653u,        1.1583213228u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5375535653u,        1.1583213228u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5750000000u,        1.1500000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5750000000u,        1.1500000000u) withcolor Yellow;
draw (        0.3646722590u,        1.0915197181u)--(        0.3400078624u,        1.1208891925u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3021129054u,        1.1149840488u)--(        0.3400078624u,        1.1208891925u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3646722590u,        1.0915197181u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3646722590u,        1.0915197181u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3021129054u,        1.1149840488u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3021129054u,        1.1149840488u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3400078624u,        1.1208891925u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3400078624u,        1.1208891925u) withcolor MidnightBlue;
draw (        0.3646722590u,        1.0915197181u)--(        0.3400348713u,        1.1209118533u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3435696844u,        1.1591009039u)--(        0.3400348713u,        1.1209118533u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3646722590u,        1.0915197181u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3646722590u,        1.0915197181u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3435696844u,        1.1591009039u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3435696844u,        1.1591009039u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3400348713u,        1.1209118533u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3400348713u,        1.1209118533u) withcolor MidnightBlue;
draw (        0.3646722590u,        1.0915197181u)--(        0.3682070721u,        1.1297087686u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3435696844u,        1.1591009039u)--(        0.3682070721u,        1.1297087686u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3646722590u,        1.0915197181u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3646722590u,        1.0915197181u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3435696844u,        1.1591009039u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3435696844u,        1.1591009039u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3682070721u,        1.1297087686u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3682070721u,        1.1297087686u) withcolor MidnightBlue;
draw (        0.3646722590u,        1.0915197181u)--(        0.3682181253u,        1.1297077439u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3973133608u,        1.1546950559u)--(        0.3682181253u,        1.1297077439u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3646722590u,        1.0915197181u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3646722590u,        1.0915197181u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3973133608u,        1.1546950559u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3973133608u,        1.1546950559u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3682181253u,        1.1297077439u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3682181253u,        1.1297077439u) withcolor MidnightBlue;
draw (        0.3646722590u,        1.0915197181u)--(        0.3937674945u,        1.1165070300u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3973133608u,        1.1546950559u)--(        0.3937674945u,        1.1165070300u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3646722590u,        1.0915197181u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3646722590u,        1.0915197181u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3973133608u,        1.1546950559u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3973133608u,        1.1546950559u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3937674945u,        1.1165070300u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3937674945u,        1.1165070300u) withcolor MidnightBlue;
draw (        0.4623000928u,        1.1602107096u)--(        0.4315138453u,        1.1373392110u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3973133608u,        1.1546950559u)--(        0.4315138453u,        1.1373392110u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4623000928u,        1.1602107096u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4623000928u,        1.1602107096u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3973133608u,        1.1546950559u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3973133608u,        1.1546950559u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4315138453u,        1.1373392110u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4315138453u,        1.1373392110u) withcolor MidnightBlue;
draw (        0.3021129054u,        1.1149840488u)--(        0.3400038733u,        1.1209147352u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3435696844u,        1.1591009039u)--(        0.3400038733u,        1.1209147352u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3021129054u,        1.1149840488u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3021129054u,        1.1149840488u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3435696844u,        1.1591009039u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3435696844u,        1.1591009039u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3400038733u,        1.1209147352u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3400038733u,        1.1209147352u) withcolor MidnightBlue;
draw (        0.3021129054u,        1.1149840488u)--(        0.2701087073u,        1.1361176664u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2347964888u,        1.1211528506u)--(        0.2701087073u,        1.1361176664u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3021129054u,        1.1149840488u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3021129054u,        1.1149840488u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2347964888u,        1.1211528506u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2347964888u,        1.1211528506u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2701087073u,        1.1361176664u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2701087073u,        1.1361176664u) withcolor MidnightBlue;
draw (        0.3021129054u,        1.1149840488u)--(        0.2701128751u,        1.1361239766u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2836255647u,        1.1720169523u)--(        0.2701128751u,        1.1361239766u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3021129054u,        1.1149840488u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3021129054u,        1.1149840488u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2836255647u,        1.1720169523u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2836255647u,        1.1720169523u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2701128751u,        1.1361239766u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2701128751u,        1.1361239766u) withcolor MidnightBlue;
draw (        0.3435696844u,        1.1591009039u)--(        0.3682129872u,        1.1297137279u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3973133608u,        1.1546950559u)--(        0.3682129872u,        1.1297137279u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3435696844u,        1.1591009039u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3435696844u,        1.1591009039u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3973133608u,        1.1546950559u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3973133608u,        1.1546950559u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3682129872u,        1.1297137279u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3682129872u,        1.1297137279u) withcolor MidnightBlue;
draw (        0.1657595975u,        1.1099734207u)--(        0.1500000000u,        1.0750000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.1657595975u,        1.1099734207u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1657595975u,        1.1099734207u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1500000000u,        1.0750000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1500000000u,        1.0750000000u) withcolor Yellow;
draw (        0.1657595975u,        1.1099734207u)--(        0.1977600318u,        1.1311127370u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2347964888u,        1.1211528506u)--(        0.1977600318u,        1.1311127370u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.1657595975u,        1.1099734207u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1657595975u,        1.1099734207u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2347964888u,        1.1211528506u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2347964888u,        1.1211528506u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1977600318u,        1.1311127370u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1977600318u,        1.1311127370u) withcolor MidnightBlue;
draw (        0.1657595975u,        1.1099734207u)--(        0.1977552435u,        1.1311199837u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1679437295u,        1.1552482269u)--(        0.1977552435u,        1.1311199837u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.1657595975u,        1.1099734207u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1657595975u,        1.1099734207u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1679437295u,        1.1552482269u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1679437295u,        1.1552482269u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1977552435u,        1.1311199837u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1977552435u,        1.1311199837u) withcolor MidnightBlue;
draw (        0.2347964888u,        1.1211528506u)--(        0.2483173214u,        1.1570427597u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2836255647u,        1.1720169523u)--(        0.2483173214u,        1.1570427597u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2347964888u,        1.1211528506u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2347964888u,        1.1211528506u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2836255647u,        1.1720169523u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2836255647u,        1.1720169523u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2483173214u,        1.1570427597u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2483173214u,        1.1570427597u) withcolor MidnightBlue;
draw (        0.2347964888u,        1.1211528506u)--(        0.2701047321u,        1.1361270432u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2836255647u,        1.1720169523u)--(        0.2701047321u,        1.1361270432u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2347964888u,        1.1211528506u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2347964888u,        1.1211528506u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2836255647u,        1.1720169523u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2836255647u,        1.1720169523u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2701047321u,        1.1361270432u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2701047321u,        1.1361270432u) withcolor MidnightBlue;
draw (        0.2347964888u,        1.1211528506u)--(        0.1977653435u,        1.1311324681u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1679437295u,        1.1552482269u)--(        0.1977653435u,        1.1311324681u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2347964888u,        1.1211528506u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2347964888u,        1.1211528506u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1679437295u,        1.1552482269u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1679437295u,        1.1552482269u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1977653435u,        1.1311324681u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1977653435u,        1.1311324681u) withcolor MidnightBlue;
draw (        0.2347964888u,        1.1211528506u)--(        0.2483049549u,        1.1570474161u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2180866253u,        1.1806641696u)--(        0.2483049549u,        1.1570474161u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2347964888u,        1.1211528506u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2347964888u,        1.1211528506u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2180866253u,        1.1806641696u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2180866253u,        1.1806641696u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2483049549u,        1.1570474161u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2483049549u,        1.1570474161u) withcolor MidnightBlue;
draw (        0.2836255647u,        1.1720169523u)--(        0.2483116924u,        1.1570560395u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2180866253u,        1.1806641696u)--(        0.2483116924u,        1.1570560395u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2836255647u,        1.1720169523u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2836255647u,        1.1720169523u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2180866253u,        1.1806641696u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2180866253u,        1.1806641696u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2483116924u,        1.1570560395u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2483116924u,        1.1570560395u) withcolor MidnightBlue;
draw (        0.2836255647u,        1.1720169523u)--(        0.2534004976u,        1.1956250824u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2180866253u,        1.1806641696u)--(        0.2534004976u,        1.1956250824u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2836255647u,        1.1720169523u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2836255647u,        1.1720169523u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2180866253u,        1.1806641696u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2180866253u,        1.1806641696u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2534004976u,        1.1956250824u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2534004976u,        1.1956250824u) withcolor Green;
% r =         0.0383522944
endfig;
end;
