prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(40);
u := 20cm;
path pols[];
pair p[];
R := .0484156015u;
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

x1 = .6031228179u;
y1 = .3079296637u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .6195605175u;
y2 = .4070965170u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .6616476571u;
y3 = .4615434969u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .7274165567u;
y4 = .5223518138u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .7393096202u;
y5 = .5939143440u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .5889949086u;
y6 = .3899795200u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .5984839760u;
y7 = .4737192742u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .6573239374u;
y8 = .5470822961u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .6871617227u;
y9 = .6279390532u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = .5789307692u;
y10 = .5449923544u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;

x11 = .6115807409u;
y11 = .6224994522u;
draw fullcircle scaled D shifted z11 withcolor RoyalBlue withpen pencircle scaled 1.0;

x12 = .6439810443u;
y12 = .6966223026u;
draw fullcircle scaled D shifted z12 withcolor RoyalBlue withpen pencircle scaled 1.0;

x13 = .5431856958u;
y13 = .6122800320u;
draw fullcircle scaled D shifted z13 withcolor RoyalBlue withpen pencircle scaled 1.0;

x14 = .5641917275u;
y14 = .6986298965u;
draw fullcircle scaled D shifted z14 withcolor RoyalBlue withpen pencircle scaled 1.0;

x15 = .4965121683u;
y15 = .7566219848u;
draw fullcircle scaled D shifted z15 withcolor RoyalBlue withpen pencircle scaled 1.0;

x16 = .5304205881u;
y16 = .8338356732u;
draw fullcircle scaled D shifted z16 withcolor RoyalBlue withpen pencircle scaled 1.0;

x17 = .5212802932u;
y17 = .9098514635u;
draw fullcircle scaled D shifted z17 withcolor RoyalBlue withpen pencircle scaled 1.0;

x18 = .5747025284u;
y18 = .9638586713u;
draw fullcircle scaled D shifted z18 withcolor RoyalBlue withpen pencircle scaled 1.0;

x19 = .6437826047u;
y19 = 1.0157010766u;
draw fullcircle scaled D shifted z19 withcolor RoyalBlue withpen pencircle scaled 1.0;

x20 = .4246904029u;
y20 = .8037876798u;
draw fullcircle scaled D shifted z20 withcolor RoyalBlue withpen pencircle scaled 1.0;

x21 = .4493930627u;
y21 = .8794380029u;
draw fullcircle scaled D shifted z21 withcolor RoyalBlue withpen pencircle scaled 1.0;

x22 = .4796016927u;
y22 = .9690559437u;
draw fullcircle scaled D shifted z22 withcolor RoyalBlue withpen pencircle scaled 1.0;

x23 = .5491006661u;
y23 = 1.0183351186u;
draw fullcircle scaled D shifted z23 withcolor RoyalBlue withpen pencircle scaled 1.0;

x24 = .6060816716u;
y24 = 1.0743713553u;
draw fullcircle scaled D shifted z24 withcolor RoyalBlue withpen pencircle scaled 1.0;

x25 = .3641042609u;
y25 = .8625167960u;
draw fullcircle scaled D shifted z25 withcolor RoyalBlue withpen pencircle scaled 1.0;

x26 = .4069169276u;
y26 = .9344484037u;
draw fullcircle scaled D shifted z26 withcolor RoyalBlue withpen pencircle scaled 1.0;

x27 = .4700378636u;
y27 = 1.0495121219u;
draw fullcircle scaled D shifted z27 withcolor RoyalBlue withpen pencircle scaled 1.0;

x28 = .4957819820u;
y28 = 1.1294722076u;
draw fullcircle scaled D shifted z28 withcolor RoyalBlue withpen pencircle scaled 1.0;

x29 = .5470599604u;
y29 = 1.1817497480u;
draw fullcircle scaled D shifted z29 withcolor RoyalBlue withpen pencircle scaled 1.0;

x30 = .3272775826u;
y30 = .9375430627u;
draw fullcircle scaled D shifted z30 withcolor RoyalBlue withpen pencircle scaled 1.0;

x31 = .4006725590u;
y31 = .9992416516u;
draw fullcircle scaled D shifted z31 withcolor RoyalBlue withpen pencircle scaled 1.0;

x32 = .4143299544u;
y32 = 1.0810460133u;
draw fullcircle scaled D shifted z32 withcolor RoyalBlue withpen pencircle scaled 1.0;

x33 = .3420411449u;
y33 = 1.0183825565u;
draw fullcircle scaled D shifted z33 withcolor RoyalBlue withpen pencircle scaled 1.0;

x34 = .3332820222u;
y34 = 1.0821095227u;
draw fullcircle scaled D shifted z34 withcolor RoyalBlue withpen pencircle scaled 1.0;

x35 = .4259844032u;
y35 = 1.1555603138u;
draw fullcircle scaled D shifted z35 withcolor RoyalBlue withpen pencircle scaled 1.0;

x36 = .1707275287u;
y36 = 1.1186267885u;
draw fullcircle scaled D shifted z36 withcolor RoyalBlue withpen pencircle scaled 1.0;

x37 = .2621094494u;
y37 = 1.1124980252u;
draw fullcircle scaled D shifted z37 withcolor RoyalBlue withpen pencircle scaled 1.0;

x38 = .3500207403u;
y38 = 1.1499461824u;
draw fullcircle scaled D shifted z38 withcolor RoyalBlue withpen pencircle scaled 1.0;

x39 = .1907049176u;
y39 = 1.1715586084u;
draw fullcircle scaled D shifted z39 withcolor RoyalBlue withpen pencircle scaled 1.0;

x40 = .2731216750u;
y40 = 1.1764665359u;
draw fullcircle scaled D shifted z40 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.6031228179u,        0.3079296637u)--(        0.6100000000u,        0.2600000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6031228179u,        0.3079296637u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6031228179u,        0.3079296637u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6100000000u,        0.2600000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6100000000u,        0.2600000000u) withcolor Yellow;
draw (        0.6616476571u,        0.4615434969u)--(        0.6750000000u,        0.4150000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6616476571u,        0.4615434969u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6616476571u,        0.4615434969u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6750000000u,        0.4150000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6750000000u,        0.4150000000u) withcolor Yellow;
draw (        0.6616476571u,        0.4615434969u)--(        0.6820461653u,        0.5054521672u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7274165567u,        0.5223518138u)--(        0.6820461653u,        0.5054521672u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6616476571u,        0.4615434969u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6616476571u,        0.4615434969u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7274165567u,        0.5223518138u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7274165567u,        0.5223518138u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6820461653u,        0.5054521672u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6820461653u,        0.5054521672u) withcolor MidnightBlue;
draw (        0.6616476571u,        0.4615434969u)--(        0.7070180484u,        0.4784431434u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7274165567u,        0.5223518138u)--(        0.7070180484u,        0.4784431434u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6616476571u,        0.4615434969u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6616476571u,        0.4615434969u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7274165567u,        0.5223518138u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7274165567u,        0.5223518138u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7070180484u,        0.4784431434u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7070180484u,        0.4784431434u) withcolor Green;
draw (        0.6616476571u,        0.4615434969u)--(        0.6369155752u,        0.5031655414u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5984839760u,        0.4737192742u)--(        0.6369155752u,        0.5031655414u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6616476571u,        0.4615434969u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6616476571u,        0.4615434969u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5984839760u,        0.4737192742u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5984839760u,        0.4737192742u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6369155752u,        0.5031655414u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6369155752u,        0.5031655414u) withcolor MidnightBlue;
draw (        0.6616476571u,        0.4615434969u)--(        0.6369275367u,        0.5031726467u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6573239374u,        0.5470822961u)--(        0.6369275367u,        0.5031726467u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6616476571u,        0.4615434969u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6616476571u,        0.4615434969u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6369275367u,        0.5031726467u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6369275367u,        0.5031726467u) withcolor MidnightBlue;
draw (        0.6616476571u,        0.4615434969u)--(        0.6820440577u,        0.5054531462u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6573239374u,        0.5470822961u)--(        0.6820440577u,        0.5054531462u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6616476571u,        0.4615434969u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6616476571u,        0.4615434969u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6820440577u,        0.5054531462u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6820440577u,        0.5054531462u) withcolor MidnightBlue;
draw (        0.7274165567u,        0.5223518138u)--(        0.7649980299u,        0.5528756297u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7393096202u,        0.5939143440u)--(        0.7649980299u,        0.5528756297u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7274165567u,        0.5223518138u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7274165567u,        0.5223518138u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7393096202u,        0.5939143440u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7393096202u,        0.5939143440u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7649980299u,        0.5528756297u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7649980299u,        0.5528756297u) withcolor Green;
draw (        0.7274165567u,        0.5223518138u)--(        0.6820454850u,        0.5054539938u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6573239374u,        0.5470822961u)--(        0.6820454850u,        0.5054539938u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7274165567u,        0.5223518138u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7274165567u,        0.5223518138u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6820454850u,        0.5054539938u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6820454850u,        0.5054539938u) withcolor MidnightBlue;
draw (        0.7393096202u,        0.5939143440u)--(        0.7850000000u,        0.6100000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7393096202u,        0.5939143440u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7393096202u,        0.5939143440u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7850000000u,        0.6100000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7850000000u,        0.6100000000u) withcolor Yellow;
draw (        0.7393096202u,        0.5939143440u)--(        0.6929897051u,        0.5798240657u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6573239374u,        0.5470822961u)--(        0.6929897051u,        0.5798240657u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7393096202u,        0.5939143440u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7393096202u,        0.5939143440u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6929897051u,        0.5798240657u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6929897051u,        0.5798240657u) withcolor MidnightBlue;
draw (        0.7393096202u,        0.5939143440u)--(        0.6929746428u,        0.5798736766u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6871617227u,        0.6279390532u)--(        0.6929746428u,        0.5798736766u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7393096202u,        0.5939143440u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7393096202u,        0.5939143440u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6871617227u,        0.6279390532u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6871617227u,        0.6279390532u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6929746428u,        0.5798736766u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6929746428u,        0.5798736766u) withcolor MidnightBlue;
draw (        0.5984839760u,        0.4737192742u)--(        0.6189074246u,        0.5176163494u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6573239374u,        0.5470822961u)--(        0.6189074246u,        0.5176163494u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5984839760u,        0.4737192742u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5984839760u,        0.4737192742u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6189074246u,        0.5176163494u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6189074246u,        0.5176163494u) withcolor MidnightBlue;
draw (        0.5984839760u,        0.4737192742u)--(        0.6369004888u,        0.5031852208u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6573239374u,        0.5470822961u)--(        0.6369004888u,        0.5031852208u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5984839760u,        0.4737192742u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5984839760u,        0.4737192742u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6369004888u,        0.5031852208u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6369004888u,        0.5031852208u) withcolor MidnightBlue;
draw (        0.5984839760u,        0.4737192742u)--(        0.6188741282u,        0.5176318255u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5789307692u,        0.5449923544u)--(        0.6188741282u,        0.5176318255u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5984839760u,        0.4737192742u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5984839760u,        0.4737192742u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5789307692u,        0.5449923544u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5789307692u,        0.5449923544u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6188741282u,        0.5176318255u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6188741282u,        0.5176318255u) withcolor MidnightBlue;
draw (        0.6573239374u,        0.5470822961u)--(        0.6515385091u,        0.5951509895u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6871617227u,        0.6279390532u)--(        0.6515385091u,        0.5951509895u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6871617227u,        0.6279390532u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6871617227u,        0.6279390532u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6515385091u,        0.5951509895u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6515385091u,        0.5951509895u) withcolor MidnightBlue;
draw (        0.6573239374u,        0.5470822961u)--(        0.6929471510u,        0.5798703598u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6871617227u,        0.6279390532u)--(        0.6929471510u,        0.5798703598u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6871617227u,        0.6279390532u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6871617227u,        0.6279390532u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6929471510u,        0.5798703598u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6929471510u,        0.5798703598u) withcolor MidnightBlue;
draw (        0.6573239374u,        0.5470822961u)--(        0.6188842428u,        0.5176465976u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5789307692u,        0.5449923544u)--(        0.6188842428u,        0.5176465976u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5789307692u,        0.5449923544u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5789307692u,        0.5449923544u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6188842428u,        0.5176465976u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6188842428u,        0.5176465976u) withcolor MidnightBlue;
draw (        0.6573239374u,        0.5470822961u)--(        0.6173704638u,        0.5744280529u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5789307692u,        0.5449923544u)--(        0.6173704638u,        0.5744280529u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5789307692u,        0.5449923544u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5789307692u,        0.5449923544u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6173704638u,        0.5744280529u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6173704638u,        0.5744280529u) withcolor MidnightBlue;
draw (        0.6573239374u,        0.5470822961u)--(        0.6173728682u,        0.5744315655u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6115807409u,        0.6224994522u)--(        0.6173728682u,        0.5744315655u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6115807409u,        0.6224994522u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6115807409u,        0.6224994522u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6173728682u,        0.5744315655u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6173728682u,        0.5744315655u) withcolor MidnightBlue;
draw (        0.6573239374u,        0.5470822961u)--(        0.6515318101u,        0.5951501828u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6115807409u,        0.6224994522u)--(        0.6515318101u,        0.5951501828u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6573239374u,        0.5470822961u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6115807409u,        0.6224994522u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6115807409u,        0.6224994522u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6515318101u,        0.5951501828u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6515318101u,        0.5951501828u) withcolor MidnightBlue;
draw (        0.6871617227u,        0.6279390532u)--(        0.6515349841u,        0.5951548198u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6115807409u,        0.6224994522u)--(        0.6515349841u,        0.5951548198u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6871617227u,        0.6279390532u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6871617227u,        0.6279390532u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6115807409u,        0.6224994522u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6115807409u,        0.6224994522u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6515349841u,        0.5951548198u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6515349841u,        0.5951548198u) withcolor MidnightBlue;
draw (        0.5789307692u,        0.5449923544u)--(        0.6173682040u,        0.5744310037u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6115807409u,        0.6224994522u)--(        0.6173682040u,        0.5744310037u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5789307692u,        0.5449923544u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5789307692u,        0.5449923544u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6115807409u,        0.6224994522u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6115807409u,        0.6224994522u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6173682040u,        0.5744310037u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6173682040u,        0.5744310037u) withcolor MidnightBlue;
draw (        0.6115807409u,        0.6224994522u)--(        0.6033987832u,        0.6702186940u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6439810443u,        0.6966223026u)--(        0.6033987832u,        0.6702186940u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6115807409u,        0.6224994522u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6115807409u,        0.6224994522u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6439810443u,        0.6966223026u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6439810443u,        0.6966223026u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6033987832u,        0.6702186940u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6033987832u,        0.6702186940u) withcolor MidnightBlue;
draw (        0.6115807409u,        0.6224994522u)--(        0.5723751441u,        0.6509070507u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5431856958u,        0.6122800320u)--(        0.5723751441u,        0.6509070507u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6115807409u,        0.6224994522u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6115807409u,        0.6224994522u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5431856958u,        0.6122800320u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5431856958u,        0.6122800320u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5723751441u,        0.6509070507u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5723751441u,        0.6509070507u) withcolor MidnightBlue;
draw (        0.6115807409u,        0.6224994522u)--(        0.5723783336u,        0.6509114520u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5641917275u,        0.6986298965u)--(        0.5723783336u,        0.6509114520u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6115807409u,        0.6224994522u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6115807409u,        0.6224994522u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5641917275u,        0.6986298965u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5641917275u,        0.6986298965u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5723783336u,        0.6509114520u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5723783336u,        0.6509114520u) withcolor MidnightBlue;
draw (        0.6115807409u,        0.6224994522u)--(        0.6033941348u,        0.6702178967u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5641917275u,        0.6986298965u)--(        0.6033941348u,        0.6702178967u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6115807409u,        0.6224994522u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6115807409u,        0.6224994522u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5641917275u,        0.6986298965u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5641917275u,        0.6986298965u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6033941348u,        0.6702178967u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6033941348u,        0.6702178967u) withcolor MidnightBlue;
draw (        0.6439810443u,        0.6966223026u)--(        0.6033968576u,        0.6702216539u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5641917275u,        0.6986298965u)--(        0.6033968576u,        0.6702216539u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6439810443u,        0.6966223026u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6439810443u,        0.6966223026u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5641917275u,        0.6986298965u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5641917275u,        0.6986298965u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6033968576u,        0.6702216539u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6033968576u,        0.6702216539u) withcolor MidnightBlue;
draw (        0.5431856958u,        0.6122800320u)--(        0.5350000000u,        0.6600000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5431856958u,        0.6122800320u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5431856958u,        0.6122800320u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5350000000u,        0.6600000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5350000000u,        0.6600000000u) withcolor Yellow;
draw (        0.5431856958u,        0.6122800320u)--(        0.5723709866u,        0.6509101921u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5641917275u,        0.6986298965u)--(        0.5723709866u,        0.6509101921u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5431856958u,        0.6122800320u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5431856958u,        0.6122800320u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5641917275u,        0.6986298965u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5641917275u,        0.6986298965u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5723709866u,        0.6509101921u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5723709866u,        0.6509101921u) withcolor MidnightBlue;
draw (        0.5641917275u,        0.6986298965u)--(        0.5350000000u,        0.6600000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5641917275u,        0.6986298965u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5641917275u,        0.6986298965u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5350000000u,        0.6600000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5350000000u,        0.6600000000u) withcolor Yellow;
draw (        0.4965121683u,        0.7566219848u)--(        0.4750000000u,        0.8000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4965121683u,        0.7566219848u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4965121683u,        0.7566219848u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4750000000u,        0.8000000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4750000000u,        0.8000000000u) withcolor Yellow;
draw (        0.5304205881u,        0.8338356732u)--(        0.5650000000u,        0.8000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5304205881u,        0.8338356732u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5304205881u,        0.8338356732u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5650000000u,        0.8000000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5650000000u,        0.8000000000u) withcolor Yellow;
draw (        0.5304205881u,        0.8338356732u)--(        0.4964215014u,        0.8683049722u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5212802932u,        0.9098514635u)--(        0.4964215014u,        0.8683049722u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5304205881u,        0.8338356732u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5304205881u,        0.8338356732u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5212802932u,        0.9098514635u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5212802932u,        0.9098514635u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4964215014u,        0.8683049722u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4964215014u,        0.8683049722u) withcolor MidnightBlue;
draw (        0.5212802932u,        0.9098514635u)--(        0.5266465258u,        0.9579687579u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5747025284u,        0.9638586713u)--(        0.5266465258u,        0.9579687579u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5212802932u,        0.9098514635u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5212802932u,        0.9098514635u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5747025284u,        0.9638586713u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5747025284u,        0.9638586713u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5266465258u,        0.9579687579u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5266465258u,        0.9579687579u) withcolor MidnightBlue;
draw (        0.5212802932u,        0.9098514635u)--(        0.4965001561u,        0.8682580112u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4493930627u,        0.8794380029u)--(        0.4965001561u,        0.8682580112u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5212802932u,        0.9098514635u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5212802932u,        0.9098514635u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4493930627u,        0.8794380029u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4493930627u,        0.8794380029u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4965001561u,        0.8682580112u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4965001561u,        0.8682580112u) withcolor MidnightBlue;
draw (        0.5212802932u,        0.9098514635u)--(        0.4741731998u,        0.9210314552u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4493930627u,        0.8794380029u)--(        0.4741731998u,        0.9210314552u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5212802932u,        0.9098514635u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5212802932u,        0.9098514635u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4493930627u,        0.8794380029u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4493930627u,        0.8794380029u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4741731998u,        0.9210314552u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4741731998u,        0.9210314552u) withcolor MidnightBlue;
draw (        0.5212802932u,        0.9098514635u)--(        0.4741534400u,        0.9209478673u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4796016927u,        0.9690559437u)--(        0.4741534400u,        0.9209478673u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5212802932u,        0.9098514635u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5212802932u,        0.9098514635u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4796016927u,        0.9690559437u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4796016927u,        0.9690559437u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4741534400u,        0.9209478673u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4741534400u,        0.9209478673u) withcolor MidnightBlue;
draw (        0.5212802932u,        0.9098514635u)--(        0.5267285459u,        0.9579595399u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4796016927u,        0.9690559437u)--(        0.5267285459u,        0.9579595399u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5212802932u,        0.9098514635u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5212802932u,        0.9098514635u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4796016927u,        0.9690559437u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4796016927u,        0.9690559437u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5267285459u,        0.9579595399u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5267285459u,        0.9579595399u) withcolor MidnightBlue;
draw (        0.5747025284u,        0.9638586713u)--(        0.5600000000u,        1.0100000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5747025284u,        0.9638586713u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5747025284u,        0.9638586713u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5600000000u,        1.0100000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5600000000u,        1.0100000000u) withcolor Yellow;
draw (        0.5747025284u,        0.9638586713u)--(        0.5275783084u,        0.9749662520u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5491006661u,        1.0183351186u)--(        0.5275783084u,        0.9749662520u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5747025284u,        0.9638586713u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5747025284u,        0.9638586713u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5491006661u,        1.0183351186u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5491006661u,        1.0183351186u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5275783084u,        0.9749662520u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5275783084u,        0.9749662520u) withcolor MidnightBlue;
draw (        0.6437826047u,        1.0157010766u)--(        0.6900000000u,        1.0300000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6437826047u,        1.0157010766u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6437826047u,        1.0157010766u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6900000000u,        1.0300000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6900000000u,        1.0300000000u) withcolor Yellow;
draw (        0.6437826047u,        1.0157010766u)--(        0.5966748735u,        1.0268783809u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6060816716u,        1.0743713553u)--(        0.5966748735u,        1.0268783809u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6437826047u,        1.0157010766u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6437826047u,        1.0157010766u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6060816716u,        1.0743713553u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6060816716u,        1.0743713553u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5966748735u,        1.0268783809u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5966748735u,        1.0268783809u) withcolor MidnightBlue;
draw (        0.4246904029u,        0.8037876798u)--(        0.4108220900u,        0.8501745353u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4493930627u,        0.8794380029u)--(        0.4108220900u,        0.8501745353u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4246904029u,        0.8037876798u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4246904029u,        0.8037876798u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4493930627u,        0.8794380029u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4493930627u,        0.8794380029u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4108220900u,        0.8501745353u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4108220900u,        0.8501745353u) withcolor MidnightBlue;
draw (        0.4493930627u,        0.8794380029u)--(        0.4548380451u,        0.9275464495u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4069169276u,        0.9344484037u)--(        0.4548380451u,        0.9275464495u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4493930627u,        0.8794380029u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4493930627u,        0.8794380029u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4069169276u,        0.9344484037u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4069169276u,        0.9344484037u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4548380451u,        0.9275464495u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4548380451u,        0.9275464495u) withcolor MidnightBlue;
draw (        0.4796016927u,        0.9690559437u)--(        0.5276601580u,        0.9749257293u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5491006661u,        1.0183351186u)--(        0.5276601580u,        0.9749257293u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4796016927u,        0.9690559437u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4796016927u,        0.9690559437u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5491006661u,        1.0183351186u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5491006661u,        1.0183351186u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5276601580u,        0.9749257293u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5276601580u,        0.9749257293u) withcolor MidnightBlue;
draw (        0.4796016927u,        0.9690559437u)--(        0.4548255322u,        0.9274601225u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4069169276u,        0.9344484037u)--(        0.4548255322u,        0.9274601225u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4796016927u,        0.9690559437u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4796016927u,        0.9690559437u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4069169276u,        0.9344484037u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4069169276u,        0.9344484037u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4548255322u,        0.9274601225u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4548255322u,        0.9274601225u) withcolor MidnightBlue;
draw (        0.4796016927u,        0.9690559437u)--(        0.4484927789u,        1.0061545414u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4700378636u,        1.0495121219u)--(        0.4484927789u,        1.0061545414u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4796016927u,        0.9690559437u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4796016927u,        0.9690559437u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4700378636u,        1.0495121219u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4700378636u,        1.0495121219u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4484927789u,        1.0061545414u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4484927789u,        1.0061545414u) withcolor MidnightBlue;
draw (        0.4796016927u,        0.9690559437u)--(        0.5011467774u,        1.0124135242u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4700378636u,        1.0495121219u)--(        0.5011467774u,        1.0124135242u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4796016927u,        0.9690559437u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4796016927u,        0.9690559437u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4700378636u,        1.0495121219u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4700378636u,        1.0495121219u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5011467774u,        1.0124135242u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5011467774u,        1.0124135242u) withcolor MidnightBlue;
draw (        0.5491006661u,        1.0183351186u)--(        0.5967595378u,        1.0268616900u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6060816716u,        1.0743713553u)--(        0.5967595378u,        1.0268616900u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5491006661u,        1.0183351186u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5491006661u,        1.0183351186u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6060816716u,        1.0743713553u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6060816716u,        1.0743713553u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5967595378u,        1.0268616900u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5967595378u,        1.0268616900u) withcolor MidnightBlue;
draw (        0.6060816716u,        1.0743713553u)--(        0.5900000000u,        1.1200000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6060816716u,        1.0743713553u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6060816716u,        1.0743713553u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5900000000u,        1.1200000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5900000000u,        1.1200000000u) withcolor Yellow;
draw (        0.3641042609u,        0.8625167960u)--(        0.3550000000u,        0.8150000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3641042609u,        0.8625167960u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3641042609u,        0.8625167960u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3550000000u,        0.8150000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3550000000u,        0.8150000000u) withcolor Yellow;
draw (        0.4069169276u,        0.9344484037u)--(        0.3681649286u,        0.9634717226u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3272775826u,        0.9375430627u)--(        0.3681649286u,        0.9634717226u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4069169276u,        0.9344484037u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4069169276u,        0.9344484037u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3272775826u,        0.9375430627u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3272775826u,        0.9375430627u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3681649286u,        0.9634717226u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3681649286u,        0.9634717226u) withcolor MidnightBlue;
draw (        0.4069169276u,        0.9344484037u)--(        0.3681162052u,        0.9634065530u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4006725590u,        0.9992416516u)--(        0.3681162052u,        0.9634065530u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4069169276u,        0.9344484037u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4069169276u,        0.9344484037u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4006725590u,        0.9992416516u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4006725590u,        0.9992416516u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3681162052u,        0.9634065530u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3681162052u,        0.9634065530u) withcolor MidnightBlue;
draw (        0.4700378636u,        1.0495121219u)--(        0.4599860647u,        1.0968727786u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4957819820u,        1.1294722076u)--(        0.4599860647u,        1.0968727786u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4700378636u,        1.0495121219u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4700378636u,        1.0495121219u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4957819820u,        1.1294722076u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4957819820u,        1.1294722076u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4599860647u,        1.0968727786u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4599860647u,        1.0968727786u) withcolor MidnightBlue;
draw (        0.4700378636u,        1.0495121219u)--(        0.4600788683u,        1.0968923802u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4143299544u,        1.0810460133u)--(        0.4600788683u,        1.0968923802u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4700378636u,        1.0495121219u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4700378636u,        1.0495121219u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4143299544u,        1.0810460133u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4143299544u,        1.0810460133u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4600788683u,        1.0968923802u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4600788683u,        1.0968923802u) withcolor MidnightBlue;
draw (        0.5470599604u,        1.1817497480u)--(        0.5500000000u,        1.2300000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5470599604u,        1.1817497480u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5470599604u,        1.1817497480u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5500000000u,        1.2300000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5500000000u,        1.2300000000u) withcolor Yellow;
draw (        0.3272775826u,        0.9375430627u)--(        0.3598517462u,        0.9733619730u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3420411449u,        1.0183825565u)--(        0.3598517462u,        0.9733619730u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3272775826u,        0.9375430627u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3272775826u,        0.9375430627u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3420411449u,        1.0183825565u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3420411449u,        1.0183825565u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3598517462u,        0.9733619730u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3598517462u,        0.9733619730u) withcolor MidnightBlue;
draw (        0.4006725590u,        0.9992416516u)--(        0.3828533836u,        1.0442588422u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4143299544u,        1.0810460133u)--(        0.3828533836u,        1.0442588422u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4006725590u,        0.9992416516u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4006725590u,        0.9992416516u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4143299544u,        1.0810460133u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4143299544u,        1.0810460133u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3828533836u,        1.0442588422u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3828533836u,        1.0442588422u) withcolor MidnightBlue;
draw (        0.4006725590u,        0.9992416516u)--(        0.3597736901u,        0.9733311713u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3420411449u,        1.0183825565u)--(        0.3597736901u,        0.9733311713u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4006725590u,        0.9992416516u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4006725590u,        0.9992416516u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3420411449u,        1.0183825565u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3420411449u,        1.0183825565u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3597736901u,        0.9733311713u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3597736901u,        0.9733311713u) withcolor MidnightBlue;
draw (        0.4006725590u,        0.9992416516u)--(        0.3829400138u,        1.0442930368u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3420411449u,        1.0183825565u)--(        0.3829400138u,        1.0442930368u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4006725590u,        0.9992416516u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4006725590u,        0.9992416516u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3420411449u,        1.0183825565u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3420411449u,        1.0183825565u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3829400138u,        1.0442930368u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3829400138u,        1.0442930368u) withcolor MidnightBlue;
draw (        0.4143299544u,        1.0810460133u)--(        0.3734584478u,        1.0550923929u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3332820222u,        1.0821095227u)--(        0.3734584478u,        1.0550923929u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4143299544u,        1.0810460133u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4143299544u,        1.0810460133u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3332820222u,        1.0821095227u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3332820222u,        1.0821095227u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3734584478u,        1.0550923929u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3734584478u,        1.0550923929u) withcolor MidnightBlue;
draw (        0.4143299544u,        1.0810460133u)--(        0.3741535289u,        1.1080631430u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3332820222u,        1.0821095227u)--(        0.3741535289u,        1.1080631430u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4143299544u,        1.0810460133u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4143299544u,        1.0810460133u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3332820222u,        1.0821095227u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3332820222u,        1.0821095227u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3741535289u,        1.1080631430u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3741535289u,        1.1080631430u) withcolor MidnightBlue;
draw (        0.4143299544u,        1.0810460133u)--(        0.3901571420u,        1.1229953352u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4259844032u,        1.1555603138u)--(        0.3901571420u,        1.1229953352u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4143299544u,        1.0810460133u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4143299544u,        1.0810460133u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4259844032u,        1.1555603138u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4259844032u,        1.1555603138u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3901571420u,        1.1229953352u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3901571420u,        1.1229953352u) withcolor MidnightBlue;
draw (        0.3420411449u,        1.0183825565u)--(        0.3735130292u,        1.0551737371u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3332820222u,        1.0821095227u)--(        0.3735130292u,        1.0551737371u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3420411449u,        1.0183825565u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3420411449u,        1.0183825565u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3332820222u,        1.0821095227u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3332820222u,        1.0821095227u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3735130292u,        1.0551737371u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3735130292u,        1.0551737371u) withcolor MidnightBlue;
draw (        0.3332820222u,        1.0821095227u)--(        0.3091226297u,        1.1240665748u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2621094494u,        1.1124980252u)--(        0.3091226297u,        1.1240665748u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3332820222u,        1.0821095227u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3332820222u,        1.0821095227u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2621094494u,        1.1124980252u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2621094494u,        1.1124980252u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3091226297u,        1.1240665748u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3091226297u,        1.1240665748u) withcolor MidnightBlue;
draw (        0.3332820222u,        1.0821095227u)--(        0.3091077511u,        1.1240580039u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3500207403u,        1.1499461824u)--(        0.3091077511u,        1.1240580039u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3332820222u,        1.0821095227u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3332820222u,        1.0821095227u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3500207403u,        1.1499461824u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3500207403u,        1.1499461824u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3091077511u,        1.1240580039u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3091077511u,        1.1240580039u) withcolor MidnightBlue;
draw (        0.3332820222u,        1.0821095227u)--(        0.3741950114u,        1.1079977011u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3500207403u,        1.1499461824u)--(        0.3741950114u,        1.1079977011u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3332820222u,        1.0821095227u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3332820222u,        1.0821095227u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3500207403u,        1.1499461824u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3500207403u,        1.1499461824u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3741950114u,        1.1079977011u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3741950114u,        1.1079977011u) withcolor MidnightBlue;
draw (        0.4259844032u,        1.1555603138u)--(        0.3902057947u,        1.1229418888u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3500207403u,        1.1499461824u)--(        0.3902057947u,        1.1229418888u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4259844032u,        1.1555603138u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4259844032u,        1.1555603138u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3500207403u,        1.1499461824u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3500207403u,        1.1499461824u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3902057947u,        1.1229418888u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3902057947u,        1.1229418888u) withcolor MidnightBlue;
draw (        0.1707275287u,        1.1186267885u)--(        0.2174702104u,        1.1312439265u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2621094494u,        1.1124980252u)--(        0.2174702104u,        1.1312439265u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.1707275287u,        1.1186267885u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1707275287u,        1.1186267885u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2621094494u,        1.1124980252u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2621094494u,        1.1124980252u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2174702104u,        1.1312439265u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2174702104u,        1.1312439265u) withcolor MidnightBlue;
draw (        0.1707275287u,        1.1186267885u)--(        0.2174770612u,        1.1312185192u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1907049176u,        1.1715586084u)--(        0.2174770612u,        1.1312185192u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.1707275287u,        1.1186267885u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1707275287u,        1.1186267885u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1907049176u,        1.1715586084u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1907049176u,        1.1715586084u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2174770612u,        1.1312185192u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2174770612u,        1.1312185192u) withcolor MidnightBlue;
draw (        0.2621094494u,        1.1124980252u)--(        0.3029947122u,        1.1384299698u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3500207403u,        1.1499461824u)--(        0.3029947122u,        1.1384299698u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2621094494u,        1.1124980252u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2621094494u,        1.1124980252u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3500207403u,        1.1499461824u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3500207403u,        1.1499461824u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3029947122u,        1.1384299698u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3029947122u,        1.1384299698u) withcolor MidnightBlue;
draw (        0.2621094494u,        1.1124980252u)--(        0.3091354774u,        1.1240142378u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3500207403u,        1.1499461824u)--(        0.3091354774u,        1.1240142378u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2621094494u,        1.1124980252u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2621094494u,        1.1124980252u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3500207403u,        1.1499461824u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3500207403u,        1.1499461824u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3091354774u,        1.1240142378u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3091354774u,        1.1240142378u) withcolor MidnightBlue;
draw (        0.2621094494u,        1.1124980252u)--(        0.2174528088u,        1.1312024343u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1907049176u,        1.1715586084u)--(        0.2174528088u,        1.1312024343u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2621094494u,        1.1124980252u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2621094494u,        1.1124980252u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1907049176u,        1.1715586084u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1907049176u,        1.1715586084u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2174528088u,        1.1312024343u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2174528088u,        1.1312024343u) withcolor MidnightBlue;
draw (        0.2621094494u,        1.1124980252u)--(        0.3030219007u,        1.1383870536u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2731216750u,        1.1764665359u)--(        0.3030219007u,        1.1383870536u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2621094494u,        1.1124980252u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2621094494u,        1.1124980252u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2731216750u,        1.1764665359u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2731216750u,        1.1764665359u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3030219007u,        1.1383870536u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3030219007u,        1.1383870536u) withcolor MidnightBlue;
draw (        0.3500207403u,        1.1499461824u)--(        0.3030079385u,        1.1383760945u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2731216750u,        1.1764665359u)--(        0.3030079385u,        1.1383760945u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3500207403u,        1.1499461824u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3500207403u,        1.1499461824u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2731216750u,        1.1764665359u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2731216750u,        1.1764665359u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3030079385u,        1.1383760945u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3030079385u,        1.1383760945u) withcolor MidnightBlue;
% r =         0.0484156015
endfig;
end;
