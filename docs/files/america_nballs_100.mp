prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(100);
u := 20cm;
path pols[];
pair p[];
R := .0291721038u;
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

x1 = .5997441687u;
y1 = .2796007626u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .5890920185u;
y2 = .3194557745u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .5974000058u;
y3 = .3568902026u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .6327026846u;
y4 = .3952892494u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .6668063752u;
y5 = .4399177212u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .6931424786u;
y6 = .4792642245u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .7305543567u;
y7 = .5164293290u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .7452099634u;
y8 = .5493636085u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .7608512575u;
y9 = .5936337999u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = .5865417564u;
y10 = .3922151310u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;

x11 = .6192201461u;
y11 = .4385953434u;
draw fullcircle scaled D shifted z11 withcolor RoyalBlue withpen pencircle scaled 1.0;

x12 = .6403541248u;
y12 = .4833237349u;
draw fullcircle scaled D shifted z12 withcolor RoyalBlue withpen pencircle scaled 1.0;

x13 = .6764840238u;
y13 = .5187505238u;
draw fullcircle scaled D shifted z13 withcolor RoyalBlue withpen pencircle scaled 1.0;

x14 = .6980126540u;
y14 = .5570317379u;
draw fullcircle scaled D shifted z14 withcolor RoyalBlue withpen pencircle scaled 1.0;

x15 = .7157844052u;
y15 = .5976011166u;
draw fullcircle scaled D shifted z15 withcolor RoyalBlue withpen pencircle scaled 1.0;

x16 = .7362248211u;
y16 = .6226591835u;
draw fullcircle scaled D shifted z16 withcolor RoyalBlue withpen pencircle scaled 1.0;

x17 = .5807790418u;
y17 = .4433515528u;
draw fullcircle scaled D shifted z17 withcolor RoyalBlue withpen pencircle scaled 1.0;

x18 = .5940369321u;
y18 = .4869891718u;
draw fullcircle scaled D shifted z18 withcolor RoyalBlue withpen pencircle scaled 1.0;

x19 = .6252797313u;
y19 = .5293876300u;
draw fullcircle scaled D shifted z19 withcolor RoyalBlue withpen pencircle scaled 1.0;

x20 = .6475384636u;
y20 = .5652854306u;
draw fullcircle scaled D shifted z20 withcolor RoyalBlue withpen pencircle scaled 1.0;

x21 = .6755140958u;
y21 = .6041810397u;
draw fullcircle scaled D shifted z21 withcolor RoyalBlue withpen pencircle scaled 1.0;

x22 = .6921305001u;
y22 = .6498182465u;
draw fullcircle scaled D shifted z22 withcolor RoyalBlue withpen pencircle scaled 1.0;

x23 = .5755157707u;
y23 = .5321651252u;
draw fullcircle scaled D shifted z23 withcolor RoyalBlue withpen pencircle scaled 1.0;

x24 = .5952690094u;
y24 = .5744715760u;
draw fullcircle scaled D shifted z24 withcolor RoyalBlue withpen pencircle scaled 1.0;

x25 = .6239746084u;
y25 = .6109837654u;
draw fullcircle scaled D shifted z25 withcolor RoyalBlue withpen pencircle scaled 1.0;

x26 = .6501069626u;
y26 = .6407568395u;
draw fullcircle scaled D shifted z26 withcolor RoyalBlue withpen pencircle scaled 1.0;

x27 = .6561892216u;
y27 = .6807257247u;
draw fullcircle scaled D shifted z27 withcolor RoyalBlue withpen pencircle scaled 1.0;

x28 = .5524199618u;
y28 = .5768978324u;
draw fullcircle scaled D shifted z28 withcolor RoyalBlue withpen pencircle scaled 1.0;

x29 = .5778526313u;
y29 = .6244435963u;
draw fullcircle scaled D shifted z29 withcolor RoyalBlue withpen pencircle scaled 1.0;

x30 = .6032885711u;
y30 = .6626856257u;
draw fullcircle scaled D shifted z30 withcolor RoyalBlue withpen pencircle scaled 1.0;

x31 = .6140175721u;
y31 = .7057317745u;
draw fullcircle scaled D shifted z31 withcolor RoyalBlue withpen pencircle scaled 1.0;

x32 = .5352615775u;
y32 = .6239150450u;
draw fullcircle scaled D shifted z32 withcolor RoyalBlue withpen pencircle scaled 1.0;

x33 = .5540480399u;
y33 = .6731980856u;
draw fullcircle scaled D shifted z33 withcolor RoyalBlue withpen pencircle scaled 1.0;

x34 = .5710107301u;
y34 = .7094808588u;
draw fullcircle scaled D shifted z34 withcolor RoyalBlue withpen pencircle scaled 1.0;

x35 = .5198362864u;
y35 = .7176396000u;
draw fullcircle scaled D shifted z35 withcolor RoyalBlue withpen pencircle scaled 1.0;

x36 = .4825763447u;
y36 = .7550513521u;
draw fullcircle scaled D shifted z36 withcolor RoyalBlue withpen pencircle scaled 1.0;

x37 = .5421056577u;
y37 = .8180848237u;
draw fullcircle scaled D shifted z37 withcolor RoyalBlue withpen pencircle scaled 1.0;

x38 = .5299513321u;
y38 = .8687716044u;
draw fullcircle scaled D shifted z38 withcolor RoyalBlue withpen pencircle scaled 1.0;

x39 = .5533000523u;
y39 = .9125635899u;
draw fullcircle scaled D shifted z39 withcolor RoyalBlue withpen pencircle scaled 1.0;

x40 = .5852875073u;
y40 = .9559954384u;
draw fullcircle scaled D shifted z40 withcolor RoyalBlue withpen pencircle scaled 1.0;

x41 = .6289565053u;
y41 = .9882172714u;
draw fullcircle scaled D shifted z41 withcolor RoyalBlue withpen pencircle scaled 1.0;

x42 = .6611669001u;
y42 = 1.0255369832u;
draw fullcircle scaled D shifted z42 withcolor RoyalBlue withpen pencircle scaled 1.0;

x43 = .4280920558u;
y43 = .7709921534u;
draw fullcircle scaled D shifted z43 withcolor RoyalBlue withpen pencircle scaled 1.0;

x44 = .4696277763u;
y44 = .7941263670u;
draw fullcircle scaled D shifted z44 withcolor RoyalBlue withpen pencircle scaled 1.0;

x45 = .4871977584u;
y45 = .8801053301u;
draw fullcircle scaled D shifted z45 withcolor RoyalBlue withpen pencircle scaled 1.0;

x46 = .5115937696u;
y46 = .9219448456u;
draw fullcircle scaled D shifted z46 withcolor RoyalBlue withpen pencircle scaled 1.0;

x47 = .5438378022u;
y47 = .9649398475u;
draw fullcircle scaled D shifted z47 withcolor RoyalBlue withpen pencircle scaled 1.0;

x48 = .5777054129u;
y48 = .9986010128u;
draw fullcircle scaled D shifted z48 withcolor RoyalBlue withpen pencircle scaled 1.0;

x49 = .6101948595u;
y49 = 1.0349547361u;
draw fullcircle scaled D shifted z49 withcolor RoyalBlue withpen pencircle scaled 1.0;

x50 = .6339862043u;
y50 = 1.0604411385u;
draw fullcircle scaled D shifted z50 withcolor RoyalBlue withpen pencircle scaled 1.0;

x51 = .3973659469u;
y51 = .8093476485u;
draw fullcircle scaled D shifted z51 withcolor RoyalBlue withpen pencircle scaled 1.0;

x52 = .4154396362u;
y52 = .8494899169u;
draw fullcircle scaled D shifted z52 withcolor RoyalBlue withpen pencircle scaled 1.0;

x53 = .4412178471u;
y53 = .8752994513u;
draw fullcircle scaled D shifted z53 withcolor RoyalBlue withpen pencircle scaled 1.0;

x54 = .4651478441u;
y54 = .9244416388u;
draw fullcircle scaled D shifted z54 withcolor RoyalBlue withpen pencircle scaled 1.0;

x55 = .4963590693u;
y55 = .9690963120u;
draw fullcircle scaled D shifted z55 withcolor RoyalBlue withpen pencircle scaled 1.0;

x56 = .5266752110u;
y56 = 1.0125072761u;
draw fullcircle scaled D shifted z56 withcolor RoyalBlue withpen pencircle scaled 1.0;

x57 = .5799642533u;
y57 = 1.0430584069u;
draw fullcircle scaled D shifted z57 withcolor RoyalBlue withpen pencircle scaled 1.0;

x58 = .5962564363u;
y58 = 1.0914932443u;
draw fullcircle scaled D shifted z58 withcolor RoyalBlue withpen pencircle scaled 1.0;

x59 = .3604781949u;
y59 = .8436554871u;
draw fullcircle scaled D shifted z59 withcolor RoyalBlue withpen pencircle scaled 1.0;

x60 = .3908073442u;
y60 = .8847964234u;
draw fullcircle scaled D shifted z60 withcolor RoyalBlue withpen pencircle scaled 1.0;

x61 = .4213655894u;
y61 = .9227545450u;
draw fullcircle scaled D shifted z61 withcolor RoyalBlue withpen pencircle scaled 1.0;

x62 = .4492648580u;
y62 = .9707889022u;
draw fullcircle scaled D shifted z62 withcolor RoyalBlue withpen pencircle scaled 1.0;

x63 = .4785308756u;
y63 = 1.0158542810u;
draw fullcircle scaled D shifted z63 withcolor RoyalBlue withpen pencircle scaled 1.0;

x64 = .5081104058u;
y64 = 1.0593808224u;
draw fullcircle scaled D shifted z64 withcolor RoyalBlue withpen pencircle scaled 1.0;

x65 = .5534906856u;
y65 = 1.1611442623u;
draw fullcircle scaled D shifted z65 withcolor RoyalBlue withpen pencircle scaled 1.0;

x66 = .3419840060u;
y66 = .8869404202u;
draw fullcircle scaled D shifted z66 withcolor RoyalBlue withpen pencircle scaled 1.0;

x67 = .3757244663u;
y67 = .9296651377u;
draw fullcircle scaled D shifted z67 withcolor RoyalBlue withpen pencircle scaled 1.0;

x68 = .4074541434u;
y68 = .9723593407u;
draw fullcircle scaled D shifted z68 withcolor RoyalBlue withpen pencircle scaled 1.0;

x69 = .4342184410u;
y69 = 1.0206151214u;
draw fullcircle scaled D shifted z69 withcolor RoyalBlue withpen pencircle scaled 1.0;

x70 = .4670658524u;
y70 = 1.0648141581u;
draw fullcircle scaled D shifted z70 withcolor RoyalBlue withpen pencircle scaled 1.0;

x71 = .4961628118u;
y71 = 1.1110599668u;
draw fullcircle scaled D shifted z71 withcolor RoyalBlue withpen pencircle scaled 1.0;

x72 = .5210640139u;
y72 = 1.1514325278u;
draw fullcircle scaled D shifted z72 withcolor RoyalBlue withpen pencircle scaled 1.0;

x73 = .5376929472u;
y73 = 1.2034909989u;
draw fullcircle scaled D shifted z73 withcolor RoyalBlue withpen pencircle scaled 1.0;

x74 = .3257749524u;
y74 = .9282285454u;
draw fullcircle scaled D shifted z74 withcolor RoyalBlue withpen pencircle scaled 1.0;

x75 = .3586846939u;
y75 = .9719569461u;
draw fullcircle scaled D shifted z75 withcolor RoyalBlue withpen pencircle scaled 1.0;

x76 = .3880162345u;
y76 = 1.0168985547u;
draw fullcircle scaled D shifted z76 withcolor RoyalBlue withpen pencircle scaled 1.0;

x77 = .4183806882u;
y77 = 1.0642703841u;
draw fullcircle scaled D shifted z77 withcolor RoyalBlue withpen pencircle scaled 1.0;

x78 = .4484299957u;
y78 = 1.1091930801u;
draw fullcircle scaled D shifted z78 withcolor RoyalBlue withpen pencircle scaled 1.0;

x79 = .4696310722u;
y79 = 1.1559758622u;
draw fullcircle scaled D shifted z79 withcolor RoyalBlue withpen pencircle scaled 1.0;

x80 = .3141122531u;
y80 = .9699081033u;
draw fullcircle scaled D shifted z80 withcolor RoyalBlue withpen pencircle scaled 1.0;

x81 = .3406062193u;
y81 = 1.0184380392u;
draw fullcircle scaled D shifted z81 withcolor RoyalBlue withpen pencircle scaled 1.0;

x82 = .3735509561u;
y82 = 1.0625023839u;
draw fullcircle scaled D shifted z82 withcolor RoyalBlue withpen pencircle scaled 1.0;

x83 = .4011196204u;
y83 = 1.1105983728u;
draw fullcircle scaled D shifted z83 withcolor RoyalBlue withpen pencircle scaled 1.0;

x84 = .4256985658u;
y84 = 1.1561216126u;
draw fullcircle scaled D shifted z84 withcolor RoyalBlue withpen pencircle scaled 1.0;

x85 = .3246138680u;
y85 = 1.0472130062u;
draw fullcircle scaled D shifted z85 withcolor RoyalBlue withpen pencircle scaled 1.0;

x86 = .3508171938u;
y86 = 1.0967265447u;
draw fullcircle scaled D shifted z86 withcolor RoyalBlue withpen pencircle scaled 1.0;

x87 = .3723548004u;
y87 = 1.1396527353u;
draw fullcircle scaled D shifted z87 withcolor RoyalBlue withpen pencircle scaled 1.0;

x88 = .3798961620u;
y88 = 1.1754340114u;
draw fullcircle scaled D shifted z88 withcolor RoyalBlue withpen pencircle scaled 1.0;

x89 = .3056893420u;
y89 = 1.0840608762u;
draw fullcircle scaled D shifted z89 withcolor RoyalBlue withpen pencircle scaled 1.0;

x90 = .3195129567u;
y90 = 1.1360122267u;
draw fullcircle scaled D shifted z90 withcolor RoyalBlue withpen pencircle scaled 1.0;

x91 = .3266222792u;
y91 = 1.1726221110u;
draw fullcircle scaled D shifted z91 withcolor RoyalBlue withpen pencircle scaled 1.0;

x92 = .2089220351u;
y92 = 1.1040138520u;
draw fullcircle scaled D shifted z92 withcolor RoyalBlue withpen pencircle scaled 1.0;

x93 = .2471931987u;
y93 = 1.1275297971u;
draw fullcircle scaled D shifted z93 withcolor RoyalBlue withpen pencircle scaled 1.0;

x94 = .2751812924u;
y94 = 1.1227148179u;
draw fullcircle scaled D shifted z94 withcolor RoyalBlue withpen pencircle scaled 1.0;

x95 = .2787893569u;
y95 = 1.1722745455u;
draw fullcircle scaled D shifted z95 withcolor RoyalBlue withpen pencircle scaled 1.0;

x96 = .1562217028u;
y96 = 1.1035126334u;
draw fullcircle scaled D shifted z96 withcolor RoyalBlue withpen pencircle scaled 1.0;

x97 = .2016957702u;
y97 = 1.1381189544u;
draw fullcircle scaled D shifted z97 withcolor RoyalBlue withpen pencircle scaled 1.0;

x98 = .1879599508u;
y98 = 1.1861760767u;
draw fullcircle scaled D shifted z98 withcolor RoyalBlue withpen pencircle scaled 1.0;

x99 = .2360015261u;
y99 = 1.1780278944u;
draw fullcircle scaled D shifted z99 withcolor RoyalBlue withpen pencircle scaled 1.0;

x100 = .1556386735u;
y100 = 1.1469878990u;
draw fullcircle scaled D shifted z100 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.5997441687u,        0.2796007626u)--(        0.6250000000u,        0.2650000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5997441687u,        0.2796007626u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5997441687u,        0.2796007626u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6250000000u,        0.2650000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6250000000u,        0.2650000000u) withcolor Yellow;
draw (        0.5997441687u,        0.2796007626u)--(        0.5744890304u,        0.2942017773u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5890920185u,        0.3194557745u)--(        0.5744890304u,        0.2942017773u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5997441687u,        0.2796007626u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5997441687u,        0.2796007626u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5890920185u,        0.3194557745u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5890920185u,        0.3194557745u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5744890304u,        0.2942017773u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5744890304u,        0.2942017773u) withcolor Green;
draw (        0.5890920185u,        0.3194557745u)--(        0.5717814195u,        0.3429367200u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5974000058u,        0.3568902026u)--(        0.5717814195u,        0.3429367200u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5890920185u,        0.3194557745u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5890920185u,        0.3194557745u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5974000058u,        0.3568902026u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5974000058u,        0.3568902026u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5717814195u,        0.3429367200u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5717814195u,        0.3429367200u) withcolor Green;
draw (        0.5890920185u,        0.3194557745u)--(        0.6147106048u,        0.3334092571u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5974000058u,        0.3568902026u)--(        0.6147106048u,        0.3334092571u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5890920185u,        0.3194557745u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5890920185u,        0.3194557745u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5974000058u,        0.3568902026u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5974000058u,        0.3568902026u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6147106048u,        0.3334092571u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6147106048u,        0.3334092571u) withcolor Green;
draw (        0.5974000058u,        0.3568902026u)--(        0.6246728813u,        0.3672440378u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6327026846u,        0.3952892494u)--(        0.6246728813u,        0.3672440378u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5974000058u,        0.3568902026u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5974000058u,        0.3568902026u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6327026846u,        0.3952892494u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6327026846u,        0.3952892494u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6246728813u,        0.3672440378u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6246728813u,        0.3672440378u) withcolor Green;
draw (        0.5974000058u,        0.3568902026u)--(        0.5703934964u,        0.3679201657u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5865417564u,        0.3922151310u)--(        0.5703934964u,        0.3679201657u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5974000058u,        0.3568902026u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5974000058u,        0.3568902026u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5865417564u,        0.3922151310u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5865417564u,        0.3922151310u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5703934964u,        0.3679201657u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5703934964u,        0.3679201657u) withcolor Green;
draw (        0.6327026846u,        0.3952892494u)--(        0.6600000000u,        0.3850000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6327026846u,        0.3952892494u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6327026846u,        0.3952892494u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6600000000u,        0.3850000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6600000000u,        0.3850000000u) withcolor Yellow;
draw (        0.6327026846u,        0.3952892494u)--(        0.6434819499u,        0.4223967959u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6668063752u,        0.4399177212u)--(        0.6434819499u,        0.4223967959u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6327026846u,        0.3952892494u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6327026846u,        0.3952892494u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6668063752u,        0.4399177212u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6668063752u,        0.4399177212u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6434819499u,        0.4223967959u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6434819499u,        0.4223967959u) withcolor MidnightBlue;
draw (        0.6327026846u,        0.3952892494u)--(        0.6560271099u,        0.4128101747u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6668063752u,        0.4399177212u)--(        0.6560271099u,        0.4128101747u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6327026846u,        0.3952892494u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6327026846u,        0.3952892494u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6668063752u,        0.4399177212u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6668063752u,        0.4399177212u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6560271099u,        0.4128101747u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6560271099u,        0.4128101747u) withcolor Green;
draw (        0.6327026846u,        0.3952892494u)--(        0.6084411160u,        0.4114876422u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5865417564u,        0.3922151310u)--(        0.6084411160u,        0.4114876422u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6327026846u,        0.3952892494u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6327026846u,        0.3952892494u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5865417564u,        0.3922151310u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5865417564u,        0.3922151310u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6084411160u,        0.4114876422u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6084411160u,        0.4114876422u) withcolor MidnightBlue;
draw (        0.6327026846u,        0.3952892494u)--(        0.6084411483u,        0.4114876905u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6192201461u,        0.4385953434u)--(        0.6084411483u,        0.4114876905u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6327026846u,        0.3952892494u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6327026846u,        0.3952892494u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6192201461u,        0.4385953434u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6192201461u,        0.4385953434u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6084411483u,        0.4114876905u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6084411483u,        0.4114876905u) withcolor MidnightBlue;
draw (        0.6327026846u,        0.3952892494u)--(        0.6434816824u,        0.4223969023u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6192201461u,        0.4385953434u)--(        0.6434816824u,        0.4223969023u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6327026846u,        0.3952892494u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6327026846u,        0.3952892494u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6192201461u,        0.4385953434u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6192201461u,        0.4385953434u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6434816824u,        0.4223969023u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6434816824u,        0.4223969023u) withcolor MidnightBlue;
draw (        0.6668063752u,        0.4399177212u)--(        0.6658085151u,        0.4690727536u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6931424786u,        0.4792642245u)--(        0.6658085151u,        0.4690727536u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6668063752u,        0.4399177212u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6668063752u,        0.4399177212u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6931424786u,        0.4792642245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6931424786u,        0.4792642245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6658085151u,        0.4690727536u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6658085151u,        0.4690727536u) withcolor MidnightBlue;
draw (        0.6668063752u,        0.4399177212u)--(        0.6941403387u,        0.4501091921u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6931424786u,        0.4792642245u)--(        0.6941403387u,        0.4501091921u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6668063752u,        0.4399177212u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6668063752u,        0.4399177212u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6931424786u,        0.4792642245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6931424786u,        0.4792642245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6941403387u,        0.4501091921u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6941403387u,        0.4501091921u) withcolor Green;
draw (        0.6668063752u,        0.4399177212u)--(        0.6434817707u,        0.4223970345u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6192201461u,        0.4385953434u)--(        0.6434817707u,        0.4223970345u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6668063752u,        0.4399177212u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6668063752u,        0.4399177212u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6192201461u,        0.4385953434u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6192201461u,        0.4385953434u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6434817707u,        0.4223970345u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6434817707u,        0.4223970345u) withcolor MidnightBlue;
draw (        0.6668063752u,        0.4399177212u)--(        0.6658084187u,        0.4690727503u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6403541248u,        0.4833237349u)--(        0.6658084187u,        0.4690727503u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6668063752u,        0.4399177212u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6668063752u,        0.4399177212u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6403541248u,        0.4833237349u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6403541248u,        0.4833237349u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6658084187u,        0.4690727503u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6658084187u,        0.4690727503u) withcolor MidnightBlue;
draw (        0.6931424786u,        0.4792642245u)--(        0.7030519057u,        0.5067016969u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7305543567u,        0.5164293290u)--(        0.7030519057u,        0.5067016969u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6931424786u,        0.4792642245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6931424786u,        0.4792642245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7305543567u,        0.5164293290u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7305543567u,        0.5164293290u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7030519057u,        0.5067016969u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7030519057u,        0.5067016969u) withcolor MidnightBlue;
draw (        0.6931424786u,        0.4792642245u)--(        0.7206449296u,        0.4889918566u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7305543567u,        0.5164293290u)--(        0.7206449296u,        0.4889918566u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6931424786u,        0.4792642245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6931424786u,        0.4792642245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7305543567u,        0.5164293290u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7305543567u,        0.5164293290u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7206449296u,        0.4889918566u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7206449296u,        0.4889918566u) withcolor Green;
draw (        0.6931424786u,        0.4792642245u)--(        0.6658084773u,        0.4690728550u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6403541248u,        0.4833237349u)--(        0.6658084773u,        0.4690728550u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6931424786u,        0.4792642245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6931424786u,        0.4792642245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6403541248u,        0.4833237349u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6403541248u,        0.4833237349u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6658084773u,        0.4690728550u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6658084773u,        0.4690728550u) withcolor MidnightBlue;
draw (        0.6931424786u,        0.4792642245u)--(        0.7030516641u,        0.5067017841u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6764840238u,        0.5187505238u)--(        0.7030516641u,        0.5067017841u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6931424786u,        0.4792642245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6931424786u,        0.4792642245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6764840238u,        0.5187505238u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6764840238u,        0.5187505238u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7030516641u,        0.5067017841u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7030516641u,        0.5067017841u) withcolor MidnightBlue;
draw (        0.7305543567u,        0.5164293290u)--(        0.7588388093u,        0.5235708538u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7452099634u,        0.5493636085u)--(        0.7588388093u,        0.5235708538u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7305543567u,        0.5164293290u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7305543567u,        0.5164293290u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7452099634u,        0.5493636085u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7452099634u,        0.5493636085u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7588388093u,        0.5235708538u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7588388093u,        0.5235708538u) withcolor Green;
draw (        0.7305543567u,        0.5164293290u)--(        0.7030517825u,        0.5067020452u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6764840238u,        0.5187505238u)--(        0.7030517825u,        0.5067020452u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7305543567u,        0.5164293290u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7305543567u,        0.5164293290u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6764840238u,        0.5187505238u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6764840238u,        0.5187505238u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7030517825u,        0.5067020452u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7030517825u,        0.5067020452u) withcolor MidnightBlue;
draw (        0.7305543567u,        0.5164293290u)--(        0.7039865980u,        0.5284778075u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6764840238u,        0.5187505238u)--(        0.7039865980u,        0.5284778075u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7305543567u,        0.5164293290u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7305543567u,        0.5164293290u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6764840238u,        0.5187505238u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6764840238u,        0.5187505238u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7039865980u,        0.5284778075u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7039865980u,        0.5284778075u) withcolor MidnightBlue;
draw (        0.7305543567u,        0.5164293290u)--(        0.7039866270u,        0.5284778716u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6980126540u,        0.5570317379u)--(        0.7039866270u,        0.5284778716u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7305543567u,        0.5164293290u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7305543567u,        0.5164293290u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6980126540u,        0.5570317379u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6980126540u,        0.5570317379u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7039866270u,        0.5284778716u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7039866270u,        0.5284778716u) withcolor MidnightBlue;
draw (        0.7452099634u,        0.5493636085u)--(        0.7367025377u,        0.5772676467u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7608512575u,        0.5936337999u)--(        0.7367025377u,        0.5772676467u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7452099634u,        0.5493636085u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7452099634u,        0.5493636085u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7608512575u,        0.5936337999u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7608512575u,        0.5936337999u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7367025377u,        0.5772676467u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7367025377u,        0.5772676467u) withcolor MidnightBlue;
draw (        0.7452099634u,        0.5493636085u)--(        0.7242919552u,        0.5696970434u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6980126540u,        0.5570317379u)--(        0.7242919552u,        0.5696970434u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7452099634u,        0.5493636085u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7452099634u,        0.5493636085u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6980126540u,        0.5570317379u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6980126540u,        0.5570317379u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7242919552u,        0.5696970434u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7242919552u,        0.5696970434u) withcolor MidnightBlue;
draw (        0.7452099634u,        0.5493636085u)--(        0.7242920560u,        0.5696971471u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7157844052u,        0.5976011166u)--(        0.7242920560u,        0.5696971471u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7452099634u,        0.5493636085u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7452099634u,        0.5493636085u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7157844052u,        0.5976011166u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7157844052u,        0.5976011166u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7242920560u,        0.5696971471u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7242920560u,        0.5696971471u) withcolor MidnightBlue;
draw (        0.7452099634u,        0.5493636085u)--(        0.7367023125u,        0.5772675781u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7157844052u,        0.5976011166u)--(        0.7367023125u,        0.5772675781u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7452099634u,        0.5493636085u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7452099634u,        0.5493636085u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7157844052u,        0.5976011166u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7157844052u,        0.5976011166u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7367023125u,        0.5772675781u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7367023125u,        0.5772675781u) withcolor MidnightBlue;
draw (        0.7608512575u,        0.5936337999u)--(        0.7850000000u,        0.6100000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7608512575u,        0.5936337999u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7608512575u,        0.5936337999u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7850000000u,        0.6100000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7850000000u,        0.6100000000u) withcolor Yellow;
draw (        0.7608512575u,        0.5936337999u)--(        0.7367024726u,        0.5772677428u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7157844052u,        0.5976011166u)--(        0.7367024726u,        0.5772677428u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7608512575u,        0.5936337999u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7608512575u,        0.5936337999u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7157844052u,        0.5976011166u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7157844052u,        0.5976011166u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7367024726u,        0.5772677428u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7367024726u,        0.5772677428u) withcolor MidnightBlue;
draw (        0.5865417564u,        0.3922151310u)--(        0.5973208694u,        0.4193227380u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6192201461u,        0.4385953434u)--(        0.5973208694u,        0.4193227380u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5865417564u,        0.3922151310u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5865417564u,        0.3922151310u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6192201461u,        0.4385953434u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6192201461u,        0.4385953434u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5973208694u,        0.4193227380u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5973208694u,        0.4193227380u) withcolor MidnightBlue;
draw (        0.5865417564u,        0.3922151310u)--(        0.6084410331u,        0.4114877363u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6192201461u,        0.4385953434u)--(        0.6084410331u,        0.4114877363u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5865417564u,        0.3922151310u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5865417564u,        0.3922151310u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6192201461u,        0.4385953434u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6192201461u,        0.4385953434u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6084410331u,        0.4114877363u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6084410331u,        0.4114877363u) withcolor MidnightBlue;
draw (        0.5865417564u,        0.3922151310u)--(        0.5700000124u,        0.4162439125u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5807790418u,        0.4433515528u)--(        0.5700000124u,        0.4162439125u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5865417564u,        0.3922151310u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5865417564u,        0.3922151310u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5807790418u,        0.4433515528u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5807790418u,        0.4433515528u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5700000124u,        0.4162439125u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5700000124u,        0.4162439125u) withcolor Green;
draw (        0.5865417564u,        0.3922151310u)--(        0.5973207858u,        0.4193227713u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5807790418u,        0.4433515528u)--(        0.5973207858u,        0.4193227713u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5865417564u,        0.3922151310u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5865417564u,        0.3922151310u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5807790418u,        0.4433515528u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5807790418u,        0.4433515528u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5973207858u,        0.4193227713u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5973207858u,        0.4193227713u) withcolor MidnightBlue;
draw (        0.6192201461u,        0.4385953434u)--(        0.6158036214u,        0.4675666913u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6403541248u,        0.4833237349u)--(        0.6158036214u,        0.4675666913u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6192201461u,        0.4385953434u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6192201461u,        0.4385953434u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6403541248u,        0.4833237349u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6403541248u,        0.4833237349u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6158036214u,        0.4675666913u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6158036214u,        0.4675666913u) withcolor MidnightBlue;
draw (        0.6192201461u,        0.4385953434u)--(        0.5973208197u,        0.4193227946u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5807790418u,        0.4433515528u)--(        0.5973208197u,        0.4193227946u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6192201461u,        0.4385953434u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6192201461u,        0.4385953434u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5807790418u,        0.4433515528u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5807790418u,        0.4433515528u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5973208197u,        0.4193227946u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5973208197u,        0.4193227946u) withcolor MidnightBlue;
draw (        0.6192201461u,        0.4385953434u)--(        0.6158033680u,        0.4675666614u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5940369321u,        0.4869891718u)--(        0.6158033680u,        0.4675666614u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6192201461u,        0.4385953434u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6192201461u,        0.4385953434u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5940369321u,        0.4869891718u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5940369321u,        0.4869891718u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6158033680u,        0.4675666614u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6158033680u,        0.4675666614u) withcolor MidnightBlue;
draw (        0.6403541248u,        0.4833237349u)--(        0.6482514633u,        0.5114065353u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6764840238u,        0.5187505238u)--(        0.6482514633u,        0.5114065353u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6403541248u,        0.4833237349u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6403541248u,        0.4833237349u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6764840238u,        0.5187505238u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6764840238u,        0.5187505238u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6482514633u,        0.5114065353u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6482514633u,        0.5114065353u) withcolor MidnightBlue;
draw (        0.6403541248u,        0.4833237349u)--(        0.6158035266u,        0.4675668391u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5940369321u,        0.4869891718u)--(        0.6158035266u,        0.4675668391u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6403541248u,        0.4833237349u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6403541248u,        0.4833237349u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5940369321u,        0.4869891718u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5940369321u,        0.4869891718u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6158035266u,        0.4675668391u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6158035266u,        0.4675668391u) withcolor MidnightBlue;
draw (        0.6403541248u,        0.4833237349u)--(        0.6482513099u,        0.5114065785u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6252797313u,        0.5293876300u)--(        0.6482513099u,        0.5114065785u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6403541248u,        0.4833237349u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6403541248u,        0.4833237349u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6252797313u,        0.5293876300u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6252797313u,        0.5293876300u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6482513099u,        0.5114065785u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6482513099u,        0.5114065785u) withcolor MidnightBlue;
draw (        0.6764840238u,        0.5187505238u)--(        0.6705100989u,        0.5473044002u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6980126540u,        0.5570317379u)--(        0.6705100989u,        0.5473044002u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6764840238u,        0.5187505238u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6764840238u,        0.5187505238u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6980126540u,        0.5570317379u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6980126540u,        0.5570317379u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6705100989u,        0.5473044002u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6705100989u,        0.5473044002u) withcolor MidnightBlue;
draw (        0.6764840238u,        0.5187505238u)--(        0.7039865789u,        0.5284778615u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6980126540u,        0.5570317379u)--(        0.7039865789u,        0.5284778615u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6764840238u,        0.5187505238u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6764840238u,        0.5187505238u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6980126540u,        0.5570317379u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6980126540u,        0.5570317379u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7039865789u,        0.5284778615u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7039865789u,        0.5284778615u) withcolor MidnightBlue;
draw (        0.6764840238u,        0.5187505238u)--(        0.6482514166u,        0.5114067148u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6252797313u,        0.5293876300u)--(        0.6482514166u,        0.5114067148u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6764840238u,        0.5187505238u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6764840238u,        0.5187505238u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6252797313u,        0.5293876300u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6252797313u,        0.5293876300u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6482514166u,        0.5114067148u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6482514166u,        0.5114067148u) withcolor MidnightBlue;
draw (        0.6764840238u,        0.5187505238u)--(        0.6535123385u,        0.5367314390u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6252797313u,        0.5293876300u)--(        0.6535123385u,        0.5367314390u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6764840238u,        0.5187505238u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6764840238u,        0.5187505238u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6252797313u,        0.5293876300u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6252797313u,        0.5293876300u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6535123385u,        0.5367314390u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6535123385u,        0.5367314390u) withcolor MidnightBlue;
draw (        0.6764840238u,        0.5187505238u)--(        0.6535124365u,        0.5367315642u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6475384636u,        0.5652854306u)--(        0.6535124365u,        0.5367315642u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6764840238u,        0.5187505238u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6764840238u,        0.5187505238u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6535124365u,        0.5367315642u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6535124365u,        0.5367315642u) withcolor MidnightBlue;
draw (        0.6764840238u,        0.5187505238u)--(        0.6705100509u,        0.5473043902u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6475384636u,        0.5652854306u)--(        0.6705100509u,        0.5473043902u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6764840238u,        0.5187505238u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6764840238u,        0.5187505238u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6705100509u,        0.5473043902u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6705100509u,        0.5473043902u) withcolor MidnightBlue;
draw (        0.6980126540u,        0.5570317379u)--(        0.7242919246u,        0.5696971070u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7157844052u,        0.5976011166u)--(        0.7242919246u,        0.5696971070u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6980126540u,        0.5570317379u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6980126540u,        0.5570317379u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7157844052u,        0.5976011166u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7157844052u,        0.5976011166u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7242919246u,        0.5696971070u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7242919246u,        0.5696971070u) withcolor MidnightBlue;
draw (        0.6980126540u,        0.5570317379u)--(        0.6705100864u,        0.5473044355u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6475384636u,        0.5652854306u)--(        0.6705100864u,        0.5473044355u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6980126540u,        0.5570317379u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6980126540u,        0.5570317379u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6705100864u,        0.5473044355u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6705100864u,        0.5473044355u) withcolor MidnightBlue;
draw (        0.6980126540u,        0.5570317379u)--(        0.6750410312u,        0.5750127330u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6475384636u,        0.5652854306u)--(        0.6750410312u,        0.5750127330u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6980126540u,        0.5570317379u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6980126540u,        0.5570317379u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6750410312u,        0.5750127330u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6750410312u,        0.5750127330u) withcolor MidnightBlue;
draw (        0.6980126540u,        0.5570317379u)--(        0.6750410612u,        0.5750127714u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6755140958u,        0.6041810397u)--(        0.6750410612u,        0.5750127714u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6980126540u,        0.5570317379u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6980126540u,        0.5570317379u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6755140958u,        0.6041810397u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6755140958u,        0.6041810397u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6750410612u,        0.5750127714u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6750410612u,        0.5750127714u) withcolor MidnightBlue;
draw (        0.7157844052u,        0.5976011166u)--(        0.6990115859u,        0.6214691742u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6755140958u,        0.6041810397u)--(        0.6990115859u,        0.6214691742u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7157844052u,        0.5976011166u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7157844052u,        0.5976011166u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6755140958u,        0.6041810397u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6755140958u,        0.6041810397u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6990115859u,        0.6214691742u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6990115859u,        0.6214691742u) withcolor MidnightBlue;
draw (        0.7157844052u,        0.5976011166u)--(        0.6990118774u,        0.6214693790u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6921305001u,        0.6498182465u)--(        0.6990118774u,        0.6214693790u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7157844052u,        0.5976011166u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7157844052u,        0.5976011166u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6921305001u,        0.6498182465u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6921305001u,        0.6498182465u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6990118774u,        0.6214693790u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6990118774u,        0.6214693790u) withcolor MidnightBlue;
draw (        0.7362248211u,        0.6226591835u)--(        0.7212240570u,        0.6476789499u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6921305001u,        0.6498182465u)--(        0.7212240570u,        0.6476789499u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7362248211u,        0.6226591835u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7362248211u,        0.6226591835u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6921305001u,        0.6498182465u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6921305001u,        0.6498182465u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7212240570u,        0.6476789499u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7212240570u,        0.6476789499u) withcolor Green;
draw (        0.5807790418u,        0.4433515528u)--(        0.5700000462u,        0.4704592065u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5940369321u,        0.4869891718u)--(        0.5700000462u,        0.4704592065u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5807790418u,        0.4433515528u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5807790418u,        0.4433515528u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5940369321u,        0.4869891718u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5940369321u,        0.4869891718u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5700000462u,        0.4704592065u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5700000462u,        0.4704592065u) withcolor Green;
draw (        0.5940369321u,        0.4869891718u)--(        0.5995527257u,        0.5156350730u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6252797313u,        0.5293876300u)--(        0.5995527257u,        0.5156350730u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5940369321u,        0.4869891718u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5940369321u,        0.4869891718u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6252797313u,        0.5293876300u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6252797313u,        0.5293876300u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5995527257u,        0.5156350730u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5995527257u,        0.5156350730u) withcolor MidnightBlue;
draw (        0.5940369321u,        0.4869891718u)--(        0.5700000909u,        0.5035192020u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5755157707u,        0.5321651252u)--(        0.5700000909u,        0.5035192020u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5940369321u,        0.4869891718u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5940369321u,        0.4869891718u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5755157707u,        0.5321651252u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5755157707u,        0.5321651252u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5700000909u,        0.5035192020u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5700000909u,        0.5035192020u) withcolor Green;
draw (        0.5940369321u,        0.4869891718u)--(        0.5995526120u,        0.5156350949u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5755157707u,        0.5321651252u)--(        0.5995526120u,        0.5156350949u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5940369321u,        0.4869891718u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5940369321u,        0.4869891718u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5755157707u,        0.5321651252u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5755157707u,        0.5321651252u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5995526120u,        0.5156350949u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5995526120u,        0.5156350949u) withcolor MidnightBlue;
draw (        0.6252797313u,        0.5293876300u)--(        0.6193058823u,        0.5579415223u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6475384636u,        0.5652854306u)--(        0.6193058823u,        0.5579415223u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6252797313u,        0.5293876300u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6252797313u,        0.5293876300u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6193058823u,        0.5579415223u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6193058823u,        0.5579415223u) withcolor MidnightBlue;
draw (        0.6252797313u,        0.5293876300u)--(        0.6535123126u,        0.5367315383u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6475384636u,        0.5652854306u)--(        0.6535123126u,        0.5367315383u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6252797313u,        0.5293876300u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6252797313u,        0.5293876300u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6535123126u,        0.5367315383u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6535123126u,        0.5367315383u) withcolor MidnightBlue;
draw (        0.6252797313u,        0.5293876300u)--(        0.5995526694u,        0.5156351784u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5755157707u,        0.5321651252u)--(        0.5995526694u,        0.5156351784u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6252797313u,        0.5293876300u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6252797313u,        0.5293876300u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5755157707u,        0.5321651252u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5755157707u,        0.5321651252u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5995526694u,        0.5156351784u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5995526694u,        0.5156351784u) withcolor MidnightBlue;
draw (        0.6252797313u,        0.5293876300u)--(        0.6012428326u,        0.5459175767u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5755157707u,        0.5321651252u)--(        0.6012428326u,        0.5459175767u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6252797313u,        0.5293876300u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6252797313u,        0.5293876300u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5755157707u,        0.5321651252u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5755157707u,        0.5321651252u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6012428326u,        0.5459175767u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6012428326u,        0.5459175767u) withcolor MidnightBlue;
draw (        0.6252797313u,        0.5293876300u)--(        0.6012429142u,        0.5459176953u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5952690094u,        0.5744715760u)--(        0.6012429142u,        0.5459176953u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6252797313u,        0.5293876300u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6252797313u,        0.5293876300u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5952690094u,        0.5744715760u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5952690094u,        0.5744715760u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6012429142u,        0.5459176953u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6012429142u,        0.5459176953u) withcolor MidnightBlue;
draw (        0.6252797313u,        0.5293876300u)--(        0.6193058266u,        0.5579415106u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5952690094u,        0.5744715760u)--(        0.6193058266u,        0.5579415106u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6252797313u,        0.5293876300u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6252797313u,        0.5293876300u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5952690094u,        0.5744715760u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5952690094u,        0.5744715760u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6193058266u,        0.5579415106u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6193058266u,        0.5579415106u) withcolor MidnightBlue;
draw (        0.6475384636u,        0.5652854306u)--(        0.6480115420u,        0.5944536982u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6755140958u,        0.6041810397u)--(        0.6480115420u,        0.5944536982u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6755140958u,        0.6041810397u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6755140958u,        0.6041810397u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6480115420u,        0.5944536982u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6480115420u,        0.5944536982u) withcolor MidnightBlue;
draw (        0.6475384636u,        0.5652854306u)--(        0.6750410174u,        0.5750127721u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6755140958u,        0.6041810397u)--(        0.6750410174u,        0.5750127721u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6755140958u,        0.6041810397u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6755140958u,        0.6041810397u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6750410174u,        0.5750127721u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6750410174u,        0.5750127721u) withcolor MidnightBlue;
draw (        0.6475384636u,        0.5652854306u)--(        0.6193058692u,        0.5579415726u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5952690094u,        0.5744715760u)--(        0.6193058692u,        0.5579415726u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5952690094u,        0.5744715760u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5952690094u,        0.5744715760u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6193058692u,        0.5579415726u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6193058692u,        0.5579415726u) withcolor MidnightBlue;
draw (        0.6475384636u,        0.5652854306u)--(        0.6235016039u,        0.5818154339u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5952690094u,        0.5744715760u)--(        0.6235016039u,        0.5818154339u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5952690094u,        0.5744715760u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5952690094u,        0.5744715760u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6235016039u,        0.5818154339u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6235016039u,        0.5818154339u) withcolor MidnightBlue;
draw (        0.6475384636u,        0.5652854306u)--(        0.6235016465u,        0.5818154959u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6239746084u,        0.6109837654u)--(        0.6235016465u,        0.5818154959u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6239746084u,        0.6109837654u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6239746084u,        0.6109837654u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6235016465u,        0.5818154959u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6235016465u,        0.5818154959u) withcolor MidnightBlue;
draw (        0.6475384636u,        0.5652854306u)--(        0.6480114256u,        0.5944537001u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6239746084u,        0.6109837654u)--(        0.6480114256u,        0.5944537001u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6475384636u,        0.5652854306u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6239746084u,        0.6109837654u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6239746084u,        0.6109837654u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6480114256u,        0.5944537001u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6480114256u,        0.5944537001u) withcolor MidnightBlue;
draw (        0.6755140958u,        0.6041810397u)--(        0.6990115022u,        0.6214692879u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6921305001u,        0.6498182465u)--(        0.6990115022u,        0.6214692879u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6755140958u,        0.6041810397u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6755140958u,        0.6041810397u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6921305001u,        0.6498182465u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6921305001u,        0.6498182465u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6990115022u,        0.6214692879u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6990115022u,        0.6214692879u) withcolor MidnightBlue;
draw (        0.6755140958u,        0.6041810397u)--(        0.6480115020u,        0.5944538112u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6239746084u,        0.6109837654u)--(        0.6480115020u,        0.5944538112u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6755140958u,        0.6041810397u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6755140958u,        0.6041810397u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6239746084u,        0.6109837654u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6239746084u,        0.6109837654u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6480115020u,        0.5944538112u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6480115020u,        0.5944538112u) withcolor MidnightBlue;
draw (        0.6921305001u,        0.6498182465u)--(        0.6852488994u,        0.6781670596u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6561892216u,        0.6807257247u)--(        0.6852488994u,        0.6781670596u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6921305001u,        0.6498182465u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6921305001u,        0.6498182465u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6561892216u,        0.6807257247u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6561892216u,        0.6807257247u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6852488994u,        0.6781670596u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6852488994u,        0.6781670596u) withcolor Green;
draw (        0.5755157707u,        0.5321651252u)--(        0.6012427838u,        0.5459176680u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5952690094u,        0.5744715760u)--(        0.6012427838u,        0.5459176680u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5755157707u,        0.5321651252u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5755157707u,        0.5321651252u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5952690094u,        0.5744715760u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5952690094u,        0.5744715760u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6012427838u,        0.5459176680u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6012427838u,        0.5459176680u) withcolor MidnightBlue;
draw (        0.5755157707u,        0.5321651252u)--(        0.5508664442u,        0.5477671232u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5524199618u,        0.5768978324u)--(        0.5508664442u,        0.5477671232u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5755157707u,        0.5321651252u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5755157707u,        0.5321651252u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5524199618u,        0.5768978324u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5524199618u,        0.5768978324u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5508664442u,        0.5477671232u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5508664442u,        0.5477671232u) withcolor Green;
draw (        0.5952690094u,        0.5744715760u)--(        0.6235015875u,        0.5818154968u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6239746084u,        0.6109837654u)--(        0.6235015875u,        0.5818154968u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5952690094u,        0.5744715760u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5952690094u,        0.5744715760u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6239746084u,        0.6109837654u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6239746084u,        0.6109837654u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6235015875u,        0.5818154968u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6235015875u,        0.5818154968u) withcolor MidnightBlue;
draw (        0.5952690094u,        0.5744715760u)--(        0.5749616768u,        0.5954149252u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5524199618u,        0.5768978324u)--(        0.5749616768u,        0.5954149252u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5952690094u,        0.5744715760u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5952690094u,        0.5744715760u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5524199618u,        0.5768978324u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5524199618u,        0.5768978324u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5749616768u,        0.5954149252u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5749616768u,        0.5954149252u) withcolor MidnightBlue;
draw (        0.5952690094u,        0.5744715760u)--(        0.5749618334u,        0.5954150771u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5778526313u,        0.6244435963u)--(        0.5749618334u,        0.5954150771u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5952690094u,        0.5744715760u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5952690094u,        0.5744715760u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5778526313u,        0.6244435963u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5778526313u,        0.6244435963u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5749618334u,        0.5954150771u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5749618334u,        0.5954150771u) withcolor MidnightBlue;
draw (        0.6239746084u,        0.6109837654u)--(        0.6055500978u,        0.6336012154u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5778526313u,        0.6244435963u)--(        0.6055500978u,        0.6336012154u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6239746084u,        0.6109837654u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6239746084u,        0.6109837654u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5778526313u,        0.6244435963u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5778526313u,        0.6244435963u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6055500978u,        0.6336012154u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6055500978u,        0.6336012154u) withcolor MidnightBlue;
draw (        0.6239746084u,        0.6109837654u)--(        0.6055502328u,        0.6336013254u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6032885711u,        0.6626856257u)--(        0.6055502328u,        0.6336013254u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6239746084u,        0.6109837654u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6239746084u,        0.6109837654u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6032885711u,        0.6626856257u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6032885711u,        0.6626856257u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6055502328u,        0.6336013254u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6055502328u,        0.6336013254u) withcolor MidnightBlue;
draw (        0.6501069626u,        0.6407568395u)--(        0.6323545189u,        0.6639055410u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6561892216u,        0.6807257247u)--(        0.6323545189u,        0.6639055410u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6501069626u,        0.6407568395u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6501069626u,        0.6407568395u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6561892216u,        0.6807257247u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6561892216u,        0.6807257247u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6323545189u,        0.6639055410u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6323545189u,        0.6639055410u) withcolor MidnightBlue;
draw (        0.6501069626u,        0.6407568395u)--(        0.6324325980u,        0.6639652103u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6032885711u,        0.6626856257u)--(        0.6324325980u,        0.6639652103u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6501069626u,        0.6407568395u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6501069626u,        0.6407568395u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6032885711u,        0.6626856257u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6032885711u,        0.6626856257u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6324325980u,        0.6639652103u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6324325980u,        0.6639652103u) withcolor MidnightBlue;
draw (        0.6561892216u,        0.6807257247u)--(        0.6270378457u,        0.6796262158u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6032885711u,        0.6626856257u)--(        0.6270378457u,        0.6796262158u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6561892216u,        0.6807257247u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6561892216u,        0.6807257247u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6032885711u,        0.6626856257u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6032885711u,        0.6626856257u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6270378457u,        0.6796262158u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6270378457u,        0.6796262158u) withcolor MidnightBlue;
draw (        0.6561892216u,        0.6807257247u)--(        0.6270378343u,        0.6796265181u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6140175721u,        0.7057317745u)--(        0.6270378343u,        0.6796265181u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6561892216u,        0.6807257247u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6561892216u,        0.6807257247u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6140175721u,        0.7057317745u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6140175721u,        0.7057317745u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6270378343u,        0.6796265181u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6270378343u,        0.6796265181u) withcolor MidnightBlue;
draw (        0.6561892216u,        0.6807257247u)--(        0.6431689594u,        0.7068309811u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6140175721u,        0.7057317745u)--(        0.6431689594u,        0.7068309811u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6561892216u,        0.6807257247u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6561892216u,        0.6807257247u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6140175721u,        0.7057317745u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6140175721u,        0.7057317745u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6431689594u,        0.7068309811u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6431689594u,        0.7068309811u) withcolor Green;
draw (        0.5524199618u,        0.5768978324u)--(        0.5749615269u,        0.5954151076u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5778526313u,        0.6244435963u)--(        0.5749615269u,        0.5954151076u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5524199618u,        0.5768978324u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5524199618u,        0.5768978324u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5778526313u,        0.6244435963u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5778526313u,        0.6244435963u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5749615269u,        0.5954151076u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5749615269u,        0.5954151076u) withcolor MidnightBlue;
draw (        0.5524199618u,        0.5768978324u)--(        0.5297574584u,        0.5952668983u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5352615775u,        0.6239150450u)--(        0.5297574584u,        0.5952668983u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5524199618u,        0.5768978324u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5524199618u,        0.5768978324u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5352615775u,        0.6239150450u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5352615775u,        0.6239150450u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5297574584u,        0.5952668983u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5297574584u,        0.5952668983u) withcolor Green;
draw (        0.5778526313u,        0.6244435963u)--(        0.5755911367u,        0.6535279096u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6032885711u,        0.6626856257u)--(        0.5755911367u,        0.6535279096u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5778526313u,        0.6244435963u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5778526313u,        0.6244435963u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6032885711u,        0.6626856257u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6032885711u,        0.6626856257u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5755911367u,        0.6535279096u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5755911367u,        0.6535279096u) withcolor MidnightBlue;
draw (        0.5778526313u,        0.6244435963u)--(        0.6055500657u,        0.6336013124u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6032885711u,        0.6626856257u)--(        0.6055500657u,        0.6336013124u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5778526313u,        0.6244435963u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5778526313u,        0.6244435963u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6032885711u,        0.6626856257u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6032885711u,        0.6626856257u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6055500657u,        0.6336013124u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6055500657u,        0.6336013124u) withcolor MidnightBlue;
draw (        0.5778526313u,        0.6244435963u)--(        0.5563097199u,        0.6441137414u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5352615775u,        0.6239150450u)--(        0.5563097199u,        0.6441137414u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5778526313u,        0.6244435963u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5778526313u,        0.6244435963u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5352615775u,        0.6239150450u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5352615775u,        0.6239150450u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5563097199u,        0.6441137414u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5563097199u,        0.6441137414u) withcolor MidnightBlue;
draw (        0.5778526313u,        0.6244435963u)--(        0.5563097645u,        0.6441137902u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5540480399u,        0.6731980856u)--(        0.5563097645u,        0.6441137902u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5778526313u,        0.6244435963u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5778526313u,        0.6244435963u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5540480399u,        0.6731980856u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5540480399u,        0.6731980856u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5563097645u,        0.6441137902u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5563097645u,        0.6441137902u) withcolor MidnightBlue;
draw (        0.5778526313u,        0.6244435963u)--(        0.5755909067u,        0.6535278917u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5540480399u,        0.6731980856u)--(        0.5755909067u,        0.6535278917u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5778526313u,        0.6244435963u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5778526313u,        0.6244435963u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5540480399u,        0.6731980856u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5540480399u,        0.6731980856u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5755909067u,        0.6535278917u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5755909067u,        0.6535278917u) withcolor MidnightBlue;
draw (        0.6032885711u,        0.6626856257u)--(        0.6270376837u,        0.6796264430u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6140175721u,        0.7057317745u)--(        0.6270376837u,        0.6796264430u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6032885711u,        0.6626856257u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6032885711u,        0.6626856257u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6140175721u,        0.7057317745u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6140175721u,        0.7057317745u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6270376837u,        0.6796264430u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6270376837u,        0.6796264430u) withcolor MidnightBlue;
draw (        0.6032885711u,        0.6626856257u)--(        0.5755910799u,        0.6535280814u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5540480399u,        0.6731980856u)--(        0.5755910799u,        0.6535280814u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6032885711u,        0.6626856257u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6032885711u,        0.6626856257u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5540480399u,        0.6731980856u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5540480399u,        0.6731980856u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5755910799u,        0.6535280814u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5755910799u,        0.6535280814u) withcolor MidnightBlue;
draw (        0.6032885711u,        0.6626856257u)--(        0.5817455311u,        0.6823556299u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5540480399u,        0.6731980856u)--(        0.5817455311u,        0.6823556299u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6032885711u,        0.6626856257u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6032885711u,        0.6626856257u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5540480399u,        0.6731980856u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5540480399u,        0.6731980856u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5817455311u,        0.6823556299u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5817455311u,        0.6823556299u) withcolor MidnightBlue;
draw (        0.6032885711u,        0.6626856257u)--(        0.5817455932u,        0.6823556980u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5710107301u,        0.7094808588u)--(        0.5817455932u,        0.6823556980u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6032885711u,        0.6626856257u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6032885711u,        0.6626856257u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5710107301u,        0.7094808588u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5710107301u,        0.7094808588u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5817455932u,        0.6823556980u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5817455932u,        0.6823556980u) withcolor MidnightBlue;
draw (        0.6140175721u,        0.7057317745u)--(        0.5942183914u,        0.7271561550u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5710107301u,        0.7094808588u)--(        0.5942183914u,        0.7271561550u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6140175721u,        0.7057317745u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6140175721u,        0.7057317745u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5710107301u,        0.7094808588u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5710107301u,        0.7094808588u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5942183914u,        0.7271561550u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5942183914u,        0.7271561550u) withcolor Green;
draw (        0.5352615775u,        0.6239150450u)--(        0.5329999378u,        0.6529993470u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5540480399u,        0.6731980856u)--(        0.5329999378u,        0.6529993470u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5352615775u,        0.6239150450u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5352615775u,        0.6239150450u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5540480399u,        0.6731980856u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5540480399u,        0.6731980856u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5329999378u,        0.6529993470u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5329999378u,        0.6529993470u) withcolor Green;
draw (        0.5352615775u,        0.6239150450u)--(        0.5563096795u,        0.6441137836u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5540480399u,        0.6731980856u)--(        0.5563096795u,        0.6441137836u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5352615775u,        0.6239150450u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5352615775u,        0.6239150450u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5540480399u,        0.6731980856u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5540480399u,        0.6731980856u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5563096795u,        0.6441137836u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5563096795u,        0.6441137836u) withcolor MidnightBlue;
draw (        0.5540480399u,        0.6731980856u)--(        0.5433132516u,        0.7003232761u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5710107301u,        0.7094808588u)--(        0.5433132516u,        0.7003232761u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5540480399u,        0.6731980856u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5540480399u,        0.6731980856u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5710107301u,        0.7094808588u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5710107301u,        0.7094808588u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5433132516u,        0.7003232761u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5433132516u,        0.7003232761u) withcolor MidnightBlue;
draw (        0.5540480399u,        0.6731980856u)--(        0.5817455184u,        0.6823556683u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5710107301u,        0.7094808588u)--(        0.5817455184u,        0.6823556683u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5540480399u,        0.6731980856u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5540480399u,        0.6731980856u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5710107301u,        0.7094808588u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5710107301u,        0.7094808588u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5817455184u,        0.6823556683u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5817455184u,        0.6823556683u) withcolor MidnightBlue;
draw (        0.5540480399u,        0.6731980856u)--(        0.5433129108u,        0.7003231412u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5198362864u,        0.7176396000u)--(        0.5433129108u,        0.7003231412u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5540480399u,        0.6731980856u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5540480399u,        0.6731980856u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5198362864u,        0.7176396000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5198362864u,        0.7176396000u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5433129108u,        0.7003231412u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5433129108u,        0.7003231412u) withcolor MidnightBlue;
draw (        0.5710107301u,        0.7094808588u)--(        0.5433131769u,        0.7003235020u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5198362864u,        0.7176396000u)--(        0.5433131769u,        0.7003235020u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5710107301u,        0.7094808588u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5710107301u,        0.7094808588u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5198362864u,        0.7176396000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5198362864u,        0.7176396000u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5433131769u,        0.7003235020u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5433131769u,        0.7003235020u) withcolor MidnightBlue;
draw (        0.5198362864u,        0.7176396000u)--(        0.4924127866u,        0.7275876297u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4825763447u,        0.7550513521u)--(        0.4924127866u,        0.7275876297u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5198362864u,        0.7176396000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5198362864u,        0.7176396000u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4825763447u,        0.7550513521u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4825763447u,        0.7550513521u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4924127866u,        0.7275876297u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4924127866u,        0.7275876297u) withcolor Green;
draw (        0.5198362864u,        0.7176396000u)--(        0.5099998445u,        0.7451033224u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4825763447u,        0.7550513521u)--(        0.5099998445u,        0.7451033224u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5198362864u,        0.7176396000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5198362864u,        0.7176396000u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4825763447u,        0.7550513521u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4825763447u,        0.7550513521u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5099998445u,        0.7451033224u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5099998445u,        0.7451033224u) withcolor Green;
draw (        0.4825763447u,        0.7550513521u)--(        0.5100000000u,        0.7650000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4825763447u,        0.7550513521u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4825763447u,        0.7550513521u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5100000000u,        0.7650000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5100000000u,        0.7650000000u) withcolor Yellow;
draw (        0.4825763447u,        0.7550513521u)--(        0.4534431980u,        0.7565584697u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4280920558u,        0.7709921534u)--(        0.4534431980u,        0.7565584697u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4825763447u,        0.7550513521u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4825763447u,        0.7550513521u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4280920558u,        0.7709921534u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4280920558u,        0.7709921534u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4534431980u,        0.7565584697u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4534431980u,        0.7565584697u) withcolor Green;
draw (        0.5421056577u,        0.8180848237u)--(        0.5500000000u,        0.7900000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5421056577u,        0.8180848237u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5421056577u,        0.8180848237u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5500000000u,        0.7900000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5500000000u,        0.7900000000u) withcolor Yellow;
draw (        0.5421056577u,        0.8180848237u)--(        0.5650000000u,        0.8000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5421056577u,        0.8180848237u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5421056577u,        0.8180848237u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5650000000u,        0.8000000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5650000000u,        0.8000000000u) withcolor Yellow;
draw (        0.5299513321u,        0.8687716044u)--(        0.5280903374u,        0.8978842878u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5533000523u,        0.9125635899u)--(        0.5280903374u,        0.8978842878u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5299513321u,        0.8687716044u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5299513321u,        0.8687716044u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5533000523u,        0.9125635899u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5533000523u,        0.9125635899u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5280903374u,        0.8978842878u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5280903374u,        0.8978842878u) withcolor MidnightBlue;
draw (        0.5299513321u,        0.8687716044u)--(        0.5134494447u,        0.8928277751u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4871977584u,        0.8801053301u)--(        0.5134494447u,        0.8928277751u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5299513321u,        0.8687716044u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5299513321u,        0.8687716044u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4871977584u,        0.8801053301u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4871977584u,        0.8801053301u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5134494447u,        0.8928277751u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5134494447u,        0.8928277751u) withcolor MidnightBlue;
draw (        0.5299513321u,        0.8687716044u)--(        0.5134559527u,        0.8928322382u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5115937696u,        0.9219448456u)--(        0.5134559527u,        0.8928322382u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5299513321u,        0.8687716044u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5299513321u,        0.8687716044u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5115937696u,        0.9219448456u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5115937696u,        0.9219448456u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5134559527u,        0.8928322382u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5134559527u,        0.8928322382u) withcolor MidnightBlue;
draw (        0.5299513321u,        0.8687716044u)--(        0.5280891489u,        0.8978842118u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5115937696u,        0.9219448456u)--(        0.5280891489u,        0.8978842118u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5299513321u,        0.8687716044u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5299513321u,        0.8687716044u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5115937696u,        0.9219448456u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5115937696u,        0.9219448456u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5280891489u,        0.8978842118u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5280891489u,        0.8978842118u) withcolor MidnightBlue;
draw (        0.5533000523u,        0.9125635899u)--(        0.5603409353u,        0.9408732633u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5852875073u,        0.9559954384u)--(        0.5603409353u,        0.9408732633u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5533000523u,        0.9125635899u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5533000523u,        0.9125635899u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5852875073u,        0.9559954384u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5852875073u,        0.9559954384u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5603409353u,        0.9408732633u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5603409353u,        0.9408732633u) withcolor MidnightBlue;
draw (        0.5533000523u,        0.9125635899u)--(        0.5280900299u,        0.8978848158u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5115937696u,        0.9219448456u)--(        0.5280900299u,        0.8978848158u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5533000523u,        0.9125635899u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5533000523u,        0.9125635899u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5115937696u,        0.9219448456u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5115937696u,        0.9219448456u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5280900299u,        0.8978848158u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5280900299u,        0.8978848158u) withcolor MidnightBlue;
draw (        0.5533000523u,        0.9125635899u)--(        0.5368037920u,        0.9366236197u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5115937696u,        0.9219448456u)--(        0.5368037920u,        0.9366236197u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5533000523u,        0.9125635899u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5533000523u,        0.9125635899u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5115937696u,        0.9219448456u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5115937696u,        0.9219448456u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5368037920u,        0.9366236197u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5368037920u,        0.9366236197u) withcolor MidnightBlue;
draw (        0.5533000523u,        0.9125635899u)--(        0.5368089815u,        0.9366271769u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5438378022u,        0.9649398475u)--(        0.5368089815u,        0.9366271769u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5533000523u,        0.9125635899u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5533000523u,        0.9125635899u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5438378022u,        0.9649398475u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5438378022u,        0.9649398475u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5368089815u,        0.9366271769u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5368089815u,        0.9366271769u) withcolor MidnightBlue;
draw (        0.5533000523u,        0.9125635899u)--(        0.5603288730u,        0.9408762605u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5438378022u,        0.9649398475u)--(        0.5603288730u,        0.9408762605u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5533000523u,        0.9125635899u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5533000523u,        0.9125635899u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5438378022u,        0.9649398475u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5438378022u,        0.9649398475u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5603288730u,        0.9408762605u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5603288730u,        0.9408762605u) withcolor MidnightBlue;
draw (        0.5852875073u,        0.9559954384u)--(        0.6007631161u,        0.9807243116u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6289565053u,        0.9882172714u)--(        0.6007631161u,        0.9807243116u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5852875073u,        0.9559954384u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5852875073u,        0.9559954384u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6289565053u,        0.9882172714u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6289565053u,        0.9882172714u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6007631161u,        0.9807243116u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6007631161u,        0.9807243116u) withcolor MidnightBlue;
draw (        0.5852875073u,        0.9559954384u)--(        0.6134808965u,        0.9634883982u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6289565053u,        0.9882172714u)--(        0.6134808965u,        0.9634883982u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5852875073u,        0.9559954384u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5852875073u,        0.9559954384u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6289565053u,        0.9882172714u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6289565053u,        0.9882172714u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6134808965u,        0.9634883982u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6134808965u,        0.9634883982u) withcolor Green;
draw (        0.5852875073u,        0.9559954384u)--(        0.5603361119u,        0.9408812230u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5438378022u,        0.9649398475u)--(        0.5603361119u,        0.9408812230u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5852875073u,        0.9559954384u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5852875073u,        0.9559954384u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5438378022u,        0.9649398475u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5438378022u,        0.9649398475u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5603361119u,        0.9408812230u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5603361119u,        0.9408812230u) withcolor MidnightBlue;
draw (        0.5852875073u,        0.9559954384u)--(        0.6007598673u,        0.9807263444u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5777054129u,        0.9986010128u)--(        0.6007598673u,        0.9807263444u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5852875073u,        0.9559954384u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5852875073u,        0.9559954384u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5777054129u,        0.9986010128u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5777054129u,        0.9986010128u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6007598673u,        0.9807263444u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6007598673u,        0.9807263444u) withcolor MidnightBlue;
draw (        0.6289565053u,        0.9882172714u)--(        0.6550000000u,        0.9750000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6289565053u,        0.9882172714u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6289565053u,        0.9882172714u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6550000000u,        0.9750000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6550000000u,        0.9750000000u) withcolor Yellow;
draw (        0.6289565053u,        0.9882172714u)--(        0.6332500533u,        1.0170716833u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6611669001u,        1.0255369832u)--(        0.6332500533u,        1.0170716833u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6289565053u,        0.9882172714u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6289565053u,        0.9882172714u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6611669001u,        1.0255369832u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6611669001u,        1.0255369832u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6332500533u,        1.0170716833u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6332500533u,        1.0170716833u) withcolor MidnightBlue;
draw (        0.6289565053u,        0.9882172714u)--(        0.6007618847u,        0.9807289467u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5777054129u,        0.9986010128u)--(        0.6007618847u,        0.9807289467u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6289565053u,        0.9882172714u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6289565053u,        0.9882172714u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5777054129u,        0.9986010128u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5777054129u,        0.9986010128u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6007618847u,        0.9807289467u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6007618847u,        0.9807289467u) withcolor MidnightBlue;
draw (        0.6289565053u,        0.9882172714u)--(        0.6059000335u,        1.0060893375u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5777054129u,        0.9986010128u)--(        0.6059000335u,        1.0060893375u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6289565053u,        0.9882172714u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6289565053u,        0.9882172714u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5777054129u,        0.9986010128u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5777054129u,        0.9986010128u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6059000335u,        1.0060893375u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6059000335u,        1.0060893375u) withcolor MidnightBlue;
draw (        0.6289565053u,        0.9882172714u)--(        0.6059078049u,        1.0060993588u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6101948595u,        1.0349547361u)--(        0.6059078049u,        1.0060993588u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6289565053u,        0.9882172714u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6289565053u,        0.9882172714u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6101948595u,        1.0349547361u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6101948595u,        1.0349547361u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6059078049u,        1.0060993588u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6059078049u,        1.0060993588u) withcolor MidnightBlue;
draw (        0.6289565053u,        0.9882172714u)--(        0.6332435599u,        1.0170726487u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6101948595u,        1.0349547361u)--(        0.6332435599u,        1.0170726487u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6289565053u,        0.9882172714u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6289565053u,        0.9882172714u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6101948595u,        1.0349547361u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6101948595u,        1.0349547361u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6332435599u,        1.0170726487u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6332435599u,        1.0170726487u) withcolor MidnightBlue;
draw (        0.6611669001u,        1.0255369832u)--(        0.6900000000u,        1.0300000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6611669001u,        1.0255369832u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6611669001u,        1.0255369832u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6900000000u,        1.0300000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6900000000u,        1.0300000000u) withcolor Yellow;
draw (        0.6611669001u,        1.0255369832u)--(        0.6332480184u,        1.0170783968u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6101948595u,        1.0349547361u)--(        0.6332480184u,        1.0170783968u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6611669001u,        1.0255369832u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6611669001u,        1.0255369832u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6101948595u,        1.0349547361u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6101948595u,        1.0349547361u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6332480184u,        1.0170783968u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6332480184u,        1.0170783968u) withcolor MidnightBlue;
draw (        0.6611669001u,        1.0255369832u)--(        0.6625827066u,        1.0546747101u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6339862043u,        1.0604411385u)--(        0.6625827066u,        1.0546747101u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6611669001u,        1.0255369832u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6611669001u,        1.0255369832u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6339862043u,        1.0604411385u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6339862043u,        1.0604411385u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6625827066u,        1.0546747101u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6625827066u,        1.0546747101u) withcolor Green;
draw (        0.4280920558u,        0.7709921534u)--(        0.4250000000u,        0.8000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4280920558u,        0.7709921534u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4280920558u,        0.7709921534u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4250000000u,        0.8000000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4250000000u,        0.8000000000u) withcolor Yellow;
draw (        0.4280920558u,        0.7709921534u)--(        0.4004581750u,        0.7803398950u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3973659469u,        0.8093476485u)--(        0.4004581750u,        0.7803398950u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4280920558u,        0.7709921534u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4280920558u,        0.7709921534u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3973659469u,        0.8093476485u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3973659469u,        0.8093476485u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4004581750u,        0.7803398950u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4004581750u,        0.7803398950u) withcolor Green;
draw (        0.4280920558u,        0.7709921534u)--(        0.4249998276u,        0.7999999069u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3973659469u,        0.8093476485u)--(        0.4249998276u,        0.7999999069u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4280920558u,        0.7709921534u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4280920558u,        0.7709921534u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3973659469u,        0.8093476485u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3973659469u,        0.8093476485u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4249998276u,        0.7999999069u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4249998276u,        0.7999999069u) withcolor Green;
draw (        0.4696277763u,        0.7941263670u)--(        0.4700000000u,        0.7650000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4696277763u,        0.7941263670u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4696277763u,        0.7941263670u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4700000000u,        0.7650000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4700000000u,        0.7650000000u) withcolor Yellow;
draw (        0.4871977584u,        0.8801053301u)--(        0.5134475415u,        0.8928317014u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5115937696u,        0.9219448456u)--(        0.5134475415u,        0.8928317014u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4871977584u,        0.8801053301u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4871977584u,        0.8801053301u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5115937696u,        0.9219448456u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5115937696u,        0.9219448456u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5134475415u,        0.8928317014u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5134475415u,        0.8928317014u) withcolor MidnightBlue;
draw (        0.4871977584u,        0.8801053301u)--(        0.4660577909u,        0.8600027578u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4412178471u,        0.8752994513u)--(        0.4660577909u,        0.8600027578u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4871977584u,        0.8801053301u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4871977584u,        0.8801053301u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4412178471u,        0.8752994513u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4412178471u,        0.8752994513u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4660577909u,        0.8600027578u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4660577909u,        0.8600027578u) withcolor Green;
draw (        0.4871977584u,        0.8801053301u)--(        0.4623578146u,        0.8954020236u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4412178471u,        0.8752994513u)--(        0.4623578146u,        0.8954020236u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4871977584u,        0.8801053301u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4871977584u,        0.8801053301u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4412178471u,        0.8752994513u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4412178471u,        0.8752994513u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4623578146u,        0.8954020236u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4623578146u,        0.8954020236u) withcolor MidnightBlue;
draw (        0.4871977584u,        0.8801053301u)--(        0.4623585342u,        0.8954031921u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4651478441u,        0.9244416388u)--(        0.4623585342u,        0.8954031921u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4871977584u,        0.8801053301u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4871977584u,        0.8801053301u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4651478441u,        0.9244416388u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4651478441u,        0.9244416388u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4623585342u,        0.8954031921u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4623585342u,        0.8954031921u) withcolor MidnightBlue;
draw (        0.5115937696u,        0.9219448456u)--(        0.5186310798u,        0.9502554073u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5438378022u,        0.9649398475u)--(        0.5186310798u,        0.9502554073u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5115937696u,        0.9219448456u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5115937696u,        0.9219448456u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5438378022u,        0.9649398475u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5438378022u,        0.9649398475u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5186310798u,        0.9502554073u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5186310798u,        0.9502554073u) withcolor MidnightBlue;
draw (        0.5115937696u,        0.9219448456u)--(        0.5368004920u,        0.9366292858u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5438378022u,        0.9649398475u)--(        0.5368004920u,        0.9366292858u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5115937696u,        0.9219448456u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5115937696u,        0.9219448456u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5438378022u,        0.9649398475u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5438378022u,        0.9649398475u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5368004920u,        0.9366292858u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5368004920u,        0.9366292858u) withcolor MidnightBlue;
draw (        0.5115937696u,        0.9219448456u)--(        0.4893161560u,        0.9407788460u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4651478441u,        0.9244416388u)--(        0.4893161560u,        0.9407788460u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5115937696u,        0.9219448456u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5115937696u,        0.9219448456u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4651478441u,        0.9244416388u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4651478441u,        0.9244416388u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4893161560u,        0.9407788460u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4893161560u,        0.9407788460u) withcolor MidnightBlue;
draw (        0.5115937696u,        0.9219448456u)--(        0.4893219539u,        0.9407857018u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4963590693u,        0.9690963120u)--(        0.4893219539u,        0.9407857018u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5115937696u,        0.9219448456u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5115937696u,        0.9219448456u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4893219539u,        0.9407857018u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4893219539u,        0.9407857018u) withcolor MidnightBlue;
draw (        0.5115937696u,        0.9219448456u)--(        0.5186308849u,        0.9502554557u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4963590693u,        0.9690963120u)--(        0.5186308849u,        0.9502554557u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5115937696u,        0.9219448456u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5115937696u,        0.9219448456u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5186308849u,        0.9502554557u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5186308849u,        0.9502554557u) withcolor MidnightBlue;
draw (        0.5438378022u,        0.9649398475u)--(        0.5489548540u,        0.9936596563u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5777054129u,        0.9986010128u)--(        0.5489548540u,        0.9936596563u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5438378022u,        0.9649398475u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5438378022u,        0.9649398475u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5777054129u,        0.9986010128u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5777054129u,        0.9986010128u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5489548540u,        0.9936596563u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5489548540u,        0.9936596563u) withcolor MidnightBlue;
draw (        0.5438378022u,        0.9649398475u)--(        0.5186309836u,        0.9502555723u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4963590693u,        0.9690963120u)--(        0.5186309836u,        0.9502555723u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5438378022u,        0.9649398475u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5438378022u,        0.9649398475u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5186309836u,        0.9502555723u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5186309836u,        0.9502555723u) withcolor MidnightBlue;
draw (        0.5438378022u,        0.9649398475u)--(        0.5215658878u,        0.9837805871u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4963590693u,        0.9690963120u)--(        0.5215658878u,        0.9837805871u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5438378022u,        0.9649398475u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5438378022u,        0.9649398475u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5215658878u,        0.9837805871u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5215658878u,        0.9837805871u) withcolor MidnightBlue;
draw (        0.5438378022u,        0.9649398475u)--(        0.5215699360u,        0.9837853715u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5266752110u,        1.0125072761u)--(        0.5215699360u,        0.9837853715u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5438378022u,        0.9649398475u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5438378022u,        0.9649398475u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5266752110u,        1.0125072761u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5266752110u,        1.0125072761u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5215699360u,        0.9837853715u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5215699360u,        0.9837853715u) withcolor MidnightBlue;
draw (        0.5438378022u,        0.9649398475u)--(        0.5489430772u,        0.9936617521u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5266752110u,        1.0125072761u)--(        0.5489430772u,        0.9936617521u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5438378022u,        0.9649398475u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5438378022u,        0.9649398475u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5266752110u,        1.0125072761u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5266752110u,        1.0125072761u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5489430772u,        0.9936617521u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5489430772u,        0.9936617521u) withcolor MidnightBlue;
draw (        0.5777054129u,        0.9986010128u)--(        0.6058969400u,        1.0061009751u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6101948595u,        1.0349547361u)--(        0.6058969400u,        1.0061009751u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5777054129u,        0.9986010128u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5777054129u,        0.9986010128u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6101948595u,        1.0349547361u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6101948595u,        1.0349547361u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6058969400u,        1.0061009751u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6058969400u,        1.0061009751u) withcolor MidnightBlue;
draw (        0.5777054129u,        0.9986010128u)--(        0.5489525685u,        0.9936729728u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5266752110u,        1.0125072761u)--(        0.5489525685u,        0.9936729728u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5777054129u,        0.9986010128u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5777054129u,        0.9986010128u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5266752110u,        1.0125072761u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5266752110u,        1.0125072761u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5489525685u,        0.9936729728u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5489525685u,        0.9936729728u) withcolor MidnightBlue;
draw (        0.5777054129u,        0.9986010128u)--(        0.5600012095u,        1.0217866294u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5799642533u,        1.0430584069u)--(        0.5600012095u,        1.0217866294u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5777054129u,        0.9986010128u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5777054129u,        0.9986010128u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5799642533u,        1.0430584069u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5799642533u,        1.0430584069u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5600012095u,        1.0217866294u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5600012095u,        1.0217866294u) withcolor Green;
draw (        0.6101948595u,        1.0349547361u)--(        0.6049921351u,        1.0636591486u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6339862043u,        1.0604411385u)--(        0.6049921351u,        1.0636591486u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6101948595u,        1.0349547361u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6101948595u,        1.0349547361u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6339862043u,        1.0604411385u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6339862043u,        1.0604411385u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6049921351u,        1.0636591486u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6049921351u,        1.0636591486u) withcolor MidnightBlue;
draw (        0.6101948595u,        1.0349547361u)--(        0.6014540614u,        1.0627865540u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5799642533u,        1.0430584069u)--(        0.6014540614u,        1.0627865540u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6101948595u,        1.0349547361u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6101948595u,        1.0349547361u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5799642533u,        1.0430584069u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5799642533u,        1.0430584069u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6014540614u,        1.0627865540u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6014540614u,        1.0627865540u) withcolor MidnightBlue;
draw (        0.6101948595u,        1.0349547361u)--(        0.6014642725u,        1.0627897588u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5962564363u,        1.0914932443u)--(        0.6014642725u,        1.0627897588u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6101948595u,        1.0349547361u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6101948595u,        1.0349547361u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5962564363u,        1.0914932443u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5962564363u,        1.0914932443u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6014642725u,        1.0627897588u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6014642725u,        1.0627897588u) withcolor MidnightBlue;
draw (        0.6101948595u,        1.0349547361u)--(        0.6049870233u,        1.0636582216u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5962564363u,        1.0914932443u)--(        0.6049870233u,        1.0636582216u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6101948595u,        1.0349547361u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6101948595u,        1.0349547361u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5962564363u,        1.0914932443u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5962564363u,        1.0914932443u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6049870233u,        1.0636582216u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6049870233u,        1.0636582216u) withcolor MidnightBlue;
draw (        0.6339862043u,        1.0604411385u)--(        0.6049922130u,        1.0636598499u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5962564363u,        1.0914932443u)--(        0.6049922130u,        1.0636598499u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6339862043u,        1.0604411385u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6339862043u,        1.0604411385u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5962564363u,        1.0914932443u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5962564363u,        1.0914932443u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6049922130u,        1.0636598499u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6049922130u,        1.0636598499u) withcolor MidnightBlue;
draw (        0.6339862043u,        1.0604411385u)--(        0.6252504277u,        1.0882745329u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5962564363u,        1.0914932443u)--(        0.6252504277u,        1.0882745329u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6339862043u,        1.0604411385u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6339862043u,        1.0604411385u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5962564363u,        1.0914932443u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5962564363u,        1.0914932443u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6252504277u,        1.0882745329u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6252504277u,        1.0882745329u) withcolor Green;
draw (        0.3973659469u,        0.8093476485u)--(        0.4250000000u,        0.8000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3973659469u,        0.8093476485u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3973659469u,        0.8093476485u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4250000000u,        0.8000000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4250000000u,        0.8000000000u) withcolor Yellow;
draw (        0.3973659469u,        0.8093476485u)--(        0.3889464468u,        0.8372783423u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4154396362u,        0.8494899169u)--(        0.3889464468u,        0.8372783423u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3973659469u,        0.8093476485u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3973659469u,        0.8093476485u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4154396362u,        0.8494899169u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4154396362u,        0.8494899169u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3889464468u,        0.8372783423u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3889464468u,        0.8372783423u) withcolor MidnightBlue;
draw (        0.3973659469u,        0.8093476485u)--(        0.3889446065u,        0.8372777875u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3604781949u,        0.8436554871u)--(        0.3889446065u,        0.8372777875u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3973659469u,        0.8093476485u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3973659469u,        0.8093476485u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3604781949u,        0.8436554871u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3604781949u,        0.8436554871u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3889446065u,        0.8372777875u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3889446065u,        0.8372777875u) withcolor MidnightBlue;
draw (        0.4154396362u,        0.8494899169u)--(        0.3889452920u,        0.8372808481u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3604781949u,        0.8436554871u)--(        0.3889452920u,        0.8372808481u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4154396362u,        0.8494899169u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4154396362u,        0.8494899169u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3604781949u,        0.8436554871u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3604781949u,        0.8436554871u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3889452920u,        0.8372808481u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3889452920u,        0.8372808481u) withcolor MidnightBlue;
draw (        0.4154396362u,        0.8494899169u)--(        0.3869725391u,        0.8558645559u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3604781949u,        0.8436554871u)--(        0.3869725391u,        0.8558645559u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4154396362u,        0.8494899169u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4154396362u,        0.8494899169u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3604781949u,        0.8436554871u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3604781949u,        0.8436554871u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3869725391u,        0.8558645559u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3869725391u,        0.8558645559u) withcolor MidnightBlue;
draw (        0.4154396362u,        0.8494899169u)--(        0.3869753503u,        0.8558770966u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3908073442u,        0.8847964234u)--(        0.3869753503u,        0.8558770966u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4154396362u,        0.8494899169u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4154396362u,        0.8494899169u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3908073442u,        0.8847964234u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3908073442u,        0.8847964234u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3869753503u,        0.8558770966u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3869753503u,        0.8558770966u) withcolor MidnightBlue;
draw (        0.4412178471u,        0.8752994513u)--(        0.4440091720u,        0.9043377044u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4651478441u,        0.9244416388u)--(        0.4440091720u,        0.9043377044u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4412178471u,        0.8752994513u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4412178471u,        0.8752994513u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4651478441u,        0.9244416388u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4651478441u,        0.9244416388u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4440091720u,        0.9043377044u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4440091720u,        0.9043377044u) withcolor MidnightBlue;
draw (        0.4412178471u,        0.8752994513u)--(        0.4623565193u,        0.8954033857u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4651478441u,        0.9244416388u)--(        0.4623565193u,        0.8954033857u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4412178471u,        0.8752994513u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4412178471u,        0.8752994513u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4651478441u,        0.9244416388u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4651478441u,        0.9244416388u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4623565193u,        0.8954033857u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4623565193u,        0.8954033857u) withcolor MidnightBlue;
draw (        0.4412178471u,        0.8752994513u)--(        0.4185856447u,        0.8937058374u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3908073442u,        0.8847964234u)--(        0.4185856447u,        0.8937058374u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4412178471u,        0.8752994513u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4412178471u,        0.8752994513u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3908073442u,        0.8847964234u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3908073442u,        0.8847964234u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4185856447u,        0.8937058374u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4185856447u,        0.8937058374u) withcolor MidnightBlue;
draw (        0.4412178471u,        0.8752994513u)--(        0.4185927125u,        0.8937145245u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4213655894u,        0.9227545450u)--(        0.4185927125u,        0.8937145245u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4412178471u,        0.8752994513u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4412178471u,        0.8752994513u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4213655894u,        0.9227545450u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4213655894u,        0.9227545450u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4185927125u,        0.8937145245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4185927125u,        0.8937145245u) withcolor MidnightBlue;
draw (        0.4412178471u,        0.8752994513u)--(        0.4439907240u,        0.9043394718u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4213655894u,        0.9227545450u)--(        0.4439907240u,        0.9043394718u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4412178471u,        0.8752994513u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4412178471u,        0.8752994513u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4213655894u,        0.9227545450u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4213655894u,        0.9227545450u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4439907240u,        0.9043394718u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4439907240u,        0.9043394718u) withcolor MidnightBlue;
draw (        0.4651478441u,        0.9244416388u)--(        0.4721975047u,        0.9527491277u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4963590693u,        0.9690963120u)--(        0.4721975047u,        0.9527491277u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4651478441u,        0.9244416388u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4651478441u,        0.9244416388u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4721975047u,        0.9527491277u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4721975047u,        0.9527491277u) withcolor MidnightBlue;
draw (        0.4651478441u,        0.9244416388u)--(        0.4893094087u,        0.9407888231u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4963590693u,        0.9690963120u)--(        0.4893094087u,        0.9407888231u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4651478441u,        0.9244416388u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4651478441u,        0.9244416388u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4893094087u,        0.9407888231u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4893094087u,        0.9407888231u) withcolor MidnightBlue;
draw (        0.4651478441u,        0.9244416388u)--(        0.4439984571u,        0.9043489769u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4213655894u,        0.9227545450u)--(        0.4439984571u,        0.9043489769u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4651478441u,        0.9244416388u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4651478441u,        0.9244416388u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4213655894u,        0.9227545450u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4213655894u,        0.9227545450u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4439984571u,        0.9043489769u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4439984571u,        0.9043489769u) withcolor MidnightBlue;
draw (        0.4651478441u,        0.9244416388u)--(        0.4721915216u,        0.9527506170u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4492648580u,        0.9707889022u)--(        0.4721915216u,        0.9527506170u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4651478441u,        0.9244416388u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4651478441u,        0.9244416388u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4492648580u,        0.9707889022u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4492648580u,        0.9707889022u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4721915216u,        0.9527506170u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4721915216u,        0.9527506170u) withcolor MidnightBlue;
draw (        0.4963590693u,        0.9690963120u)--(        0.5014719710u,        0.9978168599u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5266752110u,        1.0125072761u)--(        0.5014719710u,        0.9978168599u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5266752110u,        1.0125072761u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5266752110u,        1.0125072761u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5014719710u,        0.9978168599u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5014719710u,        0.9978168599u) withcolor MidnightBlue;
draw (        0.4963590693u,        0.9690963120u)--(        0.5215623093u,        0.9837867282u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5266752110u,        1.0125072761u)--(        0.5215623093u,        0.9837867282u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5266752110u,        1.0125072761u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5266752110u,        1.0125072761u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5215623093u,        0.9837867282u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5215623093u,        0.9837867282u) withcolor MidnightBlue;
draw (        0.4963590693u,        0.9690963120u)--(        0.4721941968u,        0.9527540178u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4492648580u,        0.9707889022u)--(        0.4721941968u,        0.9527540178u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4492648580u,        0.9707889022u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4492648580u,        0.9707889022u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4721941968u,        0.9527540178u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4721941968u,        0.9527540178u) withcolor MidnightBlue;
draw (        0.4963590693u,        0.9690963120u)--(        0.4734297304u,        0.9871311964u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4492648580u,        0.9707889022u)--(        0.4734297304u,        0.9871311964u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4492648580u,        0.9707889022u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4492648580u,        0.9707889022u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4734297304u,        0.9871311964u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4734297304u,        0.9871311964u) withcolor MidnightBlue;
draw (        0.4963590693u,        0.9690963120u)--(        0.4734300383u,        0.9871315879u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4785308756u,        1.0158542810u)--(        0.4734300383u,        0.9871315879u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4785308756u,        1.0158542810u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4785308756u,        1.0158542810u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4734300383u,        0.9871315879u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4734300383u,        0.9871315879u) withcolor MidnightBlue;
draw (        0.4963590693u,        0.9690963120u)--(        0.5014599065u,        0.9978190050u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4785308756u,        1.0158542810u)--(        0.5014599065u,        0.9978190050u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4963590693u,        0.9690963120u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4785308756u,        1.0158542810u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4785308756u,        1.0158542810u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5014599065u,        0.9978190050u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5014599065u,        0.9978190050u) withcolor MidnightBlue;
draw (        0.5266752110u,        1.0125072761u)--(        0.5014661184u,        0.9978269053u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4785308756u,        1.0158542810u)--(        0.5014661184u,        0.9978269053u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5266752110u,        1.0125072761u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5266752110u,        1.0125072761u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4785308756u,        1.0158542810u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4785308756u,        1.0158542810u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5014661184u,        0.9978269053u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5014661184u,        0.9978269053u) withcolor MidnightBlue;
draw (        0.5266752110u,        1.0125072761u)--(        0.5037399682u,        1.0305346518u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4785308756u,        1.0158542810u)--(        0.5037399682u,        1.0305346518u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5266752110u,        1.0125072761u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5266752110u,        1.0125072761u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4785308756u,        1.0158542810u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4785308756u,        1.0158542810u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5037399682u,        1.0305346518u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5037399682u,        1.0305346518u) withcolor MidnightBlue;
draw (        0.5266752110u,        1.0125072761u)--(        0.5037422896u,        1.0305376048u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5081104058u,        1.0593808224u)--(        0.5037422896u,        1.0305376048u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5266752110u,        1.0125072761u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5266752110u,        1.0125072761u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5081104058u,        1.0593808224u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5081104058u,        1.0593808224u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5037422896u,        1.0305376048u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5037422896u,        1.0305376048u) withcolor MidnightBlue;
draw (        0.5799642533u,        1.0430584069u)--(        0.5747676941u,        1.0717639362u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5962564363u,        1.0914932443u)--(        0.5747676941u,        1.0717639362u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5799642533u,        1.0430584069u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5799642533u,        1.0430584069u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5962564363u,        1.0914932443u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5962564363u,        1.0914932443u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5747676941u,        1.0717639362u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5747676941u,        1.0717639362u) withcolor Green;
draw (        0.5799642533u,        1.0430584069u)--(        0.6014529955u,        1.0627877150u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5962564363u,        1.0914932443u)--(        0.6014529955u,        1.0627877150u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5799642533u,        1.0430584069u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5799642533u,        1.0430584069u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5962564363u,        1.0914932443u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5962564363u,        1.0914932443u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6014529955u,        1.0627877150u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6014529955u,        1.0627877150u) withcolor MidnightBlue;
draw (        0.5962564363u,        1.0914932443u)--(        0.5900000000u,        1.1200000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5962564363u,        1.0914932443u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5962564363u,        1.0914932443u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5900000000u,        1.1200000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5900000000u,        1.1200000000u) withcolor Yellow;
draw (        0.3604781949u,        0.8436554871u)--(        0.3550000000u,        0.8150000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3604781949u,        0.8436554871u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3604781949u,        0.8436554871u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3550000000u,        0.8150000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3550000000u,        0.8150000000u) withcolor Yellow;
draw (        0.3604781949u,        0.8436554871u)--(        0.3869661962u,        0.8558783111u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3908073442u,        0.8847964234u)--(        0.3869661962u,        0.8558783111u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3604781949u,        0.8436554871u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3604781949u,        0.8436554871u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3908073442u,        0.8847964234u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3908073442u,        0.8847964234u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3869661962u,        0.8558783111u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3869661962u,        0.8558783111u) withcolor MidnightBlue;
draw (        0.3908073442u,        0.8847964234u)--(        0.4185825446u,        0.8937154972u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4213655894u,        0.9227545450u)--(        0.4185825446u,        0.8937154972u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3908073442u,        0.8847964234u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3908073442u,        0.8847964234u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4213655894u,        0.9227545450u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4213655894u,        0.9227545450u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4185825446u,        0.8937154972u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4185825446u,        0.8937154972u) withcolor MidnightBlue;
draw (        0.3908073442u,        0.8847964234u)--(        0.3670947726u,        0.9017883505u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3419840060u,        0.8869404202u)--(        0.3670947726u,        0.9017883505u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3908073442u,        0.8847964234u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3908073442u,        0.8847964234u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3419840060u,        0.8869404202u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3419840060u,        0.8869404202u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3670947726u,        0.9017883505u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3670947726u,        0.9017883505u) withcolor MidnightBlue;
draw (        0.3908073442u,        0.8847964234u)--(        0.3671008266u,        0.9017967959u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3757244663u,        0.9296651377u)--(        0.3671008266u,        0.9017967959u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3908073442u,        0.8847964234u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3908073442u,        0.8847964234u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3757244663u,        0.9296651377u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3757244663u,        0.9296651377u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3671008266u,        0.9017967959u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3671008266u,        0.9017967959u) withcolor MidnightBlue;
draw (        0.4213655894u,        0.9227545450u)--(        0.4276005281u,        0.9512525656u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4492648580u,        0.9707889022u)--(        0.4276005281u,        0.9512525656u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4213655894u,        0.9227545450u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4213655894u,        0.9227545450u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4492648580u,        0.9707889022u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4492648580u,        0.9707889022u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4276005281u,        0.9512525656u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4276005281u,        0.9512525656u) withcolor MidnightBlue;
draw (        0.4213655894u,        0.9227545450u)--(        0.4012158906u,        0.9438495984u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3757244663u,        0.9296651377u)--(        0.4012158906u,        0.9438495984u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4213655894u,        0.9227545450u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4213655894u,        0.9227545450u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3757244663u,        0.9296651377u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3757244663u,        0.9296651377u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4012158906u,        0.9438495984u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4012158906u,        0.9438495984u) withcolor MidnightBlue;
draw (        0.4213655894u,        0.9227545450u)--(        0.4012264980u,        0.9438597254u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4074541434u,        0.9723593407u)--(        0.4012264980u,        0.9438597254u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4213655894u,        0.9227545450u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4213655894u,        0.9227545450u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4074541434u,        0.9723593407u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4074541434u,        0.9723593407u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4012264980u,        0.9438597254u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4012264980u,        0.9438597254u) withcolor MidnightBlue;
draw (        0.4213655894u,        0.9227545450u)--(        0.4275932347u,        0.9512541603u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4074541434u,        0.9723593407u)--(        0.4275932347u,        0.9512541603u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4213655894u,        0.9227545450u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4213655894u,        0.9227545450u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4074541434u,        0.9723593407u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4074541434u,        0.9723593407u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4275932347u,        0.9512541603u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4275932347u,        0.9512541603u) withcolor MidnightBlue;
draw (        0.4492648580u,        0.9707889022u)--(        0.4543663461u,        0.9995114797u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4785308756u,        1.0158542810u)--(        0.4543663461u,        0.9995114797u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4492648580u,        0.9707889022u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4492648580u,        0.9707889022u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4785308756u,        1.0158542810u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4785308756u,        1.0158542810u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4543663461u,        0.9995114797u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4543663461u,        0.9995114797u) withcolor MidnightBlue;
draw (        0.4492648580u,        0.9707889022u)--(        0.4734293874u,        0.9871317035u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4785308756u,        1.0158542810u)--(        0.4734293874u,        0.9871317035u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4492648580u,        0.9707889022u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4492648580u,        0.9707889022u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4785308756u,        1.0158542810u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4785308756u,        1.0158542810u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4734293874u,        0.9871317035u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4734293874u,        0.9871317035u) withcolor MidnightBlue;
draw (        0.4492648580u,        0.9707889022u)--(        0.4275963822u,        0.9512571641u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4074541434u,        0.9723593407u)--(        0.4275963822u,        0.9512571641u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4492648580u,        0.9707889022u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4492648580u,        0.9707889022u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4074541434u,        0.9723593407u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4074541434u,        0.9723593407u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4275963822u,        0.9512571641u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4275963822u,        0.9512571641u) withcolor MidnightBlue;
draw (        0.4492648580u,        0.9707889022u)--(        0.4291226192u,        0.9918910788u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4074541434u,        0.9723593407u)--(        0.4291226192u,        0.9918910788u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4492648580u,        0.9707889022u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4492648580u,        0.9707889022u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4074541434u,        0.9723593407u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4074541434u,        0.9723593407u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4291226192u,        0.9918910788u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4291226192u,        0.9918910788u) withcolor MidnightBlue;
draw (        0.4492648580u,        0.9707889022u)--(        0.4291230250u,        0.9918914661u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4342184410u,        1.0206151214u)--(        0.4291230250u,        0.9918914661u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4492648580u,        0.9707889022u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4492648580u,        0.9707889022u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4291230250u,        0.9918914661u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4291230250u,        0.9918914661u) withcolor MidnightBlue;
draw (        0.4492648580u,        0.9707889022u)--(        0.4543602740u,        0.9995125575u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4342184410u,        1.0206151214u)--(        0.4543602740u,        0.9995125575u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4492648580u,        0.9707889022u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4492648580u,        0.9707889022u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4543602740u,        0.9995125575u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4543602740u,        0.9995125575u) withcolor MidnightBlue;
draw (        0.4785308756u,        1.0158542810u)--(        0.5037378570u,        1.0305382764u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5081104058u,        1.0593808224u)--(        0.5037378570u,        1.0305382764u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4785308756u,        1.0158542810u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4785308756u,        1.0158542810u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5081104058u,        1.0593808224u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5081104058u,        1.0593808224u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5037378570u,        1.0305382764u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5037378570u,        1.0305382764u) withcolor MidnightBlue;
draw (        0.4785308756u,        1.0158542810u)--(        0.4543635212u,        0.9995156574u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4342184410u,        1.0206151214u)--(        0.4543635212u,        0.9995156574u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4785308756u,        1.0158542810u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4785308756u,        1.0158542810u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4543635212u,        0.9995156574u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4543635212u,        0.9995156574u) withcolor MidnightBlue;
draw (        0.4785308756u,        1.0158542810u)--(        0.4583857953u,        1.0369537450u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4342184410u,        1.0206151214u)--(        0.4583857953u,        1.0369537450u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4785308756u,        1.0158542810u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4785308756u,        1.0158542810u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4583857953u,        1.0369537450u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4583857953u,        1.0369537450u) withcolor MidnightBlue;
draw (        0.4785308756u,        1.0158542810u)--(        0.4583933722u,        1.0369609765u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4670658524u,        1.0648141581u)--(        0.4583933722u,        1.0369609765u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4785308756u,        1.0158542810u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4785308756u,        1.0158542810u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4670658524u,        1.0648141581u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4670658524u,        1.0648141581u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4583933722u,        1.0369609765u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4583933722u,        1.0369609765u) withcolor MidnightBlue;
draw (        0.5081104058u,        1.0593808224u)--(        0.4902854675u,        1.0824737479u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4670658524u,        1.0648141581u)--(        0.4902854675u,        1.0824737479u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5081104058u,        1.0593808224u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5081104058u,        1.0593808224u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4670658524u,        1.0648141581u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4670658524u,        1.0648141581u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4902854675u,        1.0824737479u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4902854675u,        1.0824737479u) withcolor MidnightBlue;
draw (        0.5081104058u,        1.0593808224u)--(        0.4902980650u,        1.0824834662u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4961628118u,        1.1110599668u)--(        0.4902980650u,        1.0824834662u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5081104058u,        1.0593808224u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5081104058u,        1.0593808224u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4961628118u,        1.1110599668u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4961628118u,        1.1110599668u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4902980650u,        1.0824834662u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4902980650u,        1.0824834662u) withcolor MidnightBlue;
draw (        0.3419840060u,        0.8869404202u)--(        0.3670872969u,        0.9018009861u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3757244663u,        0.9296651377u)--(        0.3670872969u,        0.9018009861u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3419840060u,        0.8869404202u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3419840060u,        0.8869404202u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3757244663u,        0.9296651377u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3757244663u,        0.9296651377u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3670872969u,        0.9018009861u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3670872969u,        0.9018009861u) withcolor MidnightBlue;
draw (        0.3757244663u,        0.9296651377u)--(        0.4012079973u,        0.9438637744u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4074541434u,        0.9723593407u)--(        0.4012079973u,        0.9438637744u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3757244663u,        0.9296651377u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3757244663u,        0.9296651377u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4074541434u,        0.9723593407u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4074541434u,        0.9723593407u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4012079973u,        0.9438637744u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4012079973u,        0.9438637744u) withcolor MidnightBlue;
draw (        0.3757244663u,        0.9296651377u)--(        0.3503167937u,        0.9439990769u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3257749524u,        0.9282285454u)--(        0.3503167937u,        0.9439990769u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3757244663u,        0.9296651377u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3757244663u,        0.9296651377u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3257749524u,        0.9282285454u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3257749524u,        0.9282285454u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3503167937u,        0.9439990769u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3503167937u,        0.9439990769u) withcolor MidnightBlue;
draw (        0.3757244663u,        0.9296651377u)--(        0.3503224299u,        0.9440090628u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3586846939u,        0.9719569461u)--(        0.3503224299u,        0.9440090628u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3757244663u,        0.9296651377u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3757244663u,        0.9296651377u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3586846939u,        0.9719569461u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3586846939u,        0.9719569461u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3503224299u,        0.9440090628u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3503224299u,        0.9440090628u) withcolor MidnightBlue;
draw (        0.4074541434u,        0.9723593407u)--(        0.4125504580u,        1.0010828365u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4342184410u,        1.0206151214u)--(        0.4125504580u,        1.0010828365u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4074541434u,        0.9723593407u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4074541434u,        0.9723593407u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4125504580u,        1.0010828365u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4125504580u,        1.0010828365u) withcolor MidnightBlue;
draw (        0.4074541434u,        0.9723593407u)--(        0.4291221264u,        0.9918916255u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4342184410u,        1.0206151214u)--(        0.4291221264u,        0.9918916255u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4074541434u,        0.9723593407u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4074541434u,        0.9723593407u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4291221264u,        0.9918916255u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4291221264u,        0.9918916255u) withcolor MidnightBlue;
draw (        0.4074541434u,        0.9723593407u)--(        0.3829373160u,        0.9881687307u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3586846939u,        0.9719569461u)--(        0.3829373160u,        0.9881687307u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4074541434u,        0.9723593407u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4074541434u,        0.9723593407u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3586846939u,        0.9719569461u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3586846939u,        0.9719569461u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3829373160u,        0.9881687307u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3829373160u,        0.9881687307u) withcolor MidnightBlue;
draw (        0.4074541434u,        0.9723593407u)--(        0.3829391958u,        0.9881716454u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3880162345u,        1.0168985547u)--(        0.3829391958u,        0.9881716454u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4074541434u,        0.9723593407u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4074541434u,        0.9723593407u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3829391958u,        0.9881716454u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3829391958u,        0.9881716454u) withcolor MidnightBlue;
draw (        0.4074541434u,        0.9723593407u)--(        0.4125311821u,        1.0010862499u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3880162345u,        1.0168985547u)--(        0.4125311821u,        1.0010862499u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4074541434u,        0.9723593407u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4074541434u,        0.9723593407u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4125311821u,        1.0010862499u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4125311821u,        1.0010862499u) withcolor MidnightBlue;
draw (        0.4342184410u,        1.0206151214u)--(        0.4429067206u,        1.0484633787u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4670658524u,        1.0648141581u)--(        0.4429067206u,        1.0484633787u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4670658524u,        1.0648141581u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4670658524u,        1.0648141581u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4429067206u,        1.0484633787u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4429067206u,        1.0484633787u) withcolor MidnightBlue;
draw (        0.4342184410u,        1.0206151214u)--(        0.4583775728u,        1.0369659008u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4670658524u,        1.0648141581u)--(        0.4583775728u,        1.0369659008u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4670658524u,        1.0648141581u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4670658524u,        1.0648141581u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4583775728u,        1.0369659008u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4583775728u,        1.0369659008u) withcolor MidnightBlue;
draw (        0.4342184410u,        1.0206151214u)--(        0.4125379400u,        1.0010967322u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3880162345u,        1.0168985547u)--(        0.4125379400u,        1.0010967322u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4125379400u,        1.0010967322u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4125379400u,        1.0010967322u) withcolor MidnightBlue;
draw (        0.4342184410u,        1.0206151214u)--(        0.4096967354u,        1.0364169439u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3880162345u,        1.0168985547u)--(        0.4096967354u,        1.0364169439u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4096967354u,        1.0364169439u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4096967354u,        1.0364169439u) withcolor MidnightBlue;
draw (        0.4342184410u,        1.0206151214u)--(        0.4096987927u,        1.0364201359u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4183806882u,        1.0642703841u)--(        0.4096987927u,        1.0364201359u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4096987927u,        1.0364201359u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4096987927u,        1.0364201359u) withcolor MidnightBlue;
draw (        0.4342184410u,        1.0206151214u)--(        0.4429003365u,        1.0484653696u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4183806882u,        1.0642703841u)--(        0.4429003365u,        1.0484653696u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4342184410u,        1.0206151214u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4429003365u,        1.0484653696u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4429003365u,        1.0484653696u) withcolor MidnightBlue;
draw (        0.4670658524u,        1.0648141581u)--(        0.4729543028u,        1.0933857837u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4961628118u,        1.1110599668u)--(        0.4729543028u,        1.0933857837u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4670658524u,        1.0648141581u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4670658524u,        1.0648141581u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4961628118u,        1.1110599668u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4961628118u,        1.1110599668u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4729543028u,        1.0933857837u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4729543028u,        1.0933857837u) withcolor MidnightBlue;
draw (        0.4670658524u,        1.0648141581u)--(        0.4902743614u,        1.0824883412u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4961628118u,        1.1110599668u)--(        0.4902743614u,        1.0824883412u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4670658524u,        1.0648141581u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4670658524u,        1.0648141581u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4961628118u,        1.1110599668u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4961628118u,        1.1110599668u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4902743614u,        1.0824883412u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4902743614u,        1.0824883412u) withcolor MidnightBlue;
draw (        0.4670658524u,        1.0648141581u)--(        0.4429027937u,        1.0484691823u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4183806882u,        1.0642703841u)--(        0.4429027937u,        1.0484691823u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4670658524u,        1.0648141581u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4670658524u,        1.0648141581u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4429027937u,        1.0484691823u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4429027937u,        1.0484691823u) withcolor MidnightBlue;
draw (        0.4670658524u,        1.0648141581u)--(        0.4425437469u,        1.0806153599u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4183806882u,        1.0642703841u)--(        0.4425437469u,        1.0806153599u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4670658524u,        1.0648141581u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4670658524u,        1.0648141581u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4425437469u,        1.0806153599u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4425437469u,        1.0806153599u) withcolor MidnightBlue;
draw (        0.4670658524u,        1.0648141581u)--(        0.4425469565u,        1.0806203398u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4484299957u,        1.1091930801u)--(        0.4425469565u,        1.0806203398u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4670658524u,        1.0648141581u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4670658524u,        1.0648141581u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4484299957u,        1.1091930801u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4484299957u,        1.1091930801u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4425469565u,        1.0806203398u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4425469565u,        1.0806203398u) withcolor MidnightBlue;
draw (        0.4670658524u,        1.0648141581u)--(        0.4729488917u,        1.0933868984u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4484299957u,        1.1091930801u)--(        0.4729488917u,        1.0933868984u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4670658524u,        1.0648141581u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4670658524u,        1.0648141581u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4484299957u,        1.1091930801u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4484299957u,        1.1091930801u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4729488917u,        1.0933868984u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4729488917u,        1.0933868984u) withcolor MidnightBlue;
draw (        0.4961628118u,        1.1110599668u)--(        0.4941565903u,        1.1401630029u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5210640139u,        1.1514325278u)--(        0.4941565903u,        1.1401630029u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4961628118u,        1.1110599668u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4961628118u,        1.1110599668u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5210640139u,        1.1514325278u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5210640139u,        1.1514325278u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4941565903u,        1.1401630029u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4941565903u,        1.1401630029u) withcolor MidnightBlue;
draw (        0.4961628118u,        1.1110599668u)--(        0.4729509833u,        1.0933901435u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4484299957u,        1.1091930801u)--(        0.4729509833u,        1.0933901435u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4961628118u,        1.1110599668u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4961628118u,        1.1110599668u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4484299957u,        1.1091930801u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4484299957u,        1.1091930801u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4729509833u,        1.0933901435u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4729509833u,        1.0933901435u) withcolor MidnightBlue;
draw (        0.4961628118u,        1.1110599668u)--(        0.4716418242u,        1.1268629034u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4484299957u,        1.1091930801u)--(        0.4716418242u,        1.1268629034u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4961628118u,        1.1110599668u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4961628118u,        1.1110599668u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4484299957u,        1.1091930801u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4484299957u,        1.1091930801u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4716418242u,        1.1268629034u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4716418242u,        1.1268629034u) withcolor MidnightBlue;
draw (        0.4961628118u,        1.1110599668u)--(        0.4716487320u,        1.1268736168u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4696310722u,        1.1559758622u)--(        0.4716487320u,        1.1268736168u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4961628118u,        1.1110599668u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4961628118u,        1.1110599668u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4696310722u,        1.1559758622u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4696310722u,        1.1559758622u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4716487320u,        1.1268736168u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4716487320u,        1.1268736168u) withcolor MidnightBlue;
draw (        0.4961628118u,        1.1110599668u)--(        0.4941451520u,        1.1401622122u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4696310722u,        1.1559758622u)--(        0.4941451520u,        1.1401622122u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4961628118u,        1.1110599668u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4961628118u,        1.1110599668u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4696310722u,        1.1559758622u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4696310722u,        1.1559758622u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4941451520u,        1.1401622122u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4941451520u,        1.1401622122u) withcolor MidnightBlue;
draw (        0.5210640139u,        1.1514325278u)--(        0.4941522882u,        1.1401732801u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4696310722u,        1.1559758622u)--(        0.4941522882u,        1.1401732801u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5210640139u,        1.1514325278u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5210640139u,        1.1514325278u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4696310722u,        1.1559758622u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4696310722u,        1.1559758622u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4941522882u,        1.1401732801u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4941522882u,        1.1401732801u) withcolor MidnightBlue;
draw (        0.5376929472u,        1.2034909989u)--(        0.5500000000u,        1.2300000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5376929472u,        1.2034909989u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5376929472u,        1.2034909989u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5500000000u,        1.2300000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5500000000u,        1.2300000000u) withcolor Yellow;
draw (        0.3257749524u,        0.9282285454u)--(        0.3503074967u,        0.9440135353u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3586846939u,        0.9719569461u)--(        0.3503074967u,        0.9440135353u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3257749524u,        0.9282285454u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3257749524u,        0.9282285454u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3586846939u,        0.9719569461u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3586846939u,        0.9719569461u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3503074967u,        0.9440135353u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3503074967u,        0.9440135353u) withcolor MidnightBlue;
draw (        0.3586846939u,        0.9719569461u)--(        0.3637660743u,        1.0006830877u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3880162345u,        1.0168985547u)--(        0.3637660743u,        1.0006830877u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3586846939u,        0.9719569461u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3586846939u,        0.9719569461u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3637660743u,        1.0006830877u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3637660743u,        1.0006830877u) withcolor MidnightBlue;
draw (        0.3586846939u,        0.9719569461u)--(        0.3829348541u,        0.9881724131u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3880162345u,        1.0168985547u)--(        0.3829348541u,        0.9881724131u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3586846939u,        0.9719569461u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3586846939u,        0.9719569461u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3829348541u,        0.9881724131u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3829348541u,        0.9881724131u) withcolor MidnightBlue;
draw (        0.3586846939u,        0.9719569461u)--(        0.3355354007u,        0.9897086181u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3141122531u,        0.9699081033u)--(        0.3355354007u,        0.9897086181u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3586846939u,        0.9719569461u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3586846939u,        0.9719569461u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3141122531u,        0.9699081033u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3141122531u,        0.9699081033u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3355354007u,        0.9897086181u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3355354007u,        0.9897086181u) withcolor MidnightBlue;
draw (        0.3586846939u,        0.9719569461u)--(        0.3355363553u,        0.9897098629u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3406062193u,        1.0184380392u)--(        0.3355363553u,        0.9897098629u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3586846939u,        0.9719569461u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3586846939u,        0.9719569461u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3406062193u,        1.0184380392u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3406062193u,        1.0184380392u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3355363553u,        0.9897098629u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3355363553u,        0.9897098629u) withcolor MidnightBlue;
draw (        0.3586846939u,        0.9719569461u)--(        0.3637545579u,        1.0006851224u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3406062193u,        1.0184380392u)--(        0.3637545579u,        1.0006851224u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3586846939u,        0.9719569461u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3586846939u,        0.9719569461u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3406062193u,        1.0184380392u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3406062193u,        1.0184380392u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3637545579u,        1.0006851224u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3637545579u,        1.0006851224u) withcolor MidnightBlue;
draw (        0.3880162345u,        1.0168985547u)--(        0.3967049875u,        1.0447466643u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4183806882u,        1.0642703841u)--(        0.3967049875u,        1.0447466643u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3967049875u,        1.0447466643u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3967049875u,        1.0447466643u) withcolor MidnightBlue;
draw (        0.3880162345u,        1.0168985547u)--(        0.4096919352u,        1.0364222745u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4183806882u,        1.0642703841u)--(        0.4096919352u,        1.0364222745u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4096919352u,        1.0364222745u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4096919352u,        1.0364222745u) withcolor MidnightBlue;
draw (        0.3880162345u,        1.0168985547u)--(        0.3637599839u,        1.0006921996u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3406062193u,        1.0184380392u)--(        0.3637599839u,        1.0006921996u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3406062193u,        1.0184380392u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3406062193u,        1.0184380392u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3637599839u,        1.0006921996u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3637599839u,        1.0006921996u) withcolor MidnightBlue;
draw (        0.3880162345u,        1.0168985547u)--(        0.3648624699u,        1.0346443943u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3406062193u,        1.0184380392u)--(        0.3648624699u,        1.0346443943u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3406062193u,        1.0184380392u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3406062193u,        1.0184380392u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3648624699u,        1.0346443943u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3648624699u,        1.0346443943u) withcolor MidnightBlue;
draw (        0.3880162345u,        1.0168985547u)--(        0.3648685316u,        1.0346523005u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3735509561u,        1.0625023839u)--(        0.3648685316u,        1.0346523005u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3735509561u,        1.0625023839u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3735509561u,        1.0625023839u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3648685316u,        1.0346523005u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3648685316u,        1.0346523005u) withcolor MidnightBlue;
draw (        0.3880162345u,        1.0168985547u)--(        0.3966986590u,        1.0447486380u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3735509561u,        1.0625023839u)--(        0.3966986590u,        1.0447486380u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3880162345u,        1.0168985547u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3735509561u,        1.0625023839u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3735509561u,        1.0625023839u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3966986590u,        1.0447486380u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3966986590u,        1.0447486380u) withcolor MidnightBlue;
draw (        0.4183806882u,        1.0642703841u)--(        0.4242713725u,        1.0928415493u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4484299957u,        1.1091930801u)--(        0.4242713725u,        1.0928415493u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4484299957u,        1.1091930801u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4484299957u,        1.1091930801u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4242713725u,        1.0928415493u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4242713725u,        1.0928415493u) withcolor MidnightBlue;
draw (        0.4183806882u,        1.0642703841u)--(        0.4425393115u,        1.0806219150u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4484299957u,        1.1091930801u)--(        0.4425393115u,        1.0806219150u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4484299957u,        1.1091930801u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4484299957u,        1.1091930801u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4425393115u,        1.0806219150u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4425393115u,        1.0806219150u) withcolor MidnightBlue;
draw (        0.4183806882u,        1.0642703841u)--(        0.3967007522u,        1.0447513675u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3735509561u,        1.0625023839u)--(        0.3967007522u,        1.0447513675u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3735509561u,        1.0625023839u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3735509561u,        1.0625023839u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3967007522u,        1.0447513675u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3967007522u,        1.0447513675u) withcolor MidnightBlue;
draw (        0.4183806882u,        1.0642703841u)--(        0.3952308922u,        1.0820214005u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3735509561u,        1.0625023839u)--(        0.3952308922u,        1.0820214005u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3735509561u,        1.0625023839u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3735509561u,        1.0625023839u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3952308922u,        1.0820214005u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3952308922u,        1.0820214005u) withcolor MidnightBlue;
draw (        0.4183806882u,        1.0642703841u)--(        0.3952344712u,        1.0820260670u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4011196204u,        1.1105983728u)--(        0.3952344712u,        1.0820260670u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4011196204u,        1.1105983728u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4011196204u,        1.1105983728u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3952344712u,        1.0820260670u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3952344712u,        1.0820260670u) withcolor MidnightBlue;
draw (        0.4183806882u,        1.0642703841u)--(        0.4242658374u,        1.0928426899u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4011196204u,        1.1105983728u)--(        0.4242658374u,        1.0928426899u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4183806882u,        1.0642703841u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4011196204u,        1.1105983728u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4011196204u,        1.1105983728u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4242658374u,        1.0928426899u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4242658374u,        1.0928426899u) withcolor MidnightBlue;
draw (        0.4484299957u,        1.1091930801u)--(        0.4716344183u,        1.1268726280u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4696310722u,        1.1559758622u)--(        0.4716344183u,        1.1268726280u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4484299957u,        1.1091930801u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4484299957u,        1.1091930801u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4696310722u,        1.1559758622u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4696310722u,        1.1559758622u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4716344183u,        1.1268726280u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4716344183u,        1.1268726280u) withcolor MidnightBlue;
draw (        0.4484299957u,        1.1091930801u)--(        0.4242683680u,        1.0928459892u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4011196204u,        1.1105983728u)--(        0.4242683680u,        1.0928459892u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4484299957u,        1.1091930801u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4484299957u,        1.1091930801u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4011196204u,        1.1105983728u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4011196204u,        1.1105983728u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4242683680u,        1.0928459892u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4242683680u,        1.0928459892u) withcolor MidnightBlue;
draw (        0.4484299957u,        1.1091930801u)--(        0.4252812482u,        1.1269454637u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4011196204u,        1.1105983728u)--(        0.4252812482u,        1.1269454637u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4484299957u,        1.1091930801u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4484299957u,        1.1091930801u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4011196204u,        1.1105983728u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4011196204u,        1.1105983728u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4252812482u,        1.1269454637u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4252812482u,        1.1269454637u) withcolor MidnightBlue;
draw (        0.4484299957u,        1.1091930801u)--(        0.4252865831u,        1.1269524181u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4256985658u,        1.1561216126u)--(        0.4252865831u,        1.1269524181u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4484299957u,        1.1091930801u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4484299957u,        1.1091930801u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4256985658u,        1.1561216126u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4256985658u,        1.1561216126u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4252865831u,        1.1269524181u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4252865831u,        1.1269524181u) withcolor MidnightBlue;
draw (        0.3141122531u,        0.9699081033u)--(        0.3191846327u,        0.9986358356u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3406062193u,        1.0184380392u)--(        0.3191846327u,        0.9986358356u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3141122531u,        0.9699081033u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3141122531u,        0.9699081033u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3406062193u,        1.0184380392u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3406062193u,        1.0184380392u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3191846327u,        0.9986358356u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3191846327u,        0.9986358356u) withcolor Green;
draw (        0.3141122531u,        0.9699081033u)--(        0.3355338397u,        0.9897103069u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3406062193u,        1.0184380392u)--(        0.3355338397u,        0.9897103069u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3141122531u,        0.9699081033u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3141122531u,        0.9699081033u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3406062193u,        1.0184380392u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3406062193u,        1.0184380392u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3355338397u,        0.9897103069u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3355338397u,        0.9897103069u) withcolor MidnightBlue;
draw (        0.3406062193u,        1.0184380392u)--(        0.3648541937u,        1.0346567745u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3735509561u,        1.0625023839u)--(        0.3648541937u,        1.0346567745u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3406062193u,        1.0184380392u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3406062193u,        1.0184380392u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3735509561u,        1.0625023839u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3735509561u,        1.0625023839u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3648541937u,        1.0346567745u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3648541937u,        1.0346567745u) withcolor MidnightBlue;
draw (        0.3735509561u,        1.0625023839u)--(        0.3794456611u,        1.0910727197u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4011196204u,        1.1105983728u)--(        0.3794456611u,        1.0910727197u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3735509561u,        1.0625023839u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3735509561u,        1.0625023839u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4011196204u,        1.1105983728u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4011196204u,        1.1105983728u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3794456611u,        1.0910727197u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3794456611u,        1.0910727197u) withcolor MidnightBlue;
draw (        0.3735509561u,        1.0625023839u)--(        0.3952249154u,        1.0820280370u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4011196204u,        1.1105983728u)--(        0.3952249154u,        1.0820280370u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3735509561u,        1.0625023839u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3735509561u,        1.0625023839u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4011196204u,        1.1105983728u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4011196204u,        1.1105983728u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3952249154u,        1.0820280370u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3952249154u,        1.0820280370u) withcolor MidnightBlue;
draw (        0.3735509561u,        1.0625023839u)--(        0.3449302704u,        1.0681475574u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3246138680u,        1.0472130062u)--(        0.3449302704u,        1.0681475574u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3735509561u,        1.0625023839u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3735509561u,        1.0625023839u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3246138680u,        1.0472130062u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3246138680u,        1.0472130062u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3449302704u,        1.0681475574u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3449302704u,        1.0681475574u) withcolor MidnightBlue;
draw (        0.3735509561u,        1.0625023839u)--(        0.3449316068u,        1.0681543291u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3508171938u,        1.0967265447u)--(        0.3449316068u,        1.0681543291u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3735509561u,        1.0625023839u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3735509561u,        1.0625023839u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3508171938u,        1.0967265447u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3508171938u,        1.0967265447u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3449316068u,        1.0681543291u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3449316068u,        1.0681543291u) withcolor MidnightBlue;
draw (        0.3735509561u,        1.0625023839u)--(        0.3794365431u,        1.0910745994u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3508171938u,        1.0967265447u)--(        0.3794365431u,        1.0910745994u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3735509561u,        1.0625023839u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3735509561u,        1.0625023839u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3508171938u,        1.0967265447u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3508171938u,        1.0967265447u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3794365431u,        1.0910745994u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3794365431u,        1.0910745994u) withcolor MidnightBlue;
draw (        0.4011196204u,        1.1105983728u)--(        0.4015417428u,        1.1397674223u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4256985658u,        1.1561216126u)--(        0.4015417428u,        1.1397674223u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4011196204u,        1.1105983728u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4011196204u,        1.1105983728u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4256985658u,        1.1561216126u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4256985658u,        1.1561216126u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4015417428u,        1.1397674223u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4015417428u,        1.1397674223u) withcolor MidnightBlue;
draw (        0.4011196204u,        1.1105983728u)--(        0.4252764434u,        1.1269525631u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4256985658u,        1.1561216126u)--(        0.4252764434u,        1.1269525631u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4011196204u,        1.1105983728u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4011196204u,        1.1105983728u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4256985658u,        1.1561216126u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4256985658u,        1.1561216126u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4252764434u,        1.1269525631u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4252764434u,        1.1269525631u) withcolor MidnightBlue;
draw (        0.4011196204u,        1.1105983728u)--(        0.3794378779u,        1.0910813627u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3508171938u,        1.0967265447u)--(        0.3794378779u,        1.0910813627u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4011196204u,        1.1105983728u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4011196204u,        1.1105983728u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3508171938u,        1.0967265447u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3508171938u,        1.0967265447u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3794378779u,        1.0910813627u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3794378779u,        1.0910813627u) withcolor MidnightBlue;
draw (        0.4011196204u,        1.1105983728u)--(        0.4015266779u,        1.1397676365u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3723548004u,        1.1396527353u)--(        0.4015266779u,        1.1397676365u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4011196204u,        1.1105983728u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4011196204u,        1.1105983728u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3723548004u,        1.1396527353u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3723548004u,        1.1396527353u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4015266779u,        1.1397676365u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4015266779u,        1.1397676365u) withcolor MidnightBlue;
draw (        0.4256985658u,        1.1561216126u)--(        0.4015265820u,        1.1397898388u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3723548004u,        1.1396527353u)--(        0.4015265820u,        1.1397898388u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4256985658u,        1.1561216126u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4256985658u,        1.1561216126u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3723548004u,        1.1396527353u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3723548004u,        1.1396527353u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4015265820u,        1.1397898388u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4015265820u,        1.1397898388u) withcolor MidnightBlue;
draw (        0.4256985658u,        1.1561216126u)--(        0.4087317575u,        1.1798521637u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3798961620u,        1.1754340114u)--(        0.4087317575u,        1.1798521637u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4256985658u,        1.1561216126u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4256985658u,        1.1561216126u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3798961620u,        1.1754340114u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3798961620u,        1.1754340114u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4087317575u,        1.1798521637u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4087317575u,        1.1798521637u) withcolor Green;
draw (        0.3246138680u,        1.0472130062u)--(        0.3449210446u,        1.0681565068u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3508171938u,        1.0967265447u)--(        0.3449210446u,        1.0681565068u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3246138680u,        1.0472130062u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3246138680u,        1.0472130062u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3508171938u,        1.0967265447u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3508171938u,        1.0967265447u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3449210446u,        1.0681565068u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3449210446u,        1.0681565068u) withcolor MidnightBlue;
draw (        0.3246138680u,        1.0472130062u)--(        0.2968773088u,        1.0562515305u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3056893420u,        1.0840608762u)--(        0.2968773088u,        1.0562515305u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3246138680u,        1.0472130062u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3246138680u,        1.0472130062u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3056893420u,        1.0840608762u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3056893420u,        1.0840608762u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2968773088u,        1.0562515305u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2968773088u,        1.0562515305u) withcolor Green;
draw (        0.3508171938u,        1.0967265447u)--(        0.3467806222u,        1.1256180266u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3723548004u,        1.1396527353u)--(        0.3467806222u,        1.1256180266u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3508171938u,        1.0967265447u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3508171938u,        1.0967265447u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3723548004u,        1.1396527353u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3723548004u,        1.1396527353u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3467806222u,        1.1256180266u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3467806222u,        1.1256180266u) withcolor MidnightBlue;
draw (        0.3508171938u,        1.0967265447u)--(        0.3235590719u,        1.1071191587u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3056893420u,        1.0840608762u)--(        0.3235590719u,        1.1071191587u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3508171938u,        1.0967265447u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3508171938u,        1.0967265447u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3056893420u,        1.0840608762u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3056893420u,        1.0840608762u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3235590719u,        1.1071191587u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3235590719u,        1.1071191587u) withcolor MidnightBlue;
draw (        0.3508171938u,        1.0967265447u)--(        0.3235602487u,        1.1071222446u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3195129567u,        1.1360122267u)--(        0.3235602487u,        1.1071222446u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3508171938u,        1.0967265447u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3508171938u,        1.0967265447u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3195129567u,        1.1360122267u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3195129567u,        1.1360122267u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3235602487u,        1.1071222446u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3235602487u,        1.1071222446u) withcolor MidnightBlue;
draw (        0.3508171938u,        1.0967265447u)--(        0.3467699019u,        1.1256165268u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3195129567u,        1.1360122267u)--(        0.3467699019u,        1.1256165268u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3508171938u,        1.0967265447u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3508171938u,        1.0967265447u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3195129567u,        1.1360122267u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3195129567u,        1.1360122267u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3467699019u,        1.1256165268u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3467699019u,        1.1256165268u) withcolor MidnightBlue;
draw (        0.3723548004u,        1.1396527353u)--(        0.3538826898u,        1.1622313262u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3798961620u,        1.1754340114u)--(        0.3538826898u,        1.1622313262u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3723548004u,        1.1396527353u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3723548004u,        1.1396527353u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3798961620u,        1.1754340114u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3798961620u,        1.1754340114u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3538826898u,        1.1622313262u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3538826898u,        1.1622313262u) withcolor MidnightBlue;
draw (        0.3723548004u,        1.1396527353u)--(        0.3467746343u,        1.1256289434u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3195129567u,        1.1360122267u)--(        0.3467746343u,        1.1256289434u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3723548004u,        1.1396527353u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3723548004u,        1.1396527353u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3195129567u,        1.1360122267u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3195129567u,        1.1360122267u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3467746343u,        1.1256289434u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3467746343u,        1.1256289434u) withcolor MidnightBlue;
draw (        0.3723548004u,        1.1396527353u)--(        0.3450931228u,        1.1500360186u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3195129567u,        1.1360122267u)--(        0.3450931228u,        1.1500360186u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3723548004u,        1.1396527353u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3723548004u,        1.1396527353u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3195129567u,        1.1360122267u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3195129567u,        1.1360122267u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3450931228u,        1.1500360186u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3450931228u,        1.1500360186u) withcolor MidnightBlue;
draw (        0.3723548004u,        1.1396527353u)--(        0.3450967025u,        1.1500454124u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3266222792u,        1.1726221110u)--(        0.3450967025u,        1.1500454124u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3723548004u,        1.1396527353u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3723548004u,        1.1396527353u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3266222792u,        1.1726221110u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3266222792u,        1.1726221110u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3450967025u,        1.1500454124u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3450967025u,        1.1500454124u) withcolor MidnightBlue;
draw (        0.3723548004u,        1.1396527353u)--(        0.3538803771u,        1.1622294339u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3266222792u,        1.1726221110u)--(        0.3538803771u,        1.1622294339u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3723548004u,        1.1396527353u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3723548004u,        1.1396527353u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3266222792u,        1.1726221110u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3266222792u,        1.1726221110u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3538803771u,        1.1622294339u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3538803771u,        1.1622294339u) withcolor MidnightBlue;
draw (        0.3798961620u,        1.1754340114u)--(        0.3538817810u,        1.1622331169u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3266222792u,        1.1726221110u)--(        0.3538817810u,        1.1622331169u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3798961620u,        1.1754340114u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3798961620u,        1.1754340114u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3266222792u,        1.1726221110u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3266222792u,        1.1726221110u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3538817810u,        1.1622331169u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3538817810u,        1.1622331169u) withcolor MidnightBlue;
draw (        0.3798961620u,        1.1754340114u)--(        0.3526366602u,        1.1858230055u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3266222792u,        1.1726221110u)--(        0.3526366602u,        1.1858230055u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3798961620u,        1.1754340114u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3798961620u,        1.1754340114u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3266222792u,        1.1726221110u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3266222792u,        1.1726221110u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3526366602u,        1.1858230055u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3526366602u,        1.1858230055u) withcolor Green;
draw (        0.3056893420u,        1.0840608762u)--(        0.3235558791u,        1.1071216328u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3195129567u,        1.1360122267u)--(        0.3235558791u,        1.1071216328u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3056893420u,        1.0840608762u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3056893420u,        1.0840608762u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3195129567u,        1.1360122267u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3195129567u,        1.1360122267u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3235558791u,        1.1071216328u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3235558791u,        1.1071216328u) withcolor MidnightBlue;
draw (        0.3056893420u,        1.0840608762u)--(        0.2781538404u,        1.0936945560u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2751812924u,        1.1227148179u)--(        0.2781538404u,        1.0936945560u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3056893420u,        1.0840608762u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3056893420u,        1.0840608762u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2751812924u,        1.1227148179u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2751812924u,        1.1227148179u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2781538404u,        1.0936945560u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2781538404u,        1.0936945560u) withcolor Green;
draw (        0.3195129567u,        1.1360122267u)--(        0.3450906755u,        1.1500404818u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3266222792u,        1.1726221110u)--(        0.3450906755u,        1.1500404818u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3195129567u,        1.1360122267u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3195129567u,        1.1360122267u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3266222792u,        1.1726221110u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3266222792u,        1.1726221110u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3450906755u,        1.1500404818u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3450906755u,        1.1500404818u) withcolor MidnightBlue;
draw (        0.3195129567u,        1.1360122267u)--(        0.2922440391u,        1.1463764812u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2751812924u,        1.1227148179u)--(        0.2922440391u,        1.1463764812u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3195129567u,        1.1360122267u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3195129567u,        1.1360122267u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2751812924u,        1.1227148179u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2751812924u,        1.1227148179u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2922440391u,        1.1463764812u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2922440391u,        1.1463764812u) withcolor MidnightBlue;
draw (        0.3195129567u,        1.1360122267u)--(        0.2922506936u,        1.1463939726u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2787893569u,        1.1722745455u)--(        0.2922506936u,        1.1463939726u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3195129567u,        1.1360122267u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3195129567u,        1.1360122267u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2787893569u,        1.1722745455u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2787893569u,        1.1722745455u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2922506936u,        1.1463939726u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2922506936u,        1.1463939726u) withcolor MidnightBlue;
draw (        0.2089220351u,        1.1040138520u)--(        0.1824528505u,        1.1162773712u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1562217028u,        1.1035126334u)--(        0.1824528505u,        1.1162773712u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2089220351u,        1.1040138520u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2089220351u,        1.1040138520u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1562217028u,        1.1035126334u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1562217028u,        1.1035126334u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1824528505u,        1.1162773712u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1824528505u,        1.1162773712u) withcolor MidnightBlue;
draw (        0.2089220351u,        1.1040138520u)--(        0.1824253222u,        1.1162177793u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2016957702u,        1.1381189544u)--(        0.1824253222u,        1.1162177793u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2089220351u,        1.1040138520u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2089220351u,        1.1040138520u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2016957702u,        1.1381189544u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2016957702u,        1.1381189544u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1824253222u,        1.1162177793u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1824253222u,        1.1162177793u) withcolor MidnightBlue;
draw (        0.2471931987u,        1.1275297971u)--(        0.2547855484u,        1.1556965837u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2787893569u,        1.1722745455u)--(        0.2547855484u,        1.1556965837u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2471931987u,        1.1275297971u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2471931987u,        1.1275297971u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2787893569u,        1.1722745455u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2787893569u,        1.1722745455u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2547855484u,        1.1556965837u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2547855484u,        1.1556965837u) withcolor MidnightBlue;
draw (        0.2471931987u,        1.1275297971u)--(        0.2284064196u,        1.1498472470u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2016957702u,        1.1381189544u)--(        0.2284064196u,        1.1498472470u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2471931987u,        1.1275297971u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2471931987u,        1.1275297971u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2016957702u,        1.1381189544u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2016957702u,        1.1381189544u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2284064196u,        1.1498472470u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2284064196u,        1.1498472470u) withcolor MidnightBlue;
draw (        0.2471931987u,        1.1275297971u)--(        0.2284195682u,        1.1498583088u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2360015261u,        1.1780278944u)--(        0.2284195682u,        1.1498583088u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2471931987u,        1.1275297971u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2471931987u,        1.1275297971u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2360015261u,        1.1780278944u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2360015261u,        1.1780278944u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2284195682u,        1.1498583088u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2284195682u,        1.1498583088u) withcolor MidnightBlue;
draw (        0.2471931987u,        1.1275297971u)--(        0.2547751566u,        1.1556993827u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2360015261u,        1.1780278944u)--(        0.2547751566u,        1.1556993827u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2471931987u,        1.1275297971u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2471931987u,        1.1275297971u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2360015261u,        1.1780278944u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2360015261u,        1.1780278944u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2547751566u,        1.1556993827u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2547751566u,        1.1556993827u) withcolor MidnightBlue;
draw (        0.2751812924u,        1.1227148179u)--(        0.2922327261u,        1.1463846350u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2787893569u,        1.1722745455u)--(        0.2922327261u,        1.1463846350u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2751812924u,        1.1227148179u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2751812924u,        1.1227148179u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2787893569u,        1.1722745455u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2787893569u,        1.1722745455u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2922327261u,        1.1463846350u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2922327261u,        1.1463846350u) withcolor MidnightBlue;
draw (        0.2787893569u,        1.1722745455u)--(        0.2547805084u,        1.1557038836u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2360015261u,        1.1780278944u)--(        0.2547805084u,        1.1557038836u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2787893569u,        1.1722745455u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2787893569u,        1.1722745455u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2360015261u,        1.1780278944u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2360015261u,        1.1780278944u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2547805084u,        1.1557038836u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2547805084u,        1.1557038836u) withcolor MidnightBlue;
draw (        0.1562217028u,        1.1035126334u)--(        0.1500000000u,        1.0750000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.1562217028u,        1.1035126334u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1562217028u,        1.1035126334u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1500000000u,        1.0750000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1500000000u,        1.0750000000u) withcolor Yellow;
draw (        0.1562217028u,        1.1035126334u)--(        0.1753944908u,        1.1254993530u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2016957702u,        1.1381189544u)--(        0.1753944908u,        1.1254993530u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.1562217028u,        1.1035126334u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1562217028u,        1.1035126334u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2016957702u,        1.1381189544u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2016957702u,        1.1381189544u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1753944908u,        1.1254993530u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1753944908u,        1.1254993530u) withcolor MidnightBlue;
draw (        0.1562217028u,        1.1035126334u)--(        0.1753809963u,        1.1255111131u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1556386735u,        1.1469878990u)--(        0.1753809963u,        1.1255111131u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.1562217028u,        1.1035126334u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1562217028u,        1.1035126334u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1556386735u,        1.1469878990u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1556386735u,        1.1469878990u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1753809963u,        1.1255111131u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1753809963u,        1.1255111131u) withcolor MidnightBlue;
draw (        0.2016957702u,        1.1381189544u)--(        0.2092973974u,        1.1662832386u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1879599508u,        1.1861760767u)--(        0.2092973974u,        1.1662832386u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2016957702u,        1.1381189544u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2016957702u,        1.1381189544u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1879599508u,        1.1861760767u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1879599508u,        1.1861760767u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2092973974u,        1.1662832386u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2092973974u,        1.1662832386u) withcolor MidnightBlue;
draw (        0.2016957702u,        1.1381189544u)--(        0.2092981578u,        1.1662830333u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2360015261u,        1.1780278944u)--(        0.2092981578u,        1.1662830333u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2016957702u,        1.1381189544u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2016957702u,        1.1381189544u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2360015261u,        1.1780278944u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2360015261u,        1.1780278944u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2092981578u,        1.1662830333u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2092981578u,        1.1662830333u) withcolor MidnightBlue;
draw (        0.2016957702u,        1.1381189544u)--(        0.2283991384u,        1.1498638155u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2360015261u,        1.1780278944u)--(        0.2283991384u,        1.1498638155u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2016957702u,        1.1381189544u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2016957702u,        1.1381189544u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2360015261u,        1.1780278944u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2360015261u,        1.1780278944u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2283991384u,        1.1498638155u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2283991384u,        1.1498638155u) withcolor MidnightBlue;
draw (        0.2016957702u,        1.1381189544u)--(        0.1753864491u,        1.1255161269u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1556386735u,        1.1469878990u)--(        0.1753864491u,        1.1255161269u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2016957702u,        1.1381189544u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2016957702u,        1.1381189544u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1556386735u,        1.1469878990u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1556386735u,        1.1469878990u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1753864491u,        1.1255161269u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1753864491u,        1.1255161269u) withcolor MidnightBlue;
draw (        0.1879599508u,        1.1861760767u)--(        0.2092978527u,        1.1662837270u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2360015261u,        1.1780278944u)--(        0.2092978527u,        1.1662837270u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.1879599508u,        1.1861760767u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1879599508u,        1.1861760767u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2360015261u,        1.1780278944u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2360015261u,        1.1780278944u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2092978527u,        1.1662837270u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2092978527u,        1.1662837270u) withcolor MidnightBlue;
% r =         0.0291721038
endfig;
end;
