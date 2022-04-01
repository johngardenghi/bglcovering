prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(80);
u := 20cm;
path pols[];
pair p[];
R := .0327661960u;
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

x1 = .5959549757u;
y1 = .3290531539u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .5994457994u;
y2 = .2850932442u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .5948295150u;
y3 = .3718948199u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .5807249799u;
y4 = .4235999271u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .6359517930u;
y5 = .4068251184u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .5829373034u;
y6 = .4841789766u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .6256437405u;
y7 = .4594469490u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .6797845349u;
y8 = .4468871993u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .5441414580u;
y9 = .5752474700u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = .5868904284u;
y10 = .5417897375u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;

x11 = .6320145375u;
y11 = .5154658125u;
draw fullcircle scaled D shifted z11 withcolor RoyalBlue withpen pencircle scaled 1.0;

x12 = .6772417869u;
y12 = .4895049303u;
draw fullcircle scaled D shifted z12 withcolor RoyalBlue withpen pencircle scaled 1.0;

x13 = .5500456829u;
y13 = .6279452450u;
draw fullcircle scaled D shifted z13 withcolor RoyalBlue withpen pencircle scaled 1.0;

x14 = .5972561614u;
y14 = .5977648746u;
draw fullcircle scaled D shifted z14 withcolor RoyalBlue withpen pencircle scaled 1.0;

x15 = .6399130142u;
y15 = .5704583001u;
draw fullcircle scaled D shifted z15 withcolor RoyalBlue withpen pencircle scaled 1.0;

x16 = .6850603060u;
y16 = .5441534337u;
draw fullcircle scaled D shifted z16 withcolor RoyalBlue withpen pencircle scaled 1.0;

x17 = .7301044050u;
y17 = .5080304002u;
draw fullcircle scaled D shifted z17 withcolor RoyalBlue withpen pencircle scaled 1.0;

x18 = .5604056123u;
y18 = .6750493493u;
draw fullcircle scaled D shifted z18 withcolor RoyalBlue withpen pencircle scaled 1.0;

x19 = .6069702582u;
y19 = .6479925845u;
draw fullcircle scaled D shifted z19 withcolor RoyalBlue withpen pencircle scaled 1.0;

x20 = .6542492980u;
y20 = .6220558209u;
draw fullcircle scaled D shifted z20 withcolor RoyalBlue withpen pencircle scaled 1.0;

x21 = .6973303197u;
y21 = .5889928732u;
draw fullcircle scaled D shifted z21 withcolor RoyalBlue withpen pencircle scaled 1.0;

x22 = .7447347992u;
y22 = .5527306272u;
draw fullcircle scaled D shifted z22 withcolor RoyalBlue withpen pencircle scaled 1.0;

x23 = .3973273624u;
y23 = .8110453700u;
draw fullcircle scaled D shifted z23 withcolor RoyalBlue withpen pencircle scaled 1.0;

x24 = .5092853294u;
y24 = .7216350941u;
draw fullcircle scaled D shifted z24 withcolor RoyalBlue withpen pencircle scaled 1.0;

x25 = .5691887156u;
y25 = .7130840607u;
draw fullcircle scaled D shifted z25 withcolor RoyalBlue withpen pencircle scaled 1.0;

x26 = .6169783122u;
y26 = .6991534148u;
draw fullcircle scaled D shifted z26 withcolor RoyalBlue withpen pencircle scaled 1.0;

x27 = .6634911322u;
y27 = .6728333922u;
draw fullcircle scaled D shifted z27 withcolor RoyalBlue withpen pencircle scaled 1.0;

x28 = .7106407895u;
y28 = .6408619670u;
draw fullcircle scaled D shifted z28 withcolor RoyalBlue withpen pencircle scaled 1.0;

x29 = .7532992303u;
y29 = .6026936205u;
draw fullcircle scaled D shifted z29 withcolor RoyalBlue withpen pencircle scaled 1.0;

x30 = .3581493751u;
y30 = .8476161909u;
draw fullcircle scaled D shifted z30 withcolor RoyalBlue withpen pencircle scaled 1.0;

x31 = .4194193297u;
y31 = .8549421128u;
draw fullcircle scaled D shifted z31 withcolor RoyalBlue withpen pencircle scaled 1.0;

x32 = .4319146865u;
y32 = .7673553241u;
draw fullcircle scaled D shifted z32 withcolor RoyalBlue withpen pencircle scaled 1.0;

x33 = .4850154989u;
y33 = .7693396266u;
draw fullcircle scaled D shifted z33 withcolor RoyalBlue withpen pencircle scaled 1.0;

x34 = .3353521029u;
y34 = .9002651245u;
draw fullcircle scaled D shifted z34 withcolor RoyalBlue withpen pencircle scaled 1.0;

x35 = .3927963502u;
y35 = .8892634619u;
draw fullcircle scaled D shifted z35 withcolor RoyalBlue withpen pencircle scaled 1.0;

x36 = .4542636452u;
y36 = .8802872608u;
draw fullcircle scaled D shifted z36 withcolor RoyalBlue withpen pencircle scaled 1.0;

x37 = .5107928238u;
y37 = .8712334325u;
draw fullcircle scaled D shifted z37 withcolor RoyalBlue withpen pencircle scaled 1.0;

x38 = .5411745612u;
y38 = .8183394663u;
draw fullcircle scaled D shifted z38 withcolor RoyalBlue withpen pencircle scaled 1.0;

x39 = .3098599580u;
y39 = .9496955370u;
draw fullcircle scaled D shifted z39 withcolor RoyalBlue withpen pencircle scaled 1.0;

x40 = .3695113495u;
y40 = .9419387261u;
draw fullcircle scaled D shifted z40 withcolor RoyalBlue withpen pencircle scaled 1.0;

x41 = .4276913126u;
y41 = .9277811004u;
draw fullcircle scaled D shifted z41 withcolor RoyalBlue withpen pencircle scaled 1.0;

x42 = .4876939035u;
y42 = .9241670814u;
draw fullcircle scaled D shifted z42 withcolor RoyalBlue withpen pencircle scaled 1.0;

x43 = .5459180362u;
y43 = .9081965972u;
draw fullcircle scaled D shifted z43 withcolor RoyalBlue withpen pencircle scaled 1.0;

x44 = .1648660067u;
y44 = 1.1042035197u;
draw fullcircle scaled D shifted z44 withcolor RoyalBlue withpen pencircle scaled 1.0;

x45 = .3418850762u;
y45 = .9874176953u;
draw fullcircle scaled D shifted z45 withcolor RoyalBlue withpen pencircle scaled 1.0;

x46 = .4042835061u;
y46 = .9804665101u;
draw fullcircle scaled D shifted z46 withcolor RoyalBlue withpen pencircle scaled 1.0;

x47 = .4615905712u;
y47 = .9699691750u;
draw fullcircle scaled D shifted z47 withcolor RoyalBlue withpen pencircle scaled 1.0;

x48 = .5240395312u;
y48 = .9610077011u;
draw fullcircle scaled D shifted z48 withcolor RoyalBlue withpen pencircle scaled 1.0;

x49 = .5782873344u;
y49 = .9569336314u;
draw fullcircle scaled D shifted z49 withcolor RoyalBlue withpen pencircle scaled 1.0;

x50 = .1651731031u;
y50 = 1.1372176861u;
draw fullcircle scaled D shifted z50 withcolor RoyalBlue withpen pencircle scaled 1.0;

x51 = .2212106333u;
y51 = 1.1208290280u;
draw fullcircle scaled D shifted z51 withcolor RoyalBlue withpen pencircle scaled 1.0;

x52 = .3066548997u;
y52 = 1.0808615973u;
draw fullcircle scaled D shifted z52 withcolor RoyalBlue withpen pencircle scaled 1.0;

x53 = .3313781856u;
y53 = 1.0376158828u;
draw fullcircle scaled D shifted z53 withcolor RoyalBlue withpen pencircle scaled 1.0;

x54 = .3873494378u;
y54 = 1.0232744690u;
draw fullcircle scaled D shifted z54 withcolor RoyalBlue withpen pencircle scaled 1.0;

x55 = .4449834019u;
y55 = 1.0213896339u;
draw fullcircle scaled D shifted z55 withcolor RoyalBlue withpen pencircle scaled 1.0;

x56 = .5009532695u;
y56 = 1.0063115618u;
draw fullcircle scaled D shifted z56 withcolor RoyalBlue withpen pencircle scaled 1.0;

x57 = .5631400122u;
y57 = .9990242585u;
draw fullcircle scaled D shifted z57 withcolor RoyalBlue withpen pencircle scaled 1.0;

x58 = .6257882113u;
y58 = .9899853315u;
draw fullcircle scaled D shifted z58 withcolor RoyalBlue withpen pencircle scaled 1.0;

x59 = .1638594049u;
y59 = 1.1749621321u;
draw fullcircle scaled D shifted z59 withcolor RoyalBlue withpen pencircle scaled 1.0;

x60 = .2147377429u;
y60 = 1.1711154866u;
draw fullcircle scaled D shifted z60 withcolor RoyalBlue withpen pencircle scaled 1.0;

x61 = .2699127465u;
y61 = 1.1229604776u;
draw fullcircle scaled D shifted z61 withcolor RoyalBlue withpen pencircle scaled 1.0;

x62 = .3297577457u;
y62 = 1.1292033149u;
draw fullcircle scaled D shifted z62 withcolor RoyalBlue withpen pencircle scaled 1.0;

x63 = .3621204243u;
y63 = 1.0781743644u;
draw fullcircle scaled D shifted z63 withcolor RoyalBlue withpen pencircle scaled 1.0;

x64 = .4155877989u;
y64 = 1.0706490875u;
draw fullcircle scaled D shifted z64 withcolor RoyalBlue withpen pencircle scaled 1.0;

x65 = .4738750965u;
y65 = 1.0612455013u;
draw fullcircle scaled D shifted z65 withcolor RoyalBlue withpen pencircle scaled 1.0;

x66 = .5162234888u;
y66 = 1.0396706114u;
draw fullcircle scaled D shifted z66 withcolor RoyalBlue withpen pencircle scaled 1.0;

x67 = .5946878703u;
y67 = 1.0356492141u;
draw fullcircle scaled D shifted z67 withcolor RoyalBlue withpen pencircle scaled 1.0;

x68 = .6573725793u;
y68 = 1.0265913456u;
draw fullcircle scaled D shifted z68 withcolor RoyalBlue withpen pencircle scaled 1.0;

x69 = .2545346837u;
y69 = 1.1750031374u;
draw fullcircle scaled D shifted z69 withcolor RoyalBlue withpen pencircle scaled 1.0;

x70 = .3084077687u;
y70 = 1.1702376677u;
draw fullcircle scaled D shifted z70 withcolor RoyalBlue withpen pencircle scaled 1.0;

x71 = .3657710217u;
y71 = 1.1718826150u;
draw fullcircle scaled D shifted z71 withcolor RoyalBlue withpen pencircle scaled 1.0;

x72 = .3877668800u;
y72 = 1.1207596090u;
draw fullcircle scaled D shifted z72 withcolor RoyalBlue withpen pencircle scaled 1.0;

x73 = .4469859884u;
y73 = 1.1127520740u;
draw fullcircle scaled D shifted z73 withcolor RoyalBlue withpen pencircle scaled 1.0;

x74 = .5001059769u;
y74 = 1.1136541706u;
draw fullcircle scaled D shifted z74 withcolor RoyalBlue withpen pencircle scaled 1.0;

x75 = .5809355355u;
y75 = 1.0885043402u;
draw fullcircle scaled D shifted z75 withcolor RoyalBlue withpen pencircle scaled 1.0;

x76 = .6272411329u;
y76 = 1.0722130697u;
draw fullcircle scaled D shifted z76 withcolor RoyalBlue withpen pencircle scaled 1.0;

x77 = .4240762185u;
y77 = 1.1624814437u;
draw fullcircle scaled D shifted z77 withcolor RoyalBlue withpen pencircle scaled 1.0;

x78 = .4824832725u;
y78 = 1.1568996303u;
draw fullcircle scaled D shifted z78 withcolor RoyalBlue withpen pencircle scaled 1.0;

x79 = .5428526001u;
y79 = 1.1565053156u;
draw fullcircle scaled D shifted z79 withcolor RoyalBlue withpen pencircle scaled 1.0;

x80 = .5406996162u;
y80 = 1.1987267351u;
draw fullcircle scaled D shifted z80 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.3581493751u,        0.8476161909u)--(        0.3400000000u,        0.8750000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3581493751u,        0.8476161909u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3581493751u,        0.8476161909u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3400000000u,        0.8750000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3400000000u,        0.8750000000u) withcolor Yellow;
draw (        0.3581493751u,        0.8476161909u)--(        0.3550000000u,        0.8150000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3581493751u,        0.8476161909u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3581493751u,        0.8476161909u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3550000000u,        0.8150000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3550000000u,        0.8150000000u) withcolor Yellow;
draw (        0.3581493751u,        0.8476161909u)--(        0.3600000000u,        0.8150000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3581493751u,        0.8476161909u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3581493751u,        0.8476161909u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3600000000u,        0.8150000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3600000000u,        0.8150000000u) withcolor Yellow;
draw (        0.3581493751u,        0.8476161909u)--(        0.3612810985u,        0.8802323817u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3353521029u,        0.9002651245u)--(        0.3612810985u,        0.8802323817u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3581493751u,        0.8476161909u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3581493751u,        0.8476161909u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3353521029u,        0.9002651245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3353521029u,        0.9002651245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3612810985u,        0.8802323817u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3612810985u,        0.8802323817u) withcolor MidnightBlue;
draw (        0.3581493751u,        0.8476161909u)--(        0.3612998330u,        0.8802305775u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3927963502u,        0.8892634619u)--(        0.3612998330u,        0.8802305775u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3581493751u,        0.8476161909u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3581493751u,        0.8476161909u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3927963502u,        0.8892634619u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3927963502u,        0.8892634619u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3612998330u,        0.8802305775u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3612998330u,        0.8802305775u) withcolor MidnightBlue;
draw (        0.3353521029u,        0.9002651245u)--(        0.3668540671u,        0.9092789945u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3927963502u,        0.8892634619u)--(        0.3668540671u,        0.9092789945u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3353521029u,        0.9002651245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3353521029u,        0.9002651245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3927963502u,        0.8892634619u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3927963502u,        0.8892634619u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3668540671u,        0.9092789945u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3668540671u,        0.9092789945u) withcolor MidnightBlue;
draw (        0.3353521029u,        0.9002651245u)--(        0.3612943860u,        0.8802495919u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3927963502u,        0.8892634619u)--(        0.3612943860u,        0.8802495919u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3353521029u,        0.9002651245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3353521029u,        0.9002651245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3927963502u,        0.8892634619u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3927963502u,        0.8892634619u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3612943860u,        0.8802495919u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3612943860u,        0.8802495919u) withcolor MidnightBlue;
draw (        0.3353521029u,        0.9002651245u)--(        0.3380081740u,        0.9329234908u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3098599580u,        0.9496955370u)--(        0.3380081740u,        0.9329234908u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3353521029u,        0.9002651245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3353521029u,        0.9002651245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3098599580u,        0.9496955370u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3098599580u,        0.9496955370u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3380081740u,        0.9329234908u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3380081740u,        0.9329234908u) withcolor MidnightBlue;
draw (        0.3353521029u,        0.9002651245u)--(        0.3380098142u,        0.9329233573u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3695113495u,        0.9419387261u)--(        0.3380098142u,        0.9329233573u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3353521029u,        0.9002651245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3353521029u,        0.9002651245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3695113495u,        0.9419387261u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3695113495u,        0.9419387261u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3380098142u,        0.9329233573u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3380098142u,        0.9329233573u) withcolor MidnightBlue;
draw (        0.3353521029u,        0.9002651245u)--(        0.3668536382u,        0.9092804933u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3695113495u,        0.9419387261u)--(        0.3668536382u,        0.9092804933u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3353521029u,        0.9002651245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3353521029u,        0.9002651245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3695113495u,        0.9419387261u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3695113495u,        0.9419387261u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3668536382u,        0.9092804933u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3668536382u,        0.9092804933u) withcolor MidnightBlue;
draw (        0.3927963502u,        0.8892634619u)--(        0.4250380873u,        0.8951024678u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4542636452u,        0.8802872608u)--(        0.4250380873u,        0.8951024678u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3927963502u,        0.8892634619u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3927963502u,        0.8892634619u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4542636452u,        0.8802872608u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4542636452u,        0.8802872608u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4250380873u,        0.8951024678u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4250380873u,        0.8951024678u) withcolor MidnightBlue;
draw (        0.3927963502u,        0.8892634619u)--(        0.3668551299u,        0.9092803720u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3695113495u,        0.9419387261u)--(        0.3668551299u,        0.9092803720u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3927963502u,        0.8892634619u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3927963502u,        0.8892634619u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3695113495u,        0.9419387261u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3695113495u,        0.9419387261u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3668551299u,        0.9092803720u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3668551299u,        0.9092803720u) withcolor MidnightBlue;
draw (        0.3927963502u,        0.8892634619u)--(        0.3954525697u,        0.9219218161u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3695113495u,        0.9419387261u)--(        0.3954525697u,        0.9219218161u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3927963502u,        0.8892634619u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3927963502u,        0.8892634619u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3695113495u,        0.9419387261u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3695113495u,        0.9419387261u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3954525697u,        0.9219218161u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3954525697u,        0.9219218161u) withcolor MidnightBlue;
draw (        0.3927963502u,        0.8892634619u)--(        0.3954532649u,        0.9219217595u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4276913126u,        0.9277811004u)--(        0.3954532649u,        0.9219217595u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3927963502u,        0.8892634619u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3927963502u,        0.8892634619u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4276913126u,        0.9277811004u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4276913126u,        0.9277811004u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3954532649u,        0.9219217595u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3954532649u,        0.9219217595u) withcolor MidnightBlue;
draw (        0.3927963502u,        0.8892634619u)--(        0.4250343979u,        0.8951228028u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4276913126u,        0.9277811004u)--(        0.4250343979u,        0.8951228028u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3927963502u,        0.8892634619u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3927963502u,        0.8892634619u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4276913126u,        0.9277811004u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4276913126u,        0.9277811004u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4250343979u,        0.8951228028u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4250343979u,        0.8951228028u) withcolor MidnightBlue;
draw (        0.4542636452u,        0.8802872608u)--(        0.4850498835u,        0.8915047084u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5107928238u,        0.8712334325u)--(        0.4850498835u,        0.8915047084u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4542636452u,        0.8802872608u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4542636452u,        0.8802872608u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5107928238u,        0.8712334325u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5107928238u,        0.8712334325u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4850498835u,        0.8915047084u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4850498835u,        0.8915047084u) withcolor MidnightBlue;
draw (        0.4542636452u,        0.8802872608u)--(        0.4250478502u,        0.8951217112u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4276913126u,        0.9277811004u)--(        0.4250478502u,        0.8951217112u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4542636452u,        0.8802872608u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4542636452u,        0.8802872608u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4276913126u,        0.9277811004u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4276913126u,        0.9277811004u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4250478502u,        0.8951217112u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4250478502u,        0.8951217112u) withcolor MidnightBlue;
draw (        0.4542636452u,        0.8802872608u)--(        0.4569071077u,        0.9129466500u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4276913126u,        0.9277811004u)--(        0.4569071077u,        0.9129466500u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4542636452u,        0.8802872608u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4542636452u,        0.8802872608u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4276913126u,        0.9277811004u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4276913126u,        0.9277811004u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4569071077u,        0.9129466500u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4569071077u,        0.9129466500u) withcolor MidnightBlue;
draw (        0.4542636452u,        0.8802872608u)--(        0.4569088026u,        0.9129465128u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4876939035u,        0.9241670814u)--(        0.4569088026u,        0.9129465128u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4542636452u,        0.8802872608u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4542636452u,        0.8802872608u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4876939035u,        0.9241670814u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4876939035u,        0.9241670814u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4569088026u,        0.9129465128u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4569088026u,        0.9129465128u) withcolor MidnightBlue;
draw (        0.4542636452u,        0.8802872608u)--(        0.4850487461u,        0.8915078294u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4876939035u,        0.9241670814u)--(        0.4850487461u,        0.8915078294u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4542636452u,        0.8802872608u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4542636452u,        0.8802872608u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4876939035u,        0.9241670814u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4876939035u,        0.9241670814u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4850487461u,        0.8915078294u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4850487461u,        0.8915078294u) withcolor MidnightBlue;
draw (        0.5107928238u,        0.8712334325u)--(        0.4850521261u,        0.8915075559u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4876939035u,        0.9241670814u)--(        0.4850521261u,        0.8915075559u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5107928238u,        0.8712334325u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5107928238u,        0.8712334325u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4876939035u,        0.9241670814u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4876939035u,        0.9241670814u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4850521261u,        0.8915075559u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4850521261u,        0.8915075559u) withcolor MidnightBlue;
draw (        0.5107928238u,        0.8712334325u)--(        0.5134346012u,        0.9038929580u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4876939035u,        0.9241670814u)--(        0.5134346012u,        0.9038929580u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5107928238u,        0.8712334325u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5107928238u,        0.8712334325u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4876939035u,        0.9241670814u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4876939035u,        0.9241670814u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5134346012u,        0.9038929580u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5134346012u,        0.9038929580u) withcolor MidnightBlue;
draw (        0.5107928238u,        0.8712334325u)--(        0.5134357100u,        0.9038928683u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5459180362u,        0.9081965972u)--(        0.5134357100u,        0.9038928683u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5107928238u,        0.8712334325u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5107928238u,        0.8712334325u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5459180362u,        0.9081965972u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5459180362u,        0.9081965972u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5134357100u,        0.9038928683u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5134357100u,        0.9038928683u) withcolor MidnightBlue;
draw (        0.3098599580u,        0.9496955370u)--(        0.3380092562u,        0.9329253072u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3695113495u,        0.9419387261u)--(        0.3380092562u,        0.9329253072u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3098599580u,        0.9496955370u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3098599580u,        0.9496955370u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3695113495u,        0.9419387261u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3695113495u,        0.9419387261u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3380092562u,        0.9329253072u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3380092562u,        0.9329253072u) withcolor MidnightBlue;
draw (        0.3695113495u,        0.9419387261u)--(        0.4017495380u,        0.9477972921u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4276913126u,        0.9277811004u)--(        0.4017495380u,        0.9477972921u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3695113495u,        0.9419387261u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3695113495u,        0.9419387261u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4276913126u,        0.9277811004u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4276913126u,        0.9277811004u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4017495380u,        0.9477972921u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4017495380u,        0.9477972921u) withcolor MidnightBlue;
draw (        0.3695113495u,        0.9419387261u)--(        0.3954531241u,        0.9219225344u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4276913126u,        0.9277811004u)--(        0.3954531241u,        0.9219225344u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3695113495u,        0.9419387261u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3695113495u,        0.9419387261u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4276913126u,        0.9277811004u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4276913126u,        0.9277811004u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3954531241u,        0.9219225344u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3954531241u,        0.9219225344u) withcolor MidnightBlue;
draw (        0.3695113495u,        0.9419387261u)--(        0.3720431384u,        0.9746069621u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3418850762u,        0.9874176953u)--(        0.3720431384u,        0.9746069621u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3695113495u,        0.9419387261u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3695113495u,        0.9419387261u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3418850762u,        0.9874176953u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3418850762u,        0.9874176953u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3720431384u,        0.9746069621u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3720431384u,        0.9746069621u) withcolor MidnightBlue;
draw (        0.3695113495u,        0.9419387261u)--(        0.3720455297u,        0.9746067767u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4042835061u,        0.9804665101u)--(        0.3720455297u,        0.9746067767u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3695113495u,        0.9419387261u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3695113495u,        0.9419387261u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4042835061u,        0.9804665101u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4042835061u,        0.9804665101u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3720455297u,        0.9746067767u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3720455297u,        0.9746067767u) withcolor MidnightBlue;
draw (        0.3695113495u,        0.9419387261u)--(        0.4017493259u,        0.9477984596u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4042835061u,        0.9804665101u)--(        0.4017493259u,        0.9477984596u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3695113495u,        0.9419387261u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3695113495u,        0.9419387261u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4042835061u,        0.9804665101u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4042835061u,        0.9804665101u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4017493259u,        0.9477984596u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4017493259u,        0.9477984596u) withcolor MidnightBlue;
draw (        0.4276913126u,        0.9277811004u)--(        0.4569080653u,        0.9129485361u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4876939035u,        0.9241670814u)--(        0.4569080653u,        0.9129485361u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4276913126u,        0.9277811004u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4276913126u,        0.9277811004u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4876939035u,        0.9241670814u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4876939035u,        0.9241670814u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4569080653u,        0.9129485361u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4569080653u,        0.9129485361u) withcolor MidnightBlue;
draw (        0.4276913126u,        0.9277811004u)--(        0.4017503760u,        0.9477983781u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4042835061u,        0.9804665101u)--(        0.4017503760u,        0.9477983781u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4276913126u,        0.9277811004u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4276913126u,        0.9277811004u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4042835061u,        0.9804665101u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4042835061u,        0.9804665101u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4017503760u,        0.9477983781u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4017503760u,        0.9477983781u) withcolor MidnightBlue;
draw (        0.4276913126u,        0.9277811004u)--(        0.4302244427u,        0.9604492324u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4042835061u,        0.9804665101u)--(        0.4302244427u,        0.9604492324u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4276913126u,        0.9277811004u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4276913126u,        0.9277811004u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4042835061u,        0.9804665101u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4042835061u,        0.9804665101u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4302244427u,        0.9604492324u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4302244427u,        0.9604492324u) withcolor MidnightBlue;
draw (        0.4276913126u,        0.9277811004u)--(        0.4302381552u,        0.9604481663u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4615905712u,        0.9699691750u)--(        0.4302381552u,        0.9604481663u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4276913126u,        0.9277811004u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4276913126u,        0.9277811004u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4615905712u,        0.9699691750u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4615905712u,        0.9699691750u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4302381552u,        0.9604481663u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4302381552u,        0.9604481663u) withcolor MidnightBlue;
draw (        0.4876939035u,        0.9241670814u)--(        0.5201763996u,        0.9284695285u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5459180362u,        0.9081965972u)--(        0.5201763996u,        0.9284695285u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4876939035u,        0.9241670814u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4876939035u,        0.9241670814u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5459180362u,        0.9081965972u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5459180362u,        0.9081965972u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5201763996u,        0.9284695285u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5201763996u,        0.9284695285u) withcolor MidnightBlue;
draw (        0.4876939035u,        0.9241670814u)--(        0.5134355402u,        0.9038941502u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5459180362u,        0.9081965972u)--(        0.5134355402u,        0.9038941502u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4876939035u,        0.9241670814u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4876939035u,        0.9241670814u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5459180362u,        0.9081965972u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5459180362u,        0.9081965972u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5134355402u,        0.9038941502u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5134355402u,        0.9038941502u) withcolor MidnightBlue;
draw (        0.4876939035u,        0.9241670814u)--(        0.4915521147u,        0.9567053328u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4615905712u,        0.9699691750u)--(        0.4915521147u,        0.9567053328u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4876939035u,        0.9241670814u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4876939035u,        0.9241670814u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4615905712u,        0.9699691750u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4615905712u,        0.9699691750u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4915521147u,        0.9567053328u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4915521147u,        0.9567053328u) withcolor MidnightBlue;
draw (        0.4876939035u,        0.9241670814u)--(        0.4915571032u,        0.9567047409u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5240395312u,        0.9610077011u)--(        0.4915571032u,        0.9567047409u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4876939035u,        0.9241670814u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4876939035u,        0.9241670814u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5240395312u,        0.9610077011u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5240395312u,        0.9610077011u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4915571032u,        0.9567047409u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4915571032u,        0.9567047409u) withcolor MidnightBlue;
draw (        0.4876939035u,        0.9241670814u)--(        0.5201763316u,        0.9284700416u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5240395312u,        0.9610077011u)--(        0.5201763316u,        0.9284700416u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4876939035u,        0.9241670814u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4876939035u,        0.9241670814u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5240395312u,        0.9610077011u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5240395312u,        0.9610077011u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5201763316u,        0.9284700416u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5201763316u,        0.9284700416u) withcolor MidnightBlue;
draw (        0.5459180362u,        0.9081965972u)--(        0.5201767633u,        0.9284699903u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5240395312u,        0.9610077011u)--(        0.5201767633u,        0.9284699903u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5459180362u,        0.9081965972u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5459180362u,        0.9081965972u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5240395312u,        0.9610077011u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5240395312u,        0.9610077011u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5201767633u,        0.9284699903u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5201767633u,        0.9284699903u) withcolor MidnightBlue;
draw (        0.5459180362u,        0.9081965972u)--(        0.5497808041u,        0.9407343080u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5240395312u,        0.9610077011u)--(        0.5497808041u,        0.9407343080u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5459180362u,        0.9081965972u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5459180362u,        0.9081965972u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5240395312u,        0.9610077011u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5240395312u,        0.9610077011u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5497808041u,        0.9407343080u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5497808041u,        0.9407343080u) withcolor MidnightBlue;
draw (        0.5459180362u,        0.9081965972u)--(        0.5498074717u,        0.9407311310u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5782873344u,        0.9569336314u)--(        0.5498074717u,        0.9407311310u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5459180362u,        0.9081965972u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5459180362u,        0.9081965972u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5782873344u,        0.9569336314u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5782873344u,        0.9569336314u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5498074717u,        0.9407311310u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5498074717u,        0.9407311310u) withcolor MidnightBlue;
draw (        0.5459180362u,        0.9081965972u)--(        0.5743978989u,        0.9243990976u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5782873344u,        0.9569336314u)--(        0.5743978989u,        0.9243990976u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5459180362u,        0.9081965972u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5459180362u,        0.9081965972u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5782873344u,        0.9569336314u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5782873344u,        0.9569336314u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5743978989u,        0.9243990976u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5743978989u,        0.9243990976u) withcolor Green;
draw (        0.1648660067u,        1.1042035197u)--(        0.1500000000u,        1.0750000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.1648660067u,        1.1042035197u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1648660067u,        1.1042035197u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1500000000u,        1.0750000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1500000000u,        1.0750000000u) withcolor Yellow;
draw (        0.3418850762u,        0.9874176953u)--(        0.3741237927u,        0.9932733555u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4042835061u,        0.9804665101u)--(        0.3741237927u,        0.9932733555u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3418850762u,        0.9874176953u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3418850762u,        0.9874176953u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4042835061u,        0.9804665101u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4042835061u,        0.9804665101u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3741237927u,        0.9932733555u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3741237927u,        0.9932733555u) withcolor MidnightBlue;
draw (        0.3418850762u,        0.9874176953u)--(        0.3720447896u,        0.9746108499u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4042835061u,        0.9804665101u)--(        0.3720447896u,        0.9746108499u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3418850762u,        0.9874176953u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3418850762u,        0.9874176953u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4042835061u,        0.9804665101u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4042835061u,        0.9804665101u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3720447896u,        0.9746108499u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3720447896u,        0.9746108499u) withcolor MidnightBlue;
draw (        0.3418850762u,        0.9874176953u)--(        0.3166667150u,        1.0083379691u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3313781856u,        1.0376158828u)--(        0.3166667150u,        1.0083379691u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3418850762u,        0.9874176953u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3418850762u,        0.9874176953u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3313781856u,        1.0376158828u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3313781856u,        1.0376158828u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3166667150u,        1.0083379691u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3166667150u,        1.0083379691u) withcolor Green;
draw (        0.3418850762u,        0.9874176953u)--(        0.3741193019u,        0.9932980266u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3873494378u,        1.0232744690u)--(        0.3741193019u,        0.9932980266u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3418850762u,        0.9874176953u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3418850762u,        0.9874176953u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3873494378u,        1.0232744690u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3873494378u,        1.0232744690u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3741193019u,        0.9932980266u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3741193019u,        0.9932980266u) withcolor MidnightBlue;
draw (        0.4042835061u,        0.9804665101u)--(        0.4356400260u,        0.9899739945u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4615905712u,        0.9699691750u)--(        0.4356400260u,        0.9899739945u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4042835061u,        0.9804665101u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4042835061u,        0.9804665101u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4615905712u,        0.9699691750u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4615905712u,        0.9699691750u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4356400260u,        0.9899739945u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4356400260u,        0.9899739945u) withcolor MidnightBlue;
draw (        0.4042835061u,        0.9804665101u)--(        0.4302340513u,        0.9604616907u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4615905712u,        0.9699691750u)--(        0.4302340513u,        0.9604616907u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4042835061u,        0.9804665101u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4042835061u,        0.9804665101u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4615905712u,        0.9699691750u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4615905712u,        0.9699691750u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4302340513u,        0.9604616907u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4302340513u,        0.9604616907u) withcolor MidnightBlue;
draw (        0.4042835061u,        0.9804665101u)--(        0.3741319138u,        0.9932924635u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3873494378u,        1.0232744690u)--(        0.3741319138u,        0.9932924635u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4042835061u,        0.9804665101u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4042835061u,        0.9804665101u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3873494378u,        1.0232744690u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3873494378u,        1.0232744690u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3741319138u,        0.9932924635u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3741319138u,        0.9932924635u) withcolor MidnightBlue;
draw (        0.4042835061u,        0.9804665101u)--(        0.4356367480u,        0.9899847989u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4449834019u,        1.0213896339u)--(        0.4356367480u,        0.9899847989u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4042835061u,        0.9804665101u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4042835061u,        0.9804665101u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4449834019u,        1.0213896339u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4449834019u,        1.0213896339u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4356367480u,        0.9899847989u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4356367480u,        0.9899847989u) withcolor MidnightBlue;
draw (        0.4615905712u,        0.9699691750u)--(        0.4940742079u,        0.9742630017u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5240395312u,        0.9610077011u)--(        0.4940742079u,        0.9742630017u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4615905712u,        0.9699691750u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4615905712u,        0.9699691750u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5240395312u,        0.9610077011u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5240395312u,        0.9610077011u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4940742079u,        0.9742630017u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4940742079u,        0.9742630017u) withcolor MidnightBlue;
draw (        0.4615905712u,        0.9699691750u)--(        0.4915558946u,        0.9567138744u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5240395312u,        0.9610077011u)--(        0.4915558946u,        0.9567138744u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4615905712u,        0.9699691750u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4615905712u,        0.9699691750u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5240395312u,        0.9610077011u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5240395312u,        0.9610077011u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4915558946u,        0.9567138744u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4915558946u,        0.9567138744u) withcolor MidnightBlue;
draw (        0.4615905712u,        0.9699691750u)--(        0.4356461915u,        0.9899819899u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4449834019u,        1.0213896339u)--(        0.4356461915u,        0.9899819899u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4615905712u,        0.9699691750u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4615905712u,        0.9699691750u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4449834019u,        1.0213896339u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4449834019u,        1.0213896339u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4356461915u,        0.9899819899u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4356461915u,        0.9899819899u) withcolor MidnightBlue;
draw (        0.4615905712u,        0.9699691750u)--(        0.4940724896u,        0.9742759806u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5009532695u,        1.0063115618u)--(        0.4940724896u,        0.9742759806u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4615905712u,        0.9699691750u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4615905712u,        0.9699691750u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5009532695u,        1.0063115618u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5009532695u,        1.0063115618u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4940724896u,        0.9742759806u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4940724896u,        0.9742759806u) withcolor MidnightBlue;
draw (        0.5240395312u,        0.9610077011u)--(        0.5497952398u,        0.9407526505u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5782873344u,        0.9569336314u)--(        0.5497952398u,        0.9407526505u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5240395312u,        0.9610077011u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5240395312u,        0.9610077011u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5782873344u,        0.9569336314u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5782873344u,        0.9569336314u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5497952398u,        0.9407526505u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5497952398u,        0.9407526505u) withcolor MidnightBlue;
draw (        0.5240395312u,        0.9610077011u)--(        0.4940793046u,        0.9742745176u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5009532695u,        1.0063115618u)--(        0.4940793046u,        0.9742745176u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5240395312u,        0.9610077011u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5240395312u,        0.9610077011u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5009532695u,        1.0063115618u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5009532695u,        1.0063115618u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4940793046u,        0.9742745176u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4940793046u,        0.9742745176u) withcolor MidnightBlue;
draw (        0.5240395312u,        0.9610077011u)--(        0.5309134961u,        0.9930447453u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5009532695u,        1.0063115618u)--(        0.5309134961u,        0.9930447453u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5240395312u,        0.9610077011u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5240395312u,        0.9610077011u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5009532695u,        1.0063115618u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5009532695u,        1.0063115618u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5309134961u,        0.9930447453u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5309134961u,        0.9930447453u) withcolor MidnightBlue;
draw (        0.5240395312u,        0.9610077011u)--(        0.5309244752u,        0.9930423877u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5631400122u,        0.9990242585u)--(        0.5309244752u,        0.9930423877u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5240395312u,        0.9610077011u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5240395312u,        0.9610077011u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5631400122u,        0.9990242585u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5631400122u,        0.9990242585u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5309244752u,        0.9930423877u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5309244752u,        0.9930423877u) withcolor MidnightBlue;
draw (        0.5782873344u,        0.9569336314u)--(        0.5932440289u,        0.9860870341u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5631400122u,        0.9990242585u)--(        0.5932440289u,        0.9860870341u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5782873344u,        0.9569336314u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5782873344u,        0.9569336314u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5631400122u,        0.9990242585u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5631400122u,        0.9990242585u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5932440289u,        0.9860870341u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5932440289u,        0.9860870341u) withcolor MidnightBlue;
draw (        0.5782873344u,        0.9569336314u)--(        0.5932554405u,        0.9860811767u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6257882113u,        0.9899853315u)--(        0.5932554405u,        0.9860811767u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5782873344u,        0.9569336314u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5782873344u,        0.9569336314u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6257882113u,        0.9899853315u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6257882113u,        0.9899853315u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5932554405u,        0.9860811767u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5932554405u,        0.9860811767u) withcolor MidnightBlue;
draw (        0.1651731031u,        1.1372176861u)--(        0.1350000000u,        1.1500000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.1651731031u,        1.1372176861u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1651731031u,        1.1372176861u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1350000000u,        1.1500000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1350000000u,        1.1500000000u) withcolor Yellow;
draw (        0.1651731031u,        1.1372176861u)--(        0.1973687942u,        1.1433054669u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2212106333u,        1.1208290280u)--(        0.1973687942u,        1.1433054669u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.1651731031u,        1.1372176861u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1651731031u,        1.1372176861u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2212106333u,        1.1208290280u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2212106333u,        1.1208290280u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1973687942u,        1.1433054669u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1973687942u,        1.1433054669u) withcolor MidnightBlue;
draw (        0.1651731031u,        1.1372176861u)--(        0.1973631697u,        1.1433351373u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2147377429u,        1.1711154866u)--(        0.1973631697u,        1.1433351373u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.1651731031u,        1.1372176861u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1651731031u,        1.1372176861u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2147377429u,        1.1711154866u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2147377429u,        1.1711154866u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1973631697u,        1.1433351373u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1973631697u,        1.1433351373u) withcolor MidnightBlue;
draw (        0.2212106333u,        1.1208290280u)--(        0.1973843177u,        1.1433219220u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2147377429u,        1.1711154866u)--(        0.1973843177u,        1.1433219220u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2212106333u,        1.1208290280u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2212106333u,        1.1208290280u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2147377429u,        1.1711154866u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2147377429u,        1.1711154866u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.1973843177u,        1.1433219220u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.1973843177u,        1.1433219220u) withcolor MidnightBlue;
draw (        0.2212106333u,        1.1208290280u)--(        0.2446042482u,        1.1437716168u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2699127465u,        1.1229604776u)--(        0.2446042482u,        1.1437716168u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2212106333u,        1.1208290280u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2212106333u,        1.1208290280u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2699127465u,        1.1229604776u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2699127465u,        1.1229604776u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2446042482u,        1.1437716168u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2446042482u,        1.1437716168u) withcolor MidnightBlue;
draw (        0.2212106333u,        1.1208290280u)--(        0.2445951625u,        1.1437808775u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2545346837u,        1.1750031374u)--(        0.2445951625u,        1.1437808775u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2212106333u,        1.1208290280u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2212106333u,        1.1208290280u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2545346837u,        1.1750031374u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2545346837u,        1.1750031374u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2445951625u,        1.1437808775u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2445951625u,        1.1437808775u) withcolor MidnightBlue;
draw (        0.3066548997u,        1.0808615973u)--(        0.3005337274u,        1.0486722380u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3313781856u,        1.0376158828u)--(        0.3005337274u,        1.0486722380u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3066548997u,        1.0808615973u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3066548997u,        1.0808615973u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3313781856u,        1.0376158828u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3313781856u,        1.0376158828u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3005337274u,        1.0486722380u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3005337274u,        1.0486722380u) withcolor Green;
draw (        0.3066548997u,        1.0808615973u)--(        0.3011812377u,        1.1131673649u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2699127465u,        1.1229604776u)--(        0.3011812377u,        1.1131673649u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3066548997u,        1.0808615973u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3066548997u,        1.0808615973u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2699127465u,        1.1229604776u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2699127465u,        1.1229604776u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3011812377u,        1.1131673649u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3011812377u,        1.1131673649u) withcolor MidnightBlue;
draw (        0.3066548997u,        1.0808615973u)--(        0.3011835531u,        1.1131677571u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3297577457u,        1.1292033149u)--(        0.3011835531u,        1.1131677571u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3066548997u,        1.0808615973u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3066548997u,        1.0808615973u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3297577457u,        1.1292033149u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3297577457u,        1.1292033149u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3011835531u,        1.1131677571u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3011835531u,        1.1131677571u) withcolor MidnightBlue;
draw (        0.3066548997u,        1.0808615973u)--(        0.3352290924u,        1.0968971550u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3297577457u,        1.1292033149u)--(        0.3352290924u,        1.0968971550u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3066548997u,        1.0808615973u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3066548997u,        1.0808615973u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3297577457u,        1.1292033149u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3297577457u,        1.1292033149u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3352290924u,        1.0968971550u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3352290924u,        1.0968971550u) withcolor MidnightBlue;
draw (        0.3066548997u,        1.0808615973u)--(        0.3352296154u,        1.0968962230u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3621204243u,        1.0781743644u)--(        0.3352296154u,        1.0968962230u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3066548997u,        1.0808615973u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3066548997u,        1.0808615973u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3621204243u,        1.0781743644u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3621204243u,        1.0781743644u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3352296154u,        1.0968962230u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3352296154u,        1.0968962230u) withcolor MidnightBlue;
draw (        0.3313781856u,        1.0376158828u)--(        0.3632011179u,        1.0454213027u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3873494378u,        1.0232744690u)--(        0.3632011179u,        1.0454213027u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3313781856u,        1.0376158828u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3313781856u,        1.0376158828u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3873494378u,        1.0232744690u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3873494378u,        1.0232744690u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3632011179u,        1.0454213027u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3632011179u,        1.0454213027u) withcolor MidnightBlue;
draw (        0.3313781856u,        1.0376158828u)--(        0.3631999759u,        1.0454259573u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3621204243u,        1.0781743644u)--(        0.3631999759u,        1.0454259573u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3313781856u,        1.0376158828u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3313781856u,        1.0376158828u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3621204243u,        1.0781743644u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3621204243u,        1.0781743644u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3631999759u,        1.0454259573u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3631999759u,        1.0454259573u) withcolor MidnightBlue;
draw (        0.3873494378u,        1.0232744690u)--(        0.4166752287u,        1.0378902682u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4449834019u,        1.0213896339u)--(        0.4166752287u,        1.0378902682u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3873494378u,        1.0232744690u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3873494378u,        1.0232744690u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4449834019u,        1.0213896339u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4449834019u,        1.0213896339u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4166752287u,        1.0378902682u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4166752287u,        1.0378902682u) withcolor MidnightBlue;
draw (        0.3873494378u,        1.0232744690u)--(        0.3632055567u,        1.0454261417u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3621204243u,        1.0781743644u)--(        0.3632055567u,        1.0454261417u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3873494378u,        1.0232744690u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3873494378u,        1.0232744690u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3621204243u,        1.0781743644u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3621204243u,        1.0781743644u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3632055567u,        1.0454261417u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3632055567u,        1.0454261417u) withcolor MidnightBlue;
draw (        0.3873494378u,        1.0232744690u)--(        0.3862643054u,        1.0560226917u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3621204243u,        1.0781743644u)--(        0.3862643054u,        1.0560226917u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3873494378u,        1.0232744690u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3873494378u,        1.0232744690u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3621204243u,        1.0781743644u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3621204243u,        1.0781743644u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3862643054u,        1.0560226917u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3862643054u,        1.0560226917u) withcolor MidnightBlue;
draw (        0.3873494378u,        1.0232744690u)--(        0.3862672432u,        1.0560227889u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4155877989u,        1.0706490875u)--(        0.3862672432u,        1.0560227889u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3873494378u,        1.0232744690u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3873494378u,        1.0232744690u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4155877989u,        1.0706490875u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4155877989u,        1.0706490875u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3862672432u,        1.0560227889u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3862672432u,        1.0560227889u) withcolor MidnightBlue;
draw (        0.3873494378u,        1.0232744690u)--(        0.4166699935u,        1.0379007676u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4155877989u,        1.0706490875u)--(        0.4166699935u,        1.0379007676u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3873494378u,        1.0232744690u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3873494378u,        1.0232744690u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4155877989u,        1.0706490875u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4155877989u,        1.0706490875u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4166699935u,        1.0379007676u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4166699935u,        1.0379007676u) withcolor MidnightBlue;
draw (        0.4449834019u,        1.0213896339u)--(        0.4769442970u,        1.0286093832u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5009532695u,        1.0063115618u)--(        0.4769442970u,        1.0286093832u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4449834019u,        1.0213896339u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4449834019u,        1.0213896339u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5009532695u,        1.0063115618u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5009532695u,        1.0063115618u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4769442970u,        1.0286093832u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4769442970u,        1.0286093832u) withcolor MidnightBlue;
draw (        0.4449834019u,        1.0213896339u)--(        0.4166815759u,        1.0379011524u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4155877989u,        1.0706490875u)--(        0.4166815759u,        1.0379011524u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4449834019u,        1.0213896339u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4449834019u,        1.0213896339u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4155877989u,        1.0706490875u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4155877989u,        1.0706490875u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4166815759u,        1.0379011524u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4166815759u,        1.0379011524u) withcolor MidnightBlue;
draw (        0.4449834019u,        1.0213896339u)--(        0.4769412007u,        1.0286230764u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4738750965u,        1.0612455013u)--(        0.4769412007u,        1.0286230764u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4449834019u,        1.0213896339u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4449834019u,        1.0213896339u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4738750965u,        1.0612455013u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4738750965u,        1.0612455013u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4769412007u,        1.0286230764u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4769412007u,        1.0286230764u) withcolor MidnightBlue;
draw (        0.5009532695u,        1.0063115618u)--(        0.5309209256u,        0.9930615359u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5631400122u,        0.9990242585u)--(        0.5309209256u,        0.9930615359u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5009532695u,        1.0063115618u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5009532695u,        1.0063115618u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5631400122u,        0.9990242585u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5631400122u,        0.9990242585u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5309209256u,        0.9930615359u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5309209256u,        0.9930615359u) withcolor MidnightBlue;
draw (        0.5009532695u,        1.0063115618u)--(        0.4769585414u,        1.0286247109u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4738750965u,        1.0612455013u)--(        0.4769585414u,        1.0286247109u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5009532695u,        1.0063115618u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5009532695u,        1.0063115618u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4738750965u,        1.0612455013u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4738750965u,        1.0612455013u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4769585414u,        1.0286247109u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4769585414u,        1.0286247109u) withcolor MidnightBlue;
draw (        0.5631400122u,        0.9990242585u)--(        0.5956758172u,        1.0029030451u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6257882113u,        0.9899853315u)--(        0.5956758172u,        1.0029030451u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5631400122u,        0.9990242585u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5631400122u,        0.9990242585u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6257882113u,        0.9899853315u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6257882113u,        0.9899853315u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5956758172u,        1.0029030451u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5956758172u,        1.0029030451u) withcolor MidnightBlue;
draw (        0.5631400122u,        0.9990242585u)--(        0.5932524062u,        0.9861065449u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6257882113u,        0.9899853315u)--(        0.5932524062u,        0.9861065449u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5631400122u,        0.9990242585u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5631400122u,        0.9990242585u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6257882113u,        0.9899853315u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6257882113u,        0.9899853315u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5932524062u,        0.9861065449u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5932524062u,        0.9861065449u) withcolor MidnightBlue;
draw (        0.5631400122u,        0.9990242585u)--(        0.5956764262u,        1.0028979339u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5946878703u,        1.0356492141u)--(        0.5956764262u,        1.0028979339u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5631400122u,        0.9990242585u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5631400122u,        0.9990242585u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5946878703u,        1.0356492141u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5946878703u,        1.0356492141u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5956764262u,        1.0028979339u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5956764262u,        1.0028979339u) withcolor MidnightBlue;
draw (        0.6257882113u,        0.9899853315u)--(        0.6550000000u,        0.9750000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6257882113u,        0.9899853315u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6257882113u,        0.9899853315u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6550000000u,        0.9750000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6550000000u,        0.9750000000u) withcolor Yellow;
draw (        0.6257882113u,        0.9899853315u)--(        0.5956735885u,        1.0028978484u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5946878703u,        1.0356492141u)--(        0.5956735885u,        1.0028978484u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6257882113u,        0.9899853315u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6257882113u,        0.9899853315u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5946878703u,        1.0356492141u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5946878703u,        1.0356492141u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5956735885u,        1.0028978484u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5956735885u,        1.0028978484u) withcolor MidnightBlue;
draw (        0.6257882113u,        0.9899853315u)--(        0.6248024931u,        1.0227366973u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5946878703u,        1.0356492141u)--(        0.6248024931u,        1.0227366973u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6257882113u,        0.9899853315u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6257882113u,        0.9899853315u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5946878703u,        1.0356492141u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5946878703u,        1.0356492141u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6248024931u,        1.0227366973u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6248024931u,        1.0227366973u) withcolor MidnightBlue;
draw (        0.6257882113u,        0.9899853315u)--(        0.6248337958u,        1.0227376244u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6573725793u,        1.0265913456u)--(        0.6248337958u,        1.0227376244u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6257882113u,        0.9899853315u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6257882113u,        0.9899853315u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6573725793u,        1.0265913456u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6573725793u,        1.0265913456u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6248337958u,        1.0227376244u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6248337958u,        1.0227376244u) withcolor MidnightBlue;
draw (        0.2699127465u,        1.1229604776u)--(        0.2984883464u,        1.1389935275u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3297577457u,        1.1292033149u)--(        0.2984883464u,        1.1389935275u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2699127465u,        1.1229604776u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2699127465u,        1.1229604776u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3297577457u,        1.1292033149u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3297577457u,        1.1292033149u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2984883464u,        1.1389935275u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2984883464u,        1.1389935275u) withcolor MidnightBlue;
draw (        0.2699127465u,        1.1229604776u)--(        0.3011821458u,        1.1131702650u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3297577457u,        1.1292033149u)--(        0.3011821458u,        1.1131702650u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2699127465u,        1.1229604776u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2699127465u,        1.1229604776u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3297577457u,        1.1292033149u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3297577457u,        1.1292033149u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3011821458u,        1.1131702650u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3011821458u,        1.1131702650u) withcolor MidnightBlue;
draw (        0.2699127465u,        1.1229604776u)--(        0.2446084011u,        1.1437766661u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2545346837u,        1.1750031374u)--(        0.2446084011u,        1.1437766661u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2699127465u,        1.1229604776u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2699127465u,        1.1229604776u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2545346837u,        1.1750031374u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2545346837u,        1.1750031374u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2446084011u,        1.1437766661u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2446084011u,        1.1437766661u) withcolor MidnightBlue;
draw (        0.2699127465u,        1.1229604776u)--(        0.2798390292u,        1.1541869490u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2545346837u,        1.1750031374u)--(        0.2798390292u,        1.1541869490u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2699127465u,        1.1229604776u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2699127465u,        1.1229604776u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2545346837u,        1.1750031374u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2545346837u,        1.1750031374u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2798390292u,        1.1541869490u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2798390292u,        1.1541869490u) withcolor MidnightBlue;
draw (        0.2699127465u,        1.1229604776u)--(        0.2798427553u,        1.1541857643u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3084077687u,        1.1702376677u)--(        0.2798427553u,        1.1541857643u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2699127465u,        1.1229604776u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2699127465u,        1.1229604776u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3084077687u,        1.1702376677u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3084077687u,        1.1702376677u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2798427553u,        1.1541857643u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2798427553u,        1.1541857643u) withcolor MidnightBlue;
draw (        0.2699127465u,        1.1229604776u)--(        0.2984777600u,        1.1390123811u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3084077687u,        1.1702376677u)--(        0.2984777600u,        1.1390123811u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2699127465u,        1.1229604776u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2699127465u,        1.1229604776u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3084077687u,        1.1702376677u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3084077687u,        1.1702376677u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2984777600u,        1.1390123811u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2984777600u,        1.1390123811u) withcolor MidnightBlue;
draw (        0.3297577457u,        1.1292033149u)--(        0.3566477491u,        1.1104802992u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3621204243u,        1.0781743644u)--(        0.3566477491u,        1.1104802992u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3297577457u,        1.1292033149u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3297577457u,        1.1292033149u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3621204243u,        1.0781743644u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3621204243u,        1.0781743644u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3566477491u,        1.1104802992u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3566477491u,        1.1104802992u) withcolor MidnightBlue;
draw (        0.3297577457u,        1.1292033149u)--(        0.3352304210u,        1.0968973801u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3621204243u,        1.0781743644u)--(        0.3352304210u,        1.0968973801u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3297577457u,        1.1292033149u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3297577457u,        1.1292033149u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3621204243u,        1.0781743644u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3621204243u,        1.0781743644u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3352304210u,        1.0968973801u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3352304210u,        1.0968973801u) withcolor MidnightBlue;
draw (        0.3297577457u,        1.1292033149u)--(        0.2984927605u,        1.1390076147u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3084077687u,        1.1702376677u)--(        0.2984927605u,        1.1390076147u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3297577457u,        1.1292033149u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3297577457u,        1.1292033149u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3084077687u,        1.1702376677u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3084077687u,        1.1702376677u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2984927605u,        1.1390076147u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2984927605u,        1.1390076147u) withcolor MidnightBlue;
draw (        0.3297577457u,        1.1292033149u)--(        0.3608689442u,        1.1394851887u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3657710217u,        1.1718826150u)--(        0.3608689442u,        1.1394851887u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3297577457u,        1.1292033149u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3297577457u,        1.1292033149u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3657710217u,        1.1718826150u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3657710217u,        1.1718826150u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3608689442u,        1.1394851887u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3608689442u,        1.1394851887u) withcolor MidnightBlue;
draw (        0.3297577457u,        1.1292033149u)--(        0.3608720714u,        1.1394757215u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3877668800u,        1.1207596090u)--(        0.3608720714u,        1.1394757215u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3297577457u,        1.1292033149u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3297577457u,        1.1292033149u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3877668800u,        1.1207596090u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3877668800u,        1.1207596090u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3608720714u,        1.1394757215u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3608720714u,        1.1394757215u) withcolor MidnightBlue;
draw (        0.3297577457u,        1.1292033149u)--(        0.3566525544u,        1.1104872024u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3877668800u,        1.1207596090u)--(        0.3566525544u,        1.1104872024u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3297577457u,        1.1292033149u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3297577457u,        1.1292033149u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3877668800u,        1.1207596090u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3877668800u,        1.1207596090u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3566525544u,        1.1104872024u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3566525544u,        1.1104872024u) withcolor MidnightBlue;
draw (        0.3621204243u,        1.0781743644u)--(        0.3862662399u,        1.0560248003u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4155877989u,        1.0706490875u)--(        0.3862662399u,        1.0560248003u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3621204243u,        1.0781743644u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3621204243u,        1.0781743644u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4155877989u,        1.0706490875u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4155877989u,        1.0706490875u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3862662399u,        1.0560248003u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3862662399u,        1.0560248003u) withcolor MidnightBlue;
draw (        0.3621204243u,        1.0781743644u)--(        0.3566544590u,        1.1104814351u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3877668800u,        1.1207596090u)--(        0.3566544590u,        1.1104814351u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3621204243u,        1.0781743644u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3621204243u,        1.0781743644u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3877668800u,        1.1207596090u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3877668800u,        1.1207596090u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3566544590u,        1.1104814351u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3566544590u,        1.1104814351u) withcolor MidnightBlue;
draw (        0.4155877989u,        1.0706490875u)--(        0.4469960525u,        1.0799842470u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4738750965u,        1.0612455013u)--(        0.4469960525u,        1.0799842470u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4155877989u,        1.0706490875u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4155877989u,        1.0706490875u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4738750965u,        1.0612455013u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4738750965u,        1.0612455013u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4469960525u,        1.0799842470u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4469960525u,        1.0799842470u) withcolor MidnightBlue;
draw (        0.4155877989u,        1.0706490875u)--(        0.4155661195u,        1.1034152763u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3877668800u,        1.1207596090u)--(        0.4155661195u,        1.1034152763u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4155877989u,        1.0706490875u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4155877989u,        1.0706490875u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3877668800u,        1.1207596090u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3877668800u,        1.1207596090u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4155661195u,        1.1034152763u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4155661195u,        1.1034152763u) withcolor MidnightBlue;
draw (        0.4155877989u,        1.0706490875u)--(        0.4155782200u,        1.1034152820u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4469859884u,        1.1127520740u)--(        0.4155782200u,        1.1034152820u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4155877989u,        1.0706490875u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4155877989u,        1.0706490875u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4469859884u,        1.1127520740u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4469859884u,        1.1127520740u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4155782200u,        1.1034152820u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4155782200u,        1.1034152820u) withcolor MidnightBlue;
draw (        0.4155877989u,        1.0706490875u)--(        0.4469955673u,        1.0799858794u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4469859884u,        1.1127520740u)--(        0.4469955673u,        1.0799858794u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4155877989u,        1.0706490875u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4155877989u,        1.0706490875u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4469859884u,        1.1127520740u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4469859884u,        1.1127520740u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4469955673u,        1.0799858794u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4469955673u,        1.0799858794u) withcolor MidnightBlue;
draw (        0.4738750965u,        1.0612455013u)--(        0.4469971910u,        1.0799858799u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4469859884u,        1.1127520740u)--(        0.4469971910u,        1.0799858799u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4738750965u,        1.0612455013u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4738750965u,        1.0612455013u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4469859884u,        1.1127520740u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4469859884u,        1.1127520740u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4469971910u,        1.0799858799u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4469971910u,        1.0799858799u) withcolor MidnightBlue;
draw (        0.4738750965u,        1.0612455013u)--(        0.4738638939u,        1.0940116953u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4469859884u,        1.1127520740u)--(        0.4738638939u,        1.0940116953u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4738750965u,        1.0612455013u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4738750965u,        1.0612455013u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4469859884u,        1.1127520740u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4469859884u,        1.1127520740u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4738638939u,        1.0940116953u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4738638939u,        1.0940116953u) withcolor MidnightBlue;
draw (        0.4738750965u,        1.0612455013u)--(        0.4738800962u,        1.0940116968u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5001059769u,        1.1136541706u)--(        0.4738800962u,        1.0940116968u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4738750965u,        1.0612455013u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4738750965u,        1.0612455013u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5001059769u,        1.1136541706u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5001059769u,        1.1136541706u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4738800962u,        1.0940116968u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4738800962u,        1.0940116968u) withcolor MidnightBlue;
draw (        0.5946878703u,        1.0356492141u)--(        0.6272333175u,        1.0394462474u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6573725793u,        1.0265913456u)--(        0.6272333175u,        1.0394462474u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5946878703u,        1.0356492141u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5946878703u,        1.0356492141u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6573725793u,        1.0265913456u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6573725793u,        1.0265913456u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6272333175u,        1.0394462474u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6272333175u,        1.0394462474u) withcolor MidnightBlue;
draw (        0.5946878703u,        1.0356492141u)--(        0.6248271320u,        1.0227943123u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6573725793u,        1.0265913456u)--(        0.6248271320u,        1.0227943123u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5946878703u,        1.0356492141u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5946878703u,        1.0356492141u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6573725793u,        1.0265913456u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6573725793u,        1.0265913456u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6248271320u,        1.0227943123u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6248271320u,        1.0227943123u) withcolor MidnightBlue;
draw (        0.5946878703u,        1.0356492141u)--(        0.6272332444u,        1.0394468746u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6272411329u,        1.0722130697u)--(        0.6272332444u,        1.0394468746u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5946878703u,        1.0356492141u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5946878703u,        1.0356492141u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6272411329u,        1.0722130697u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6272411329u,        1.0722130697u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6272332444u,        1.0394468746u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6272332444u,        1.0394468746u) withcolor MidnightBlue;
draw (        0.6573725793u,        1.0265913456u)--(        0.6900000000u,        1.0300000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6573725793u,        1.0265913456u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6573725793u,        1.0265913456u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6900000000u,        1.0300000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6900000000u,        1.0300000000u) withcolor Yellow;
draw (        0.6573725793u,        1.0265913456u)--(        0.6272335850u,        1.0394468746u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6272411329u,        1.0722130697u)--(        0.6272335850u,        1.0394468746u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6573725793u,        1.0265913456u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6573725793u,        1.0265913456u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6272411329u,        1.0722130697u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6272411329u,        1.0722130697u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6272335850u,        1.0394468746u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6272335850u,        1.0394468746u) withcolor MidnightBlue;
draw (        0.6573725793u,        1.0265913456u)--(        0.6573801271u,        1.0593575407u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6272411329u,        1.0722130697u)--(        0.6573801271u,        1.0593575407u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6573725793u,        1.0265913456u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6573725793u,        1.0265913456u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6272411329u,        1.0722130697u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6272411329u,        1.0722130697u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6573801271u,        1.0593575407u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6573801271u,        1.0593575407u) withcolor Green;
draw (        0.2545346837u,        1.1750031374u)--(        0.2798408475u,        1.1541891596u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3084077687u,        1.1702376677u)--(        0.2798408475u,        1.1541891596u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.2545346837u,        1.1750031374u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2545346837u,        1.1750031374u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3084077687u,        1.1702376677u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3084077687u,        1.1702376677u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.2798408475u,        1.1541891596u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.2798408475u,        1.1541891596u) withcolor MidnightBlue;
draw (        0.3657710217u,        1.1718826150u)--(        0.3926601653u,        1.1531583645u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3877668800u,        1.1207596090u)--(        0.3926601653u,        1.1531583645u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3657710217u,        1.1718826150u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3657710217u,        1.1718826150u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3877668800u,        1.1207596090u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3877668800u,        1.1207596090u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3926601653u,        1.1531583645u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3926601653u,        1.1531583645u) withcolor MidnightBlue;
draw (        0.3657710217u,        1.1718826150u)--(        0.3608777365u,        1.1394838595u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3877668800u,        1.1207596090u)--(        0.3608777365u,        1.1394838595u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3657710217u,        1.1718826150u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3657710217u,        1.1718826150u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3877668800u,        1.1207596090u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3877668800u,        1.1207596090u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3608777365u,        1.1394838595u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3608777365u,        1.1394838595u) withcolor MidnightBlue;
draw (        0.3657710217u,        1.1718826150u)--(        0.3926631184u,        1.1531626060u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4240762185u,        1.1624814437u)--(        0.3926631184u,        1.1531626060u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3657710217u,        1.1718826150u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3657710217u,        1.1718826150u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4240762185u,        1.1624814437u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4240762185u,        1.1624814437u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3926631184u,        1.1531626060u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3926631184u,        1.1531626060u) withcolor MidnightBlue;
draw (        0.3877668800u,        1.1207596090u)--(        0.4191785515u,        1.1300832613u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4469859884u,        1.1127520740u)--(        0.4191785515u,        1.1300832613u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3877668800u,        1.1207596090u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3877668800u,        1.1207596090u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4469859884u,        1.1127520740u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4469859884u,        1.1127520740u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4191785515u,        1.1300832613u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4191785515u,        1.1300832613u) withcolor MidnightBlue;
draw (        0.3877668800u,        1.1207596090u)--(        0.4155743169u,        1.1034284217u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4469859884u,        1.1127520740u)--(        0.4155743169u,        1.1034284217u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3877668800u,        1.1207596090u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3877668800u,        1.1207596090u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4469859884u,        1.1127520740u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4469859884u,        1.1127520740u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4155743169u,        1.1034284217u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4155743169u,        1.1034284217u) withcolor MidnightBlue;
draw (        0.3877668800u,        1.1207596090u)--(        0.3926645747u,        1.1531576982u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4240762185u,        1.1624814437u)--(        0.3926645747u,        1.1531576982u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3877668800u,        1.1207596090u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3877668800u,        1.1207596090u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4240762185u,        1.1624814437u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4240762185u,        1.1624814437u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.3926645747u,        1.1531576982u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3926645747u,        1.1531576982u) withcolor MidnightBlue;
draw (        0.3877668800u,        1.1207596090u)--(        0.4191785238u,        1.1300833545u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4240762185u,        1.1624814437u)--(        0.4191785238u,        1.1300833545u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.3877668800u,        1.1207596090u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.3877668800u,        1.1207596090u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4240762185u,        1.1624814437u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4240762185u,        1.1624814437u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4191785238u,        1.1300833545u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4191785238u,        1.1300833545u) withcolor MidnightBlue;
draw (        0.4469859884u,        1.1127520740u)--(        0.4738717058u,        1.0940229043u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5001059769u,        1.1136541706u)--(        0.4738717058u,        1.0940229043u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4469859884u,        1.1127520740u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4469859884u,        1.1127520740u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5001059769u,        1.1136541706u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5001059769u,        1.1136541706u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4738717058u,        1.0940229043u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4738717058u,        1.0940229043u) withcolor MidnightBlue;
draw (        0.4469859884u,        1.1127520740u)--(        0.4191786022u,        1.1300833427u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4240762185u,        1.1624814437u)--(        0.4191786022u,        1.1300833427u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4469859884u,        1.1127520740u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4469859884u,        1.1127520740u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4240762185u,        1.1624814437u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4240762185u,        1.1624814437u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4191786022u,        1.1300833427u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4191786022u,        1.1300833427u) withcolor MidnightBlue;
draw (        0.4469859884u,        1.1127520740u)--(        0.4518836047u,        1.1451501750u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4240762185u,        1.1624814437u)--(        0.4518836047u,        1.1451501750u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4469859884u,        1.1127520740u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4469859884u,        1.1127520740u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4240762185u,        1.1624814437u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4240762185u,        1.1624814437u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4518836047u,        1.1451501750u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4518836047u,        1.1451501750u) withcolor MidnightBlue;
draw (        0.4469859884u,        1.1127520740u)--(        0.4518968927u,        1.1451481635u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4824832725u,        1.1568996303u)--(        0.4518968927u,        1.1451481635u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4469859884u,        1.1127520740u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4469859884u,        1.1127520740u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4824832725u,        1.1568996303u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4824832725u,        1.1568996303u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4518968927u,        1.1451481635u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4518968927u,        1.1451481635u) withcolor MidnightBlue;
draw (        0.5001059769u,        1.1136541706u)--(        0.5125827768u,        1.1439519109u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4824832725u,        1.1568996303u)--(        0.5125827768u,        1.1439519109u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5001059769u,        1.1136541706u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5001059769u,        1.1136541706u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4824832725u,        1.1568996303u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4824832725u,        1.1568996303u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5125827768u,        1.1439519109u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5125827768u,        1.1439519109u) withcolor MidnightBlue;
draw (        0.5001059769u,        1.1136541706u)--(        0.5125872908u,        1.1439500516u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5428526001u,        1.1565053156u)--(        0.5125872908u,        1.1439500516u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5001059769u,        1.1136541706u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5001059769u,        1.1136541706u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5428526001u,        1.1565053156u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5428526001u,        1.1565053156u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5125872908u,        1.1439500516u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5125872908u,        1.1439500516u) withcolor MidnightBlue;
draw (        0.5809355355u,        1.0885043402u)--(        0.5900000000u,        1.1200000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5809355355u,        1.0885043402u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5809355355u,        1.0885043402u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5900000000u,        1.1200000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5900000000u,        1.1200000000u) withcolor Yellow;
draw (        0.5809355355u,        1.0885043402u)--(        0.6112926568u,        1.1008359575u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6272411329u,        1.0722130697u)--(        0.6112926568u,        1.1008359575u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5809355355u,        1.0885043402u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5809355355u,        1.0885043402u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6272411329u,        1.0722130697u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6272411329u,        1.0722130697u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6112926568u,        1.1008359575u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6112926568u,        1.1008359575u) withcolor Green;
draw (        0.4240762185u,        1.1624814437u)--(        0.4546681430u,        1.1742184688u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4824832725u,        1.1568996303u)--(        0.4546681430u,        1.1742184688u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4240762185u,        1.1624814437u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4240762185u,        1.1624814437u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4824832725u,        1.1568996303u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4824832725u,        1.1568996303u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4546681430u,        1.1742184688u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4546681430u,        1.1742184688u) withcolor Green;
draw (        0.4240762185u,        1.1624814437u)--(        0.4518913481u,        1.1451626052u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4824832725u,        1.1568996303u)--(        0.4518913481u,        1.1451626052u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4240762185u,        1.1624814437u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4240762185u,        1.1624814437u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4824832725u,        1.1568996303u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4824832725u,        1.1568996303u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.4518913481u,        1.1451626052u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4518913481u,        1.1451626052u) withcolor MidnightBlue;
draw (        0.4824832725u,        1.1568996303u)--(        0.5125846823u,        1.1439563415u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5428526001u,        1.1565053156u)--(        0.5125846823u,        1.1439563415u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.4824832725u,        1.1568996303u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.4824832725u,        1.1568996303u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5428526001u,        1.1565053156u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5428526001u,        1.1565053156u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5125846823u,        1.1439563415u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5125846823u,        1.1439563415u) withcolor MidnightBlue;
draw (        0.5428526001u,        1.1565053156u)--(        0.5750000000u,        1.1500000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5428526001u,        1.1565053156u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5428526001u,        1.1565053156u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5750000000u,        1.1500000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5750000000u,        1.1500000000u) withcolor Yellow;
% r =         0.0327661960
endfig;
end;
