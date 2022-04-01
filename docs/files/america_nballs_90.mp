prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(90);
u := 20cm;
path pols[];
pair p[];
R := .0309688348u;
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

x1 = .6013406827u;
y1 = .2849069132u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .5840484865u;
y2 = .3290985707u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .6163119377u;
y3 = .3569267822u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .5870051250u;
y4 = .3836550008u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .6328616584u;
y5 = .3997731280u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .5865763239u;
y6 = .4356951430u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .6322410427u;
y7 = .4518193379u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .6729846640u;
y8 = .4459031899u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .5904128523u;
y9 = .4851277115u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = .6459781071u;
y10 = .4990946821u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;

x11 = .6934321912u;
y11 = .4922578756u;
draw fullcircle scaled D shifted z11 withcolor RoyalBlue withpen pencircle scaled 1.0;

x12 = .6982155107u;
y12 = .5872194905u;
draw fullcircle scaled D shifted z12 withcolor RoyalBlue withpen pencircle scaled 1.0;

x13 = .5802077452u;
y13 = .5367241529u;
draw fullcircle scaled D shifted z13 withcolor RoyalBlue withpen pencircle scaled 1.0;

x14 = .6262631889u;
y14 = .5324953778u;
draw fullcircle scaled D shifted z14 withcolor RoyalBlue withpen pencircle scaled 1.0;

x15 = .6828065393u;
y15 = .5442861560u;
draw fullcircle scaled D shifted z15 withcolor RoyalBlue withpen pencircle scaled 1.0;

x16 = .7346495018u;
y16 = .5195828424u;
draw fullcircle scaled D shifted z16 withcolor RoyalBlue withpen pencircle scaled 1.0;

x17 = .5432420020u;
y17 = .5774231174u;
draw fullcircle scaled D shifted z17 withcolor RoyalBlue withpen pencircle scaled 1.0;

x18 = .5971919913u;
y18 = .5850596971u;
draw fullcircle scaled D shifted z18 withcolor RoyalBlue withpen pencircle scaled 1.0;

x19 = .6443253844u;
y19 = .5792739342u;
draw fullcircle scaled D shifted z19 withcolor RoyalBlue withpen pencircle scaled 1.0;

x20 = .7410417328u;
y20 = .5555878687u;
draw fullcircle scaled D shifted z20 withcolor RoyalBlue withpen pencircle scaled 1.0;

x21 = .7557860761u;
y21 = .5997246016u;
draw fullcircle scaled D shifted z21 withcolor RoyalBlue withpen pencircle scaled 1.0;

x22 = .5554474438u;
y22 = .6193418400u;
draw fullcircle scaled D shifted z22 withcolor RoyalBlue withpen pencircle scaled 1.0;

x23 = .6122210394u;
y23 = .6302869245u;
draw fullcircle scaled D shifted z23 withcolor RoyalBlue withpen pencircle scaled 1.0;

x24 = .6591014928u;
y24 = .6271318287u;
draw fullcircle scaled D shifted z24 withcolor RoyalBlue withpen pencircle scaled 1.0;

x25 = .7183405430u;
y25 = .6235044217u;
draw fullcircle scaled D shifted z25 withcolor RoyalBlue withpen pencircle scaled 1.0;

x26 = .5389761541u;
y26 = .6573459103u;
draw fullcircle scaled D shifted z26 withcolor RoyalBlue withpen pencircle scaled 1.0;

x27 = .5826625305u;
y27 = .6666451873u;
draw fullcircle scaled D shifted z27 withcolor RoyalBlue withpen pencircle scaled 1.0;

x28 = .6379881712u;
y28 = .6797469625u;
draw fullcircle scaled D shifted z28 withcolor RoyalBlue withpen pencircle scaled 1.0;

x29 = .6847608802u;
y29 = .6647994214u;
draw fullcircle scaled D shifted z29 withcolor RoyalBlue withpen pencircle scaled 1.0;

x30 = .4983253853u;
y30 = .7267565974u;
draw fullcircle scaled D shifted z30 withcolor RoyalBlue withpen pencircle scaled 1.0;

x31 = .5555210667u;
y31 = .7139382412u;
draw fullcircle scaled D shifted z31 withcolor RoyalBlue withpen pencircle scaled 1.0;

x32 = .6085728756u;
y32 = .7132003932u;
draw fullcircle scaled D shifted z32 withcolor RoyalBlue withpen pencircle scaled 1.0;

x33 = .4015589710u;
y33 = .8049487301u;
draw fullcircle scaled D shifted z33 withcolor RoyalBlue withpen pencircle scaled 1.0;

x34 = .4383219753u;
y34 = .7653461136u;
draw fullcircle scaled D shifted z34 withcolor RoyalBlue withpen pencircle scaled 1.0;

x35 = .4837562134u;
y35 = .7713895818u;
draw fullcircle scaled D shifted z35 withcolor RoyalBlue withpen pencircle scaled 1.0;

x36 = .3601114197u;
y36 = .8452050770u;
draw fullcircle scaled D shifted z36 withcolor RoyalBlue withpen pencircle scaled 1.0;

x37 = .4165513467u;
y37 = .8455446109u;
draw fullcircle scaled D shifted z37 withcolor RoyalBlue withpen pencircle scaled 1.0;

x38 = .5412341004u;
y38 = .8192683697u;
draw fullcircle scaled D shifted z38 withcolor RoyalBlue withpen pencircle scaled 1.0;

x39 = .3363451551u;
y39 = .8931081360u;
draw fullcircle scaled D shifted z39 withcolor RoyalBlue withpen pencircle scaled 1.0;

x40 = .3899487025u;
y40 = .8886963474u;
draw fullcircle scaled D shifted z40 withcolor RoyalBlue withpen pencircle scaled 1.0;

x41 = .4428745823u;
y41 = .8787455516u;
draw fullcircle scaled D shifted z41 withcolor RoyalBlue withpen pencircle scaled 1.0;

x42 = .4879141441u;
y42 = .8828277204u;
draw fullcircle scaled D shifted z42 withcolor RoyalBlue withpen pencircle scaled 1.0;

x43 = .5364456294u;
y43 = .8711574551u;
draw fullcircle scaled D shifted z43 withcolor RoyalBlue withpen pencircle scaled 1.0;

x44 = .3113577906u;
y44 = .9400246340u;
draw fullcircle scaled D shifted z44 withcolor RoyalBlue withpen pencircle scaled 1.0;

x45 = .3651920244u;
y45 = .9358664236u;
draw fullcircle scaled D shifted z45 withcolor RoyalBlue withpen pencircle scaled 1.0;

x46 = .4190345109u;
y46 = .9289572292u;
draw fullcircle scaled D shifted z46 withcolor RoyalBlue withpen pencircle scaled 1.0;

x47 = .4666860449u;
y47 = .9286039293u;
draw fullcircle scaled D shifted z47 withcolor RoyalBlue withpen pencircle scaled 1.0;

x48 = .5219432346u;
y48 = .9250041483u;
draw fullcircle scaled D shifted z48 withcolor RoyalBlue withpen pencircle scaled 1.0;

x49 = .5678032676u;
y49 = .9184100577u;
draw fullcircle scaled D shifted z49 withcolor RoyalBlue withpen pencircle scaled 1.0;

x50 = .3330409729u;
y50 = .9804481474u;
draw fullcircle scaled D shifted z50 withcolor RoyalBlue withpen pencircle scaled 1.0;

x51 = .3891017973u;
y51 = .9773058157u;
draw fullcircle scaled D shifted z51 withcolor RoyalBlue withpen pencircle scaled 1.0;

x52 = .4415580069u;
y52 = .9733885092u;
draw fullcircle scaled D shifted z52 withcolor RoyalBlue withpen pencircle scaled 1.0;

x53 = .4973862559u;
y53 = .9711985538u;
draw fullcircle scaled D shifted z53 withcolor RoyalBlue withpen pencircle scaled 1.0;

x54 = .5500493078u;
y54 = .9720730823u;
draw fullcircle scaled D shifted z54 withcolor RoyalBlue withpen pencircle scaled 1.0;

x55 = .5926306482u;
y55 = .9706954877u;
draw fullcircle scaled D shifted z55 withcolor RoyalBlue withpen pencircle scaled 1.0;

x56 = .3234370266u;
y56 = 1.0370727569u;
draw fullcircle scaled D shifted z56 withcolor RoyalBlue withpen pencircle scaled 1.0;

x57 = .3681028268u;
y57 = 1.0214525529u;
draw fullcircle scaled D shifted z57 withcolor RoyalBlue withpen pencircle scaled 1.0;

x58 = .4216063159u;
y58 = 1.0215797604u;
draw fullcircle scaled D shifted z58 withcolor RoyalBlue withpen pencircle scaled 1.0;

x59 = .4750851673u;
y59 = 1.0155248829u;
draw fullcircle scaled D shifted z59 withcolor RoyalBlue withpen pencircle scaled 1.0;

x60 = .5296621648u;
y60 = 1.0158876503u;
draw fullcircle scaled D shifted z60 withcolor RoyalBlue withpen pencircle scaled 1.0;

x61 = .5831951704u;
y61 = 1.0184505906u;
draw fullcircle scaled D shifted z61 withcolor RoyalBlue withpen pencircle scaled 1.0;

x62 = .6375703128u;
y62 = 1.0007339529u;
draw fullcircle scaled D shifted z62 withcolor RoyalBlue withpen pencircle scaled 1.0;

x63 = .6607047374u;
y63 = 1.0356856312u;
draw fullcircle scaled D shifted z63 withcolor RoyalBlue withpen pencircle scaled 1.0;

x64 = .3494547024u;
y64 = 1.0784899507u;
draw fullcircle scaled D shifted z64 withcolor RoyalBlue withpen pencircle scaled 1.0;

x65 = .3919160679u;
y65 = 1.0672600242u;
draw fullcircle scaled D shifted z65 withcolor RoyalBlue withpen pencircle scaled 1.0;

x66 = .4474515501u;
y66 = 1.0641061602u;
draw fullcircle scaled D shifted z66 withcolor RoyalBlue withpen pencircle scaled 1.0;

x67 = .4989344603u;
y67 = 1.0605534445u;
draw fullcircle scaled D shifted z67 withcolor RoyalBlue withpen pencircle scaled 1.0;

x68 = .5967242763u;
y68 = 1.0465465577u;
draw fullcircle scaled D shifted z68 withcolor RoyalBlue withpen pencircle scaled 1.0;

x69 = .6228165214u;
y69 = 1.0573376798u;
draw fullcircle scaled D shifted z69 withcolor RoyalBlue withpen pencircle scaled 1.0;

x70 = .1612995390u;
y70 = 1.1035931367u;
draw fullcircle scaled D shifted z70 withcolor RoyalBlue withpen pencircle scaled 1.0;

x71 = .2202729296u;
y71 = 1.1046765931u;
draw fullcircle scaled D shifted z71 withcolor RoyalBlue withpen pencircle scaled 1.0;

x72 = .2990518797u;
y72 = 1.0848394502u;
draw fullcircle scaled D shifted z72 withcolor RoyalBlue withpen pencircle scaled 1.0;

x73 = .3667466584u;
y73 = 1.1211944291u;
draw fullcircle scaled D shifted z73 withcolor RoyalBlue withpen pencircle scaled 1.0;

x74 = .4146764665u;
y74 = 1.1093597134u;
draw fullcircle scaled D shifted z74 withcolor RoyalBlue withpen pencircle scaled 1.0;

x75 = .4683118439u;
y75 = 1.1091635284u;
draw fullcircle scaled D shifted z75 withcolor RoyalBlue withpen pencircle scaled 1.0;

x76 = .5149952438u;
y76 = 1.1122807665u;
draw fullcircle scaled D shifted z76 withcolor RoyalBlue withpen pencircle scaled 1.0;

x77 = .5943425459u;
y77 = 1.0893154275u;
draw fullcircle scaled D shifted z77 withcolor RoyalBlue withpen pencircle scaled 1.0;

x78 = .1564686220u;
y78 = 1.1454744659u;
draw fullcircle scaled D shifted z78 withcolor RoyalBlue withpen pencircle scaled 1.0;

x79 = .2104612454u;
y79 = 1.1362983588u;
draw fullcircle scaled D shifted z79 withcolor RoyalBlue withpen pencircle scaled 1.0;

x80 = .2618244563u;
y80 = 1.1287032226u;
draw fullcircle scaled D shifted z80 withcolor RoyalBlue withpen pencircle scaled 1.0;

x81 = .3161426569u;
y81 = 1.1280149432u;
draw fullcircle scaled D shifted z81 withcolor RoyalBlue withpen pencircle scaled 1.0;

x82 = .3971005061u;
y82 = 1.1637581423u;
draw fullcircle scaled D shifted z82 withcolor RoyalBlue withpen pencircle scaled 1.0;

x83 = .4437762635u;
y83 = 1.1545040782u;
draw fullcircle scaled D shifted z83 withcolor RoyalBlue withpen pencircle scaled 1.0;

x84 = .4954895093u;
y84 = 1.1580217485u;
draw fullcircle scaled D shifted z84 withcolor RoyalBlue withpen pencircle scaled 1.0;

x85 = .5463600960u;
y85 = 1.1610126536u;
draw fullcircle scaled D shifted z85 withcolor RoyalBlue withpen pencircle scaled 1.0;

x86 = .1902373013u;
y86 = 1.1851026060u;
draw fullcircle scaled D shifted z86 withcolor RoyalBlue withpen pencircle scaled 1.0;

x87 = .2417333908u;
y87 = 1.1791821618u;
draw fullcircle scaled D shifted z87 withcolor RoyalBlue withpen pencircle scaled 1.0;

x88 = .2924457665u;
y88 = 1.1732225853u;
draw fullcircle scaled D shifted z88 withcolor RoyalBlue withpen pencircle scaled 1.0;

x89 = .3456387000u;
y89 = 1.1723255025u;
draw fullcircle scaled D shifted z89 withcolor RoyalBlue withpen pencircle scaled 1.0;

x90 = .5372660140u;
y90 = 1.2020187196u;
draw fullcircle scaled D shifted z90 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.6013406827u,        0.2849069132u)--(        0.6250000000u,        0.2650000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6013406827u,        0.2849069132u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6013406827u,        0.2849069132u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6250000000u,        0.2650000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6250000000u,        0.2650000000u) withcolor Yellow;
draw (        0.6013406827u,        0.2849069132u)--(        0.6112281698u,        0.3142549341u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5840484865u,        0.3290985707u)--(        0.6112281698u,        0.3142549341u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6013406827u,        0.2849069132u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6013406827u,        0.2849069132u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5840484865u,        0.3290985707u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5840484865u,        0.3290985707u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6112281698u,        0.3142549341u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6112281698u,        0.3142549341u) withcolor Green;
draw (        0.5870051250u,        0.3836550008u)--(        0.6035685650u,        0.4098221784u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6328616584u,        0.3997731280u)--(        0.6035685650u,        0.4098221784u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5870051250u,        0.3836550008u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5870051250u,        0.3836550008u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6328616584u,        0.3997731280u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6328616584u,        0.3997731280u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6035685650u,        0.4098221784u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6035685650u,        0.4098221784u) withcolor MidnightBlue;
draw (        0.5870051250u,        0.3836550008u)--(        0.5699990742u,        0.4095367118u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5865763239u,        0.4356951430u)--(        0.5699990742u,        0.4095367118u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5870051250u,        0.3836550008u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5870051250u,        0.3836550008u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5865763239u,        0.4356951430u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5865763239u,        0.4356951430u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5699990742u,        0.4095367118u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5699990742u,        0.4095367118u) withcolor Green;
draw (        0.5870051250u,        0.3836550008u)--(        0.6035823748u,        0.4098134319u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5865763239u,        0.4356951430u)--(        0.6035823748u,        0.4098134319u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5870051250u,        0.3836550008u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5870051250u,        0.3836550008u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5865763239u,        0.4356951430u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5865763239u,        0.4356951430u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6035823748u,        0.4098134319u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6035823748u,        0.4098134319u) withcolor MidnightBlue;
draw (        0.6328616584u,        0.3997731280u)--(        0.6600000000u,        0.3850000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6328616584u,        0.3997731280u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6328616584u,        0.3997731280u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6600000000u,        0.3850000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6600000000u,        0.3850000000u) withcolor Yellow;
draw (        0.6328616584u,        0.3997731280u)--(        0.6035606915u,        0.4097991975u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5865763239u,        0.4356951430u)--(        0.6035606915u,        0.4097991975u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6328616584u,        0.3997731280u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6328616584u,        0.3997731280u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5865763239u,        0.4356951430u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5865763239u,        0.4356951430u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6035606915u,        0.4097991975u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6035606915u,        0.4097991975u) withcolor MidnightBlue;
draw (        0.6328616584u,        0.3997731280u)--(        0.6157665244u,        0.4255960853u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6322410427u,        0.4518193379u)--(        0.6157665244u,        0.4255960853u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6328616584u,        0.3997731280u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6328616584u,        0.3997731280u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6322410427u,        0.4518193379u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6322410427u,        0.4518193379u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6157665244u,        0.4255960853u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6157665244u,        0.4255960853u) withcolor MidnightBlue;
draw (        0.6328616584u,        0.3997731280u)--(        0.6493361767u,        0.4259963806u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6322410427u,        0.4518193379u)--(        0.6493361767u,        0.4259963806u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6328616584u,        0.3997731280u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6328616584u,        0.3997731280u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6322410427u,        0.4518193379u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6322410427u,        0.4518193379u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6493361767u,        0.4259963806u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6493361767u,        0.4259963806u) withcolor MidnightBlue;
draw (        0.6328616584u,        0.3997731280u)--(        0.6566661832u,        0.4195825532u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6729846640u,        0.4459031899u)--(        0.6566661832u,        0.4195825532u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6328616584u,        0.3997731280u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6328616584u,        0.3997731280u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6729846640u,        0.4459031899u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6729846640u,        0.4459031899u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6566661832u,        0.4195825532u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6566661832u,        0.4195825532u) withcolor Green;
draw (        0.5865763239u,        0.4356951430u)--(        0.6158369583u,        0.4255519675u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6322410427u,        0.4518193379u)--(        0.6158369583u,        0.4255519675u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5865763239u,        0.4356951430u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5865763239u,        0.4356951430u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6322410427u,        0.4518193379u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6322410427u,        0.4518193379u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6158369583u,        0.4255519675u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6158369583u,        0.4255519675u) withcolor MidnightBlue;
draw (        0.5865763239u,        0.4356951430u)--(        0.5699896796u,        0.4618476182u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5904128523u,        0.4851277115u)--(        0.5699896796u,        0.4618476182u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5865763239u,        0.4356951430u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5865763239u,        0.4356951430u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5904128523u,        0.4851277115u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5904128523u,        0.4851277115u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5699896796u,        0.4618476182u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5699896796u,        0.4618476182u) withcolor Green;
draw (        0.6322410427u,        0.4518193379u)--(        0.6492881690u,        0.4259646630u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6729846640u,        0.4459031899u)--(        0.6492881690u,        0.4259646630u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6322410427u,        0.4518193379u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6322410427u,        0.4518193379u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6729846640u,        0.4459031899u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6729846640u,        0.4459031899u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6492881690u,        0.4259646630u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6492881690u,        0.4259646630u) withcolor MidnightBlue;
draw (        0.6322410427u,        0.4518193379u)--(        0.6210636351u,        0.4807007210u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5904128523u,        0.4851277115u)--(        0.6210636351u,        0.4807007210u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6322410427u,        0.4518193379u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6322410427u,        0.4518193379u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5904128523u,        0.4851277115u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5904128523u,        0.4851277115u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6210636351u,        0.4807007210u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6210636351u,        0.4807007210u) withcolor MidnightBlue;
draw (        0.6322410427u,        0.4518193379u)--(        0.6210636350u,        0.4807007209u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6459781071u,        0.4990946821u)--(        0.6210636350u,        0.4807007209u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6322410427u,        0.4518193379u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6322410427u,        0.4518193379u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6459781071u,        0.4990946821u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6459781071u,        0.4990946821u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6210636350u,        0.4807007209u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6210636350u,        0.4807007209u) withcolor MidnightBlue;
draw (        0.6729846640u,        0.4459031899u)--(        0.6750000000u,        0.4150000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6729846640u,        0.4459031899u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6729846640u,        0.4459031899u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6750000000u,        0.4150000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6750000000u,        0.4150000000u) withcolor Yellow;
draw (        0.6729846640u,        0.4459031899u)--(        0.6669092832u,        0.4762702523u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6459781071u,        0.4990946821u)--(        0.6669092832u,        0.4762702523u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6729846640u,        0.4459031899u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6729846640u,        0.4459031899u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6459781071u,        0.4990946821u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6459781071u,        0.4990946821u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6669092832u,        0.4762702523u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6669092832u,        0.4762702523u) withcolor MidnightBlue;
draw (        0.6729846640u,        0.4459031899u)--(        0.6669092835u,        0.4762702523u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6934321912u,        0.4922578756u)--(        0.6669092835u,        0.4762702523u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6729846640u,        0.4459031899u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6729846640u,        0.4459031899u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6934321912u,        0.4922578756u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6934321912u,        0.4922578756u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6669092835u,        0.4762702523u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6669092835u,        0.4762702523u) withcolor MidnightBlue;
draw (        0.5904128523u,        0.4851277115u)--(        0.6153273243u,        0.5035216728u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6459781071u,        0.4990946821u)--(        0.6153273243u,        0.5035216728u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5904128523u,        0.4851277115u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5904128523u,        0.4851277115u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6459781071u,        0.4990946821u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6459781071u,        0.4990946821u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6153273243u,        0.5035216728u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6153273243u,        0.5035216728u) withcolor MidnightBlue;
draw (        0.5904128523u,        0.4851277115u)--(        0.6210636351u,        0.4807007208u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6459781071u,        0.4990946821u)--(        0.6210636351u,        0.4807007208u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5904128523u,        0.4851277115u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5904128523u,        0.4851277115u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6459781071u,        0.4990946821u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6459781071u,        0.4990946821u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6210636351u,        0.4807007208u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6210636351u,        0.4807007208u) withcolor MidnightBlue;
draw (        0.5904128523u,        0.4851277115u)--(        0.6013544412u,        0.5140992555u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5802077452u,        0.5367241529u)--(        0.6013544412u,        0.5140992555u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5904128523u,        0.4851277115u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5904128523u,        0.4851277115u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5802077452u,        0.5367241529u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5802077452u,        0.5367241529u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6013544412u,        0.5140992555u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6013544412u,        0.5140992555u) withcolor MidnightBlue;
draw (        0.5904128523u,        0.4851277115u)--(        0.6013487165u,        0.5141014169u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6262631889u,        0.5324953778u)--(        0.6013487165u,        0.5141014169u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5904128523u,        0.4851277115u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5904128523u,        0.4851277115u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6262631889u,        0.5324953778u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6262631889u,        0.5324953778u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6013487165u,        0.5141014169u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6013487165u,        0.5141014169u) withcolor MidnightBlue;
draw (        0.5904128523u,        0.4851277115u)--(        0.6153273247u,        0.5035216724u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6262631889u,        0.5324953778u)--(        0.6153273247u,        0.5035216724u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5904128523u,        0.4851277115u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5904128523u,        0.4851277115u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6262631889u,        0.5324953778u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6262631889u,        0.5324953778u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6153273247u,        0.5035216724u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6153273247u,        0.5035216724u) withcolor MidnightBlue;
draw (        0.6459781071u,        0.4990946821u)--(        0.6725010149u,        0.5150823053u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6934321912u,        0.4922578756u)--(        0.6725010149u,        0.5150823053u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6459781071u,        0.4990946821u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6459781071u,        0.4990946821u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6934321912u,        0.4922578756u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6934321912u,        0.4922578756u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6725010149u,        0.5150823053u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6725010149u,        0.5150823053u) withcolor MidnightBlue;
draw (        0.6459781071u,        0.4990946821u)--(        0.6669092834u,        0.4762702525u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6934321912u,        0.4922578756u)--(        0.6669092834u,        0.4762702525u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6459781071u,        0.4990946821u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6459781071u,        0.4990946821u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6934321912u,        0.4922578756u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6934321912u,        0.4922578756u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6669092834u,        0.4762702525u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6669092834u,        0.4762702525u) withcolor MidnightBlue;
draw (        0.6459781071u,        0.4990946821u)--(        0.6153273243u,        0.5035216725u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6262631889u,        0.5324953778u)--(        0.6153273243u,        0.5035216725u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6459781071u,        0.4990946821u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6459781071u,        0.4990946821u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6262631889u,        0.5324953778u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6262631889u,        0.5324953778u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6153273243u,        0.5035216725u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6153273243u,        0.5035216725u) withcolor MidnightBlue;
draw (        0.6459781071u,        0.4990946821u)--(        0.6725010149u,        0.5150823053u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6828065393u,        0.5442861560u)--(        0.6725010149u,        0.5150823053u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6459781071u,        0.4990946821u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6459781071u,        0.4990946821u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6828065393u,        0.5442861560u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6828065393u,        0.5442861560u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6725010149u,        0.5150823053u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6725010149u,        0.5150823053u) withcolor MidnightBlue;
draw (        0.6934321912u,        0.4922578756u)--(        0.6725010149u,        0.5150823053u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6828065393u,        0.5442861560u)--(        0.6725010149u,        0.5150823053u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6934321912u,        0.4922578756u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6934321912u,        0.4922578756u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6828065393u,        0.5442861560u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6828065393u,        0.5442861560u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6725010149u,        0.5150823053u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6725010149u,        0.5150823053u) withcolor MidnightBlue;
draw (        0.6934321912u,        0.4922578756u)--(        0.7037377155u,        0.5214617263u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6828065393u,        0.5442861560u)--(        0.7037377155u,        0.5214617263u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6934321912u,        0.4922578756u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6934321912u,        0.4922578756u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6828065393u,        0.5442861560u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6828065393u,        0.5442861560u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7037377155u,        0.5214617263u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7037377155u,        0.5214617263u) withcolor MidnightBlue;
draw (        0.6934321912u,        0.4922578756u)--(        0.7037377156u,        0.5214617263u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7346495018u,        0.5195828424u)--(        0.7037377156u,        0.5214617263u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6934321912u,        0.4922578756u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6934321912u,        0.4922578756u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7346495018u,        0.5195828424u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7346495018u,        0.5195828424u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7037377156u,        0.5214617263u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7037377156u,        0.5214617263u) withcolor MidnightBlue;
draw (        0.6934321912u,        0.4922578756u)--(        0.7243439774u,        0.4903789917u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7346495018u,        0.5195828424u)--(        0.7243439774u,        0.4903789917u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6934321912u,        0.4922578756u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6934321912u,        0.4922578756u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7346495018u,        0.5195828424u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7346495018u,        0.5195828424u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7243439774u,        0.4903789917u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7243439774u,        0.4903789917u) withcolor Green;
draw (        0.6982155107u,        0.5872194905u)--(        0.7102292510u,        0.5586758612u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6828065393u,        0.5442861560u)--(        0.7102292510u,        0.5586758612u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6982155107u,        0.5872194905u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6982155107u,        0.5872194905u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6828065393u,        0.5442861560u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6828065393u,        0.5442861560u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7102292510u,        0.5586758612u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7102292510u,        0.5586758612u) withcolor MidnightBlue;
draw (        0.6982155107u,        0.5872194905u)--(        0.6691205169u,        0.5978284521u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6443253844u,        0.5792739342u)--(        0.6691205169u,        0.5978284521u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6982155107u,        0.5872194905u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6982155107u,        0.5872194905u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6443253844u,        0.5792739342u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6443253844u,        0.5792739342u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6691205169u,        0.5978284521u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6691205169u,        0.5978284521u) withcolor MidnightBlue;
draw (        0.6982155107u,        0.5872194905u)--(        0.7290300932u,        0.5841323821u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7410417328u,        0.5555878687u)--(        0.7290300932u,        0.5841323821u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6982155107u,        0.5872194905u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6982155107u,        0.5872194905u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7410417328u,        0.5555878687u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7410417328u,        0.5555878687u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7290300932u,        0.5841323821u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7290300932u,        0.5841323821u) withcolor MidnightBlue;
draw (        0.6982155107u,        0.5872194905u)--(        0.7102271504u,        0.5586749771u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7410417328u,        0.5555878687u)--(        0.7102271504u,        0.5586749771u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6982155107u,        0.5872194905u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6982155107u,        0.5872194905u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7410417328u,        0.5555878687u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7410417328u,        0.5555878687u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7102271504u,        0.5586749771u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7102271504u,        0.5586749771u) withcolor MidnightBlue;
draw (        0.6982155107u,        0.5872194905u)--(        0.7290299045u,        0.5841304997u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7557860761u,        0.5997246016u)--(        0.7290299045u,        0.5841304997u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6982155107u,        0.5872194905u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6982155107u,        0.5872194905u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7557860761u,        0.5997246016u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7557860761u,        0.5997246016u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7290299045u,        0.5841304997u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7290299045u,        0.5841304997u) withcolor MidnightBlue;
draw (        0.6982155107u,        0.5872194905u)--(        0.6691205169u,        0.5978284520u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6591014928u,        0.6271318287u)--(        0.6691205169u,        0.5978284520u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6982155107u,        0.5872194905u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6982155107u,        0.5872194905u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6591014928u,        0.6271318287u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6591014928u,        0.6271318287u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6691205169u,        0.5978284520u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6691205169u,        0.5978284520u) withcolor MidnightBlue;
draw (        0.6982155107u,        0.5872194905u)--(        0.6881964867u,        0.6165228672u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6591014928u,        0.6271318287u)--(        0.6881964867u,        0.6165228672u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6982155107u,        0.5872194905u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6982155107u,        0.5872194905u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6591014928u,        0.6271318287u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6591014928u,        0.6271318287u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6881964867u,        0.6165228672u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6881964867u,        0.6165228672u) withcolor MidnightBlue;
draw (        0.6982155107u,        0.5872194905u)--(        0.6881709481u,        0.6165141230u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7183405430u,        0.6235044217u)--(        0.6881709481u,        0.6165141230u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6982155107u,        0.5872194905u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6982155107u,        0.5872194905u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7183405430u,        0.6235044217u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7183405430u,        0.6235044217u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6881709481u,        0.6165141230u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6881709481u,        0.6165141230u) withcolor MidnightBlue;
draw (        0.5802077452u,        0.5367241529u)--(        0.6051189357u,        0.5551225582u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6262631889u,        0.5324953778u)--(        0.6051189357u,        0.5551225582u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5802077452u,        0.5367241529u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5802077452u,        0.5367241529u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6262631889u,        0.5324953778u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6262631889u,        0.5324953778u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6051189357u,        0.5551225582u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6051189357u,        0.5551225582u) withcolor MidnightBlue;
draw (        0.5802077452u,        0.5367241529u)--(        0.6013519983u,        0.5140969725u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6262631889u,        0.5324953778u)--(        0.6013519983u,        0.5140969725u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5802077452u,        0.5367241529u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5802077452u,        0.5367241529u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6262631889u,        0.5324953778u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6262631889u,        0.5324953778u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6013519983u,        0.5140969725u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6013519983u,        0.5140969725u) withcolor MidnightBlue;
draw (        0.5802077452u,        0.5367241529u)--(        0.5511689463u,        0.5474859785u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5432420020u,        0.5774231174u)--(        0.5511689463u,        0.5474859785u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5802077452u,        0.5367241529u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5802077452u,        0.5367241529u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5432420020u,        0.5774231174u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5432420020u,        0.5774231174u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5511689463u,        0.5474859785u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5511689463u,        0.5474859785u) withcolor Green;
draw (        0.5802077452u,        0.5367241529u)--(        0.5722808009u,        0.5666612917u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5432420020u,        0.5774231174u)--(        0.5722808009u,        0.5666612917u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5802077452u,        0.5367241529u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5802077452u,        0.5367241529u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5432420020u,        0.5774231174u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5432420020u,        0.5774231174u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5722808009u,        0.5666612917u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5722808009u,        0.5666612917u) withcolor MidnightBlue;
draw (        0.5802077452u,        0.5367241529u)--(        0.5722808009u,        0.5666612917u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5971919913u,        0.5850596971u)--(        0.5722808009u,        0.5666612917u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5802077452u,        0.5367241529u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5802077452u,        0.5367241529u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5971919913u,        0.5850596971u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5971919913u,        0.5850596971u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5722808009u,        0.5666612917u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5722808009u,        0.5666612917u) withcolor MidnightBlue;
draw (        0.5802077452u,        0.5367241529u)--(        0.6051189357u,        0.5551225582u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5971919913u,        0.5850596971u)--(        0.6051189357u,        0.5551225582u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5802077452u,        0.5367241529u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5802077452u,        0.5367241529u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5971919913u,        0.5850596971u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5971919913u,        0.5850596971u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6051189357u,        0.5551225582u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6051189357u,        0.5551225582u) withcolor MidnightBlue;
draw (        0.6262631889u,        0.5324953778u)--(        0.6522523288u,        0.5493367950u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6828065393u,        0.5442861560u)--(        0.6522523288u,        0.5493367950u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6262631889u,        0.5324953778u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6262631889u,        0.5324953778u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6828065393u,        0.5442861560u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6828065393u,        0.5442861560u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6522523288u,        0.5493367950u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6522523288u,        0.5493367950u) withcolor MidnightBlue;
draw (        0.6262631889u,        0.5324953778u)--(        0.6051189359u,        0.5551225583u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5971919913u,        0.5850596971u)--(        0.6051189359u,        0.5551225583u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6262631889u,        0.5324953778u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6262631889u,        0.5324953778u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5971919913u,        0.5850596971u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5971919913u,        0.5850596971u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6051189359u,        0.5551225583u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6051189359u,        0.5551225583u) withcolor MidnightBlue;
draw (        0.6262631889u,        0.5324953778u)--(        0.6183362443u,        0.5624325166u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5971919913u,        0.5850596971u)--(        0.6183362443u,        0.5624325166u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6262631889u,        0.5324953778u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6262631889u,        0.5324953778u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5971919913u,        0.5850596971u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5971919913u,        0.5850596971u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6183362443u,        0.5624325166u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6183362443u,        0.5624325166u) withcolor MidnightBlue;
draw (        0.6262631889u,        0.5324953778u)--(        0.6183362446u,        0.5624325167u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6443253844u,        0.5792739342u)--(        0.6183362446u,        0.5624325167u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6262631889u,        0.5324953778u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6262631889u,        0.5324953778u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6443253844u,        0.5792739342u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6443253844u,        0.5792739342u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6183362446u,        0.5624325167u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6183362446u,        0.5624325167u) withcolor MidnightBlue;
draw (        0.6262631889u,        0.5324953778u)--(        0.6522523287u,        0.5493367953u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6443253844u,        0.5792739342u)--(        0.6522523287u,        0.5493367953u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6262631889u,        0.5324953778u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6262631889u,        0.5324953778u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6443253844u,        0.5792739342u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6443253844u,        0.5792739342u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6522523287u,        0.5493367953u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6522523287u,        0.5493367953u) withcolor MidnightBlue;
draw (        0.6828065393u,        0.5442861560u)--(        0.7037377156u,        0.5214617264u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7346495018u,        0.5195828424u)--(        0.7037377156u,        0.5214617264u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6828065393u,        0.5442861560u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6828065393u,        0.5442861560u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7346495018u,        0.5195828424u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7346495018u,        0.5195828424u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7037377156u,        0.5214617264u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7037377156u,        0.5214617264u) withcolor MidnightBlue;
draw (        0.6828065393u,        0.5442861560u)--(        0.6522523289u,        0.5493367953u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6443253844u,        0.5792739342u)--(        0.6522523289u,        0.5493367953u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6828065393u,        0.5442861560u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6828065393u,        0.5442861560u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6443253844u,        0.5792739342u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6443253844u,        0.5792739342u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6522523289u,        0.5493367953u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6522523289u,        0.5493367953u) withcolor MidnightBlue;
draw (        0.6828065393u,        0.5442861560u)--(        0.7102275617u,        0.5586790800u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7410417328u,        0.5555878687u)--(        0.7102275617u,        0.5586790800u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6828065393u,        0.5442861560u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6828065393u,        0.5442861560u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7410417328u,        0.5555878687u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7410417328u,        0.5555878687u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7102275617u,        0.5586790800u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7102275617u,        0.5586790800u) withcolor MidnightBlue;
draw (        0.7346495018u,        0.5195828424u)--(        0.7624560343u,        0.5332160914u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7410417328u,        0.5555878687u)--(        0.7624560343u,        0.5332160914u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7346495018u,        0.5195828424u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7346495018u,        0.5195828424u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7410417328u,        0.5555878687u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7410417328u,        0.5555878687u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7624560343u,        0.5332160914u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7624560343u,        0.5332160914u) withcolor Green;
draw (        0.5432420020u,        0.5774231174u)--(        0.5722808009u,        0.5666612917u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5971919913u,        0.5850596971u)--(        0.5722808009u,        0.5666612917u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5432420020u,        0.5774231174u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5432420020u,        0.5774231174u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5971919913u,        0.5850596971u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5971919913u,        0.5850596971u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5722808009u,        0.5666612917u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5722808009u,        0.5666612917u) withcolor MidnightBlue;
draw (        0.5971919913u,        0.5850596971u)--(        0.6183362445u,        0.5624325168u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6443253844u,        0.5792739342u)--(        0.6183362445u,        0.5624325168u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5971919913u,        0.5850596971u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5971919913u,        0.5850596971u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6443253844u,        0.5792739342u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6443253844u,        0.5792739342u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6183362445u,        0.5624325168u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6183362445u,        0.5624325168u) withcolor MidnightBlue;
draw (        0.5971919913u,        0.5850596971u)--(        0.5859363120u,        0.6139106661u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5554474438u,        0.6193418400u)--(        0.5859363120u,        0.6139106661u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5971919913u,        0.5850596971u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5971919913u,        0.5850596971u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5554474438u,        0.6193418400u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5554474438u,        0.6193418400u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5859363120u,        0.6139106661u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5859363120u,        0.6139106661u) withcolor MidnightBlue;
draw (        0.5971919913u,        0.5850596971u)--(        0.5859363120u,        0.6139106661u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6122210394u,        0.6302869245u)--(        0.5859363120u,        0.6139106661u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5971919913u,        0.5850596971u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5971919913u,        0.5850596971u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6122210394u,        0.6302869245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6122210394u,        0.6302869245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5859363120u,        0.6139106661u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5859363120u,        0.6139106661u) withcolor MidnightBlue;
draw (        0.6443253844u,        0.5792739342u)--(        0.6343063606u,        0.6085773110u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6122210394u,        0.6302869245u)--(        0.6343063606u,        0.6085773110u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6443253844u,        0.5792739342u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6443253844u,        0.5792739342u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6122210394u,        0.6302869245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6122210394u,        0.6302869245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6343063606u,        0.6085773110u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6343063606u,        0.6085773110u) withcolor MidnightBlue;
draw (        0.6443253844u,        0.5792739342u)--(        0.6343063603u,        0.6085773109u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6591014928u,        0.6271318287u)--(        0.6343063603u,        0.6085773109u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6443253844u,        0.5792739342u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6443253844u,        0.5792739342u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6591014928u,        0.6271318287u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6591014928u,        0.6271318287u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6343063603u,        0.6085773109u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6343063603u,        0.6085773109u) withcolor MidnightBlue;
draw (        0.6443253844u,        0.5792739342u)--(        0.6691205169u,        0.5978284520u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6591014928u,        0.6271318287u)--(        0.6691205169u,        0.5978284520u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6443253844u,        0.5792739342u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6443253844u,        0.5792739342u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6591014928u,        0.6271318287u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6591014928u,        0.6271318287u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6691205169u,        0.5978284520u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6691205169u,        0.5978284520u) withcolor MidnightBlue;
draw (        0.7410417328u,        0.5555878687u)--(        0.7290290607u,        0.5841319476u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7557860761u,        0.5997246016u)--(        0.7290290607u,        0.5841319476u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7410417328u,        0.5555878687u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7410417328u,        0.5555878687u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7557860761u,        0.5997246016u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7557860761u,        0.5997246016u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7290290607u,        0.5841319476u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7290290607u,        0.5841319476u) withcolor MidnightBlue;
draw (        0.7557860761u,        0.5997246016u)--(        0.7850000000u,        0.6100000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7557860761u,        0.5997246016u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7557860761u,        0.5997246016u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7850000000u,        0.6100000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7850000000u,        0.6100000000u) withcolor Yellow;
draw (        0.5554474438u,        0.6193418400u)--(        0.5250000000u,        0.6250000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5554474438u,        0.6193418400u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5554474438u,        0.6193418400u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5250000000u,        0.6250000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5250000000u,        0.6250000000u) withcolor Yellow;
draw (        0.5554474438u,        0.6193418400u)--(        0.5817321713u,        0.6357180983u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6122210394u,        0.6302869245u)--(        0.5817321713u,        0.6357180983u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5554474438u,        0.6193418400u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5554474438u,        0.6193418400u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6122210394u,        0.6302869245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6122210394u,        0.6302869245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5817321713u,        0.6357180983u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5817321713u,        0.6357180983u) withcolor MidnightBlue;
draw (        0.5554474438u,        0.6193418400u)--(        0.5859363120u,        0.6139106662u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6122210394u,        0.6302869245u)--(        0.5859363120u,        0.6139106662u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5554474438u,        0.6193418400u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5554474438u,        0.6193418400u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6122210394u,        0.6302869245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6122210394u,        0.6302869245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5859363120u,        0.6139106662u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5859363120u,        0.6139106662u) withcolor MidnightBlue;
draw (        0.5554474438u,        0.6193418400u)--(        0.5817497770u,        0.6356898063u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5826625305u,        0.6666451873u)--(        0.5817497770u,        0.6356898063u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5554474438u,        0.6193418400u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5554474438u,        0.6193418400u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5826625305u,        0.6666451873u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5826625305u,        0.6666451873u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5817497770u,        0.6356898063u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5817497770u,        0.6356898063u) withcolor MidnightBlue;
draw (        0.6122210394u,        0.6302869245u)--(        0.6370161718u,        0.6488414425u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6591014928u,        0.6271318287u)--(        0.6370161718u,        0.6488414425u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6122210394u,        0.6302869245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6122210394u,        0.6302869245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6591014928u,        0.6271318287u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6591014928u,        0.6271318287u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6370161718u,        0.6488414425u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6370161718u,        0.6488414425u) withcolor MidnightBlue;
draw (        0.6122210394u,        0.6302869245u)--(        0.6343063604u,        0.6085773107u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6591014928u,        0.6271318287u)--(        0.6343063604u,        0.6085773107u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6122210394u,        0.6302869245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6122210394u,        0.6302869245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6591014928u,        0.6271318287u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6591014928u,        0.6271318287u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6343063604u,        0.6085773107u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6343063604u,        0.6085773107u) withcolor MidnightBlue;
draw (        0.6122210394u,        0.6302869245u)--(        0.5817272641u,        0.6356904783u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5826625305u,        0.6666451873u)--(        0.5817272641u,        0.6356904783u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6122210394u,        0.6302869245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6122210394u,        0.6302869245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5826625305u,        0.6666451873u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5826625305u,        0.6666451873u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5817272641u,        0.6356904783u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5817272641u,        0.6356904783u) withcolor MidnightBlue;
draw (        0.6122210394u,        0.6302869245u)--(        0.6131563059u,        0.6612416334u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5826625305u,        0.6666451873u)--(        0.6131563059u,        0.6612416334u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6122210394u,        0.6302869245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6122210394u,        0.6302869245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5826625305u,        0.6666451873u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5826625305u,        0.6666451873u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6131563059u,        0.6612416334u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6131563059u,        0.6612416334u) withcolor MidnightBlue;
draw (        0.6122210394u,        0.6302869245u)--(        0.6131563061u,        0.6612416334u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6379881712u,        0.6797469625u)--(        0.6131563061u,        0.6612416334u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6122210394u,        0.6302869245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6122210394u,        0.6302869245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6379881712u,        0.6797469625u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6379881712u,        0.6797469625u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6131563061u,        0.6612416334u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6131563061u,        0.6612416334u) withcolor MidnightBlue;
draw (        0.6122210394u,        0.6302869245u)--(        0.6370529046u,        0.6487922535u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6379881712u,        0.6797469625u)--(        0.6370529046u,        0.6487922535u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6122210394u,        0.6302869245u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6122210394u,        0.6302869245u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6379881712u,        0.6797469625u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6379881712u,        0.6797469625u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6370529046u,        0.6487922535u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6370529046u,        0.6487922535u) withcolor MidnightBlue;
draw (        0.6591014928u,        0.6271318287u)--(        0.6892624049u,        0.6341594964u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7183405430u,        0.6235044217u)--(        0.6892624049u,        0.6341594964u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6591014928u,        0.6271318287u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6591014928u,        0.6271318287u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7183405430u,        0.6235044217u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7183405430u,        0.6235044217u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6892624049u,        0.6341594964u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6892624049u,        0.6341594964u) withcolor MidnightBlue;
draw (        0.6591014928u,        0.6271318287u)--(        0.6881796309u,        0.6164767540u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7183405430u,        0.6235044217u)--(        0.6881796309u,        0.6164767540u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6591014928u,        0.6271318287u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6591014928u,        0.6271318287u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7183405430u,        0.6235044217u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7183405430u,        0.6235044217u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6881796309u,        0.6164767540u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6881796309u,        0.6164767540u) withcolor MidnightBlue;
draw (        0.6591014928u,        0.6271318287u)--(        0.6369704723u,        0.6487948540u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6379881712u,        0.6797469625u)--(        0.6369704723u,        0.6487948540u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6591014928u,        0.6271318287u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6591014928u,        0.6271318287u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6379881712u,        0.6797469625u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6379881712u,        0.6797469625u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6369704723u,        0.6487948540u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6369704723u,        0.6487948540u) withcolor MidnightBlue;
draw (        0.6591014928u,        0.6271318287u)--(        0.6892624049u,        0.6341594964u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6847608802u,        0.6647994214u)--(        0.6892624049u,        0.6341594964u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6591014928u,        0.6271318287u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6591014928u,        0.6271318287u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6847608802u,        0.6647994214u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6847608802u,        0.6647994214u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6892624049u,        0.6341594964u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6892624049u,        0.6341594964u) withcolor MidnightBlue;
draw (        0.7183405430u,        0.6235044217u)--(        0.6892624049u,        0.6341594964u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6847608802u,        0.6647994214u)--(        0.6892624049u,        0.6341594964u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7183405430u,        0.6235044217u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7183405430u,        0.6235044217u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6847608802u,        0.6647994214u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6847608802u,        0.6647994214u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6892624049u,        0.6341594964u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6892624049u,        0.6341594964u) withcolor MidnightBlue;
draw (        0.7183405430u,        0.6235044217u)--(        0.7138390183u,        0.6541443467u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6847608802u,        0.6647994214u)--(        0.7138390183u,        0.6541443467u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.7183405430u,        0.6235044217u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7183405430u,        0.6235044217u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6847608802u,        0.6647994214u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6847608802u,        0.6647994214u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.7138390183u,        0.6541443467u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.7138390183u,        0.6541443467u) withcolor Green;
draw (        0.5389761541u,        0.6573459103u)--(        0.5563523779u,        0.6829805660u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5826625305u,        0.6666451873u)--(        0.5563523779u,        0.6829805660u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5389761541u,        0.6573459103u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5389761541u,        0.6573459103u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5826625305u,        0.6666451873u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5826625305u,        0.6666451873u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5563523779u,        0.6829805660u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5563523779u,        0.6829805660u) withcolor MidnightBlue;
draw (        0.5389761541u,        0.6573459103u)--(        0.5381448428u,        0.6883035854u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5555210667u,        0.7139382412u)--(        0.5381448428u,        0.6883035854u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5389761541u,        0.6573459103u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5389761541u,        0.6573459103u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5555210667u,        0.7139382412u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5555210667u,        0.7139382412u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5381448428u,        0.6883035854u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5381448428u,        0.6883035854u) withcolor Green;
draw (        0.5389761541u,        0.6573459103u)--(        0.5563523780u,        0.6829805660u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5555210667u,        0.7139382412u)--(        0.5563523780u,        0.6829805660u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5389761541u,        0.6573459103u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5389761541u,        0.6573459103u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5555210667u,        0.7139382412u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5555210667u,        0.7139382412u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5563523780u,        0.6829805660u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5563523780u,        0.6829805660u) withcolor MidnightBlue;
draw (        0.5826625305u,        0.6666451873u)--(        0.6131563059u,        0.6612416336u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6379881712u,        0.6797469625u)--(        0.6131563059u,        0.6612416336u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5826625305u,        0.6666451873u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5826625305u,        0.6666451873u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6379881712u,        0.6797469625u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6379881712u,        0.6797469625u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6131563059u,        0.6612416336u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6131563059u,        0.6612416336u) withcolor MidnightBlue;
draw (        0.5826625305u,        0.6666451873u)--(        0.5563523779u,        0.6829805660u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5555210667u,        0.7139382412u)--(        0.5563523779u,        0.6829805660u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5826625305u,        0.6666451873u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5826625305u,        0.6666451873u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5555210667u,        0.7139382412u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5555210667u,        0.7139382412u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5563523779u,        0.6829805660u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5563523779u,        0.6829805660u) withcolor MidnightBlue;
draw (        0.5826625305u,        0.6666451873u)--(        0.5818312194u,        0.6976028624u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5555210667u,        0.7139382412u)--(        0.5818312194u,        0.6976028624u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5826625305u,        0.6666451873u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5826625305u,        0.6666451873u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5555210667u,        0.7139382412u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5555210667u,        0.7139382412u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5818312194u,        0.6976028624u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5818312194u,        0.6976028624u) withcolor MidnightBlue;
draw (        0.5826625305u,        0.6666451873u)--(        0.5818188971u,        0.6976025291u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6085728756u,        0.7132003932u)--(        0.5818188971u,        0.6976025291u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5826625305u,        0.6666451873u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5826625305u,        0.6666451873u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6085728756u,        0.7132003932u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6085728756u,        0.7132003932u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5818188971u,        0.6976025291u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5818188971u,        0.6976025291u) withcolor MidnightBlue;
draw (        0.6379881712u,        0.6797469625u)--(        0.6671204848u,        0.6902530102u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6847608802u,        0.6647994214u)--(        0.6671204848u,        0.6902530102u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6379881712u,        0.6797469625u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6379881712u,        0.6797469625u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6847608802u,        0.6647994214u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6847608802u,        0.6647994214u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6671204848u,        0.6902530102u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6671204848u,        0.6902530102u) withcolor Green;
draw (        0.5555210667u,        0.7139382412u)--(        0.5822691781u,        0.7295461643u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6085728756u,        0.7132003932u)--(        0.5822691781u,        0.7295461643u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5555210667u,        0.7139382412u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5555210667u,        0.7139382412u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6085728756u,        0.7132003932u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6085728756u,        0.7132003932u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5822691781u,        0.7295461643u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5822691781u,        0.7295461643u) withcolor Green;
draw (        0.5555210667u,        0.7139382412u)--(        0.5818247641u,        0.6975924701u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6085728756u,        0.7132003932u)--(        0.5818247641u,        0.6975924701u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5555210667u,        0.7139382412u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5555210667u,        0.7139382412u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6085728756u,        0.7132003932u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6085728756u,        0.7132003932u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5818247641u,        0.6975924701u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5818247641u,        0.6975924701u) withcolor MidnightBlue;
draw (        0.5926306482u,        0.9706954877u)--(        0.6066984826u,        0.9982847023u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5831951704u,        1.0184505906u)--(        0.6066984826u,        0.9982847023u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5926306482u,        0.9706954877u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5926306482u,        0.9706954877u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5831951704u,        1.0184505906u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5831951704u,        1.0184505906u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6066984826u,        0.9982847023u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6066984826u,        0.9982847023u) withcolor MidnightBlue;
draw (        0.5926306482u,        0.9706954877u)--(        0.6066984825u,        0.9982847024u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6375703128u,        1.0007339529u)--(        0.6066984825u,        0.9982847024u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5926306482u,        0.9706954877u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5926306482u,        0.9706954877u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6375703128u,        1.0007339529u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6375703128u,        1.0007339529u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6066984825u,        0.9982847024u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6066984825u,        0.9982847024u) withcolor MidnightBlue;
draw (        0.5296621648u,        1.0158876503u)--(        0.5600000000u,        1.0100000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5296621648u,        1.0158876503u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5296621648u,        1.0158876503u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5600000000u,        1.0100000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5600000000u,        1.0100000000u) withcolor Yellow;
draw (        0.5831951704u,        1.0184505906u)--(        0.6140670007u,        1.0208998412u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6375703128u,        1.0007339529u)--(        0.6140670007u,        1.0208998412u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5831951704u,        1.0184505906u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5831951704u,        1.0184505906u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6375703128u,        1.0007339529u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6375703128u,        1.0007339529u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6140670007u,        1.0208998412u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6140670007u,        1.0208998412u) withcolor MidnightBlue;
draw (        0.5831951704u,        1.0184505906u)--(        0.6066984825u,        0.9982847023u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6375703128u,        1.0007339529u)--(        0.6066984825u,        0.9982847023u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5831951704u,        1.0184505906u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5831951704u,        1.0184505906u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6375703128u,        1.0007339529u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6375703128u,        1.0007339529u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6066984825u,        0.9982847023u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6066984825u,        0.9982847023u) withcolor MidnightBlue;
draw (        0.5831951704u,        1.0184505906u)--(        0.6140677984u,        1.0208897657u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5967242763u,        1.0465465577u)--(        0.6140677984u,        1.0208897657u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5831951704u,        1.0184505906u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5831951704u,        1.0184505906u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5967242763u,        1.0465465577u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5967242763u,        1.0465465577u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6140677984u,        1.0208897657u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6140677984u,        1.0208897657u) withcolor MidnightBlue;
draw (        0.6375703128u,        1.0007339529u)--(        0.6140454015u,        1.0208746401u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5967242763u,        1.0465465577u)--(        0.6140454015u,        1.0208746401u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6375703128u,        1.0007339529u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6375703128u,        1.0007339529u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5967242763u,        1.0465465577u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5967242763u,        1.0465465577u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6140454015u,        1.0208746401u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6140454015u,        1.0208746401u) withcolor MidnightBlue;
draw (        0.6375703128u,        1.0007339529u)--(        0.6202491876u,        1.0264058705u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5967242763u,        1.0465465577u)--(        0.6202491876u,        1.0264058705u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6375703128u,        1.0007339529u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6375703128u,        1.0007339529u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5967242763u,        1.0465465577u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5967242763u,        1.0465465577u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6202491876u,        1.0264058705u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6202491876u,        1.0264058705u) withcolor MidnightBlue;
draw (        0.6375703128u,        1.0007339529u)--(        0.6203414882u,        1.0264679058u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6228165214u,        1.0573376798u)--(        0.6203414882u,        1.0264679058u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.6375703128u,        1.0007339529u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6375703128u,        1.0007339529u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6228165214u,        1.0573376798u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6228165214u,        1.0573376798u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6203414882u,        1.0264679058u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6203414882u,        1.0264679058u) withcolor MidnightBlue;
draw (        0.5967242763u,        1.0465465577u)--(        0.6203046884u,        1.0264708784u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6228165214u,        1.0573376798u)--(        0.6203046884u,        1.0264708784u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5967242763u,        1.0465465577u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5967242763u,        1.0465465577u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6228165214u,        1.0573376798u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6228165214u,        1.0573376798u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.6203046884u,        1.0264708784u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.6203046884u,        1.0264708784u) withcolor MidnightBlue;
draw (        0.5943425459u,        1.0893154275u)--(        0.5900000000u,        1.1200000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 4.0 shifted (        0.5943425459u,        1.0893154275u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5943425459u,        1.0893154275u) withcolor 0.9red;
draw fullcircle scaled 4.0 shifted (        0.5900000000u,        1.1200000000u) withcolor 0.4white;
fill fullcircle scaled 4.0 shifted (        0.5900000000u,        1.1200000000u) withcolor Yellow;
% r =         0.0309688348
endfig;
end;
