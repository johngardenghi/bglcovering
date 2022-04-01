prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(40);
u := 20cm;
path pols[];
pair p[];
R := .0921104165u;
D := 2*R;

pols[1] = 
(        0.0000000000u,        0.6666666667u)--
(        0.0000000000u,        0.4666666667u)--
(        0.1333333333u,        0.4666666667u)--
(        0.1333333333u,        0.6666666667u)--
(        0.0000000000u,        0.6666666667u)--
cycle;
pols[2] = 
(        0.2333333333u,        0.3333333333u)--
(        0.1333333333u,        0.4666666667u)--
(        0.0000000000u,        0.4666666667u)--
(        0.0000000000u,        0.2000000000u)--
(        0.1333333333u,        0.2000000000u)--
(        0.2333333333u,        0.3333333333u)--
cycle;
pols[3] = 
(        0.0000000000u,        0.2000000000u)--
(        0.0000000000u,       -0.2666666667u)--
(        0.1333333333u,       -0.2666666667u)--
(        0.1333333333u,        0.2000000000u)--
(        0.0000000000u,        0.2000000000u)--
cycle;
pols[4] = 
(        0.0000000000u,       -0.2666666667u)--
(        0.0000000000u,       -0.3333333333u)--
(        0.1333333333u,       -0.3333333333u)--
(        0.1333333333u,       -0.2666666667u)--
(        0.0000000000u,       -0.2666666667u)--
cycle;
pols[5] = 
(        0.1333333333u,        0.6666666667u)--
(        0.1333333333u,        0.4666666667u)--
(        0.4666666667u,        0.4666666667u)--
(        0.4666666667u,        0.6666666667u)--
(        0.1333333333u,        0.6666666667u)--
cycle;
pols[6] = 
(        0.4666666667u,        0.4666666667u)--
(        0.5333333333u,        0.4666666667u)--
(        0.6000000000u,        0.6666666667u)--
(        0.4666666667u,        0.6666666667u)--
(        0.4666666667u,        0.4666666667u)--
cycle;
pols[7] = 
(        0.5333333333u,        0.4666666667u)--
(        0.4666666667u,        0.4666666667u)--
(        0.3000000000u,        0.3333333333u)--
(        0.4666666667u,        0.2000000000u)--
(        0.6000000000u,        0.3333333333u)--
(        0.5333333333u,        0.4666666667u)--
cycle;
pols[8] = 
(        0.4666666667u,        0.2000000000u)--
(        0.1333333333u,        0.2000000000u)--
(        0.1333333333u,       -0.2666666667u)--
(        0.4000000000u,        0.0000000000u)--
(        0.4666666667u,        0.2000000000u)--
cycle;
pols[9] = 
(        0.7333333333u,        0.1333333333u)--
(        0.6000000000u,        0.3333333333u)--
(        0.4666666667u,        0.2000000000u)--
(        0.4000000000u,        0.0000000000u)--
(        0.6000000000u,        0.0000000000u)--
(        0.7333333333u,        0.1333333333u)--
cycle;
pols[10] = 
(        0.8666666667u,       -0.2666666667u)--
(        0.8666666667u,       -0.3333333333u)--
(        1.0000000000u,       -0.3333333333u)--
(        1.0000000000u,       -0.2666666667u)--
(        0.8666666667u,       -0.2666666667u)--
cycle;
pols[11] = 
(        0.8666666667u,        0.3333333333u)--
(        0.7333333333u,        0.1333333333u)--
(        0.8666666667u,       -0.2666666667u)--
(        1.0000000000u,       -0.2666666667u)--
(        1.0000000000u,        0.3333333333u)--
(        0.8666666667u,        0.3333333333u)--
cycle;
pols[12] = 
(        0.8666666667u,        0.6666666667u)--
(        0.8000000000u,        0.6666666667u)--
(        0.7333333333u,        0.5333333333u)--
(        0.8666666667u,        0.3333333333u)--
(        1.0000000000u,        0.3333333333u)--
(        1.0000000000u,        0.5333333333u)--
(        0.8666666667u,        0.6666666667u)--
cycle;
pols[13] = 
(        0.7333333333u,        0.5333333333u)--
(        0.8000000000u,        0.6666666667u)--
(        0.5333333333u,        0.4666666667u)--
(        0.6000000000u,        0.3333333333u)--
(        0.7333333333u,        0.5333333333u)--
cycle;
pols[14] = 
(        0.7333333333u,        0.1333333333u)--
(        0.6000000000u,        0.0000000000u)--
(        0.8666666667u,       -0.2666666667u)--
(        0.7333333333u,        0.1333333333u)--
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
draw (        0.0000000000u,        0.4666666667u)--(        0.1333333333u,        0.4666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1333333333u,        0.4666666667u)--(        0.1333333333u,        0.6666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1333333333u,        0.4666666667u)--(        0.0000000000u,        0.4666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.0000000000u,        0.2000000000u)--(        0.1333333333u,        0.2000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.0000000000u,       -0.2666666667u)--(        0.1333333333u,       -0.2666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1333333333u,       -0.2666666667u)--(        0.1333333333u,        0.2000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1333333333u,        0.2000000000u)--(        0.0000000000u,        0.2000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1333333333u,       -0.2666666667u)--(        0.0000000000u,       -0.2666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1333333333u,        0.6666666667u)--(        0.1333333333u,        0.4666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4666666667u,        0.4666666667u)--(        0.4666666667u,        0.6666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4666666667u,        0.4666666667u)--(        0.5333333333u,        0.4666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4666666667u,        0.6666666667u)--(        0.4666666667u,        0.4666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5333333333u,        0.4666666667u)--(        0.4666666667u,        0.4666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4666666667u,        0.2000000000u)--(        0.6000000000u,        0.3333333333u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.6000000000u,        0.3333333333u)--(        0.5333333333u,        0.4666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1333333333u,        0.2000000000u)--(        0.1333333333u,       -0.2666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4000000000u,        0.0000000000u)--(        0.4666666667u,        0.2000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.6000000000u,        0.3333333333u)--(        0.4666666667u,        0.2000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.4666666667u,        0.2000000000u)--(        0.4000000000u,        0.0000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.6000000000u,        0.0000000000u)--(        0.7333333333u,        0.1333333333u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        1.0000000000u,       -0.2666666667u)--(        0.8666666667u,       -0.2666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.7333333333u,        0.1333333333u)--(        0.8666666667u,       -0.2666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.8666666667u,       -0.2666666667u)--(        1.0000000000u,       -0.2666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        1.0000000000u,        0.3333333333u)--(        0.8666666667u,        0.3333333333u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.8000000000u,        0.6666666667u)--(        0.7333333333u,        0.5333333333u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.8666666667u,        0.3333333333u)--(        1.0000000000u,        0.3333333333u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.7333333333u,        0.5333333333u)--(        0.8000000000u,        0.6666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5333333333u,        0.4666666667u)--(        0.6000000000u,        0.3333333333u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.7333333333u,        0.1333333333u)--(        0.6000000000u,        0.0000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.8666666667u,       -0.2666666667u)--(        0.7333333333u,        0.1333333333u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.0000000000u,        0.6666666667u)--(        0.0000000000u,        0.4666666667u) withpen pencircle scaled 1.0;
draw (        0.1333333333u,        0.6666666667u)--(        0.0000000000u,        0.6666666667u) withpen pencircle scaled 1.0;
draw (        0.2333333333u,        0.3333333333u)--(        0.1333333333u,        0.4666666667u) withpen pencircle scaled 1.0;
draw (        0.0000000000u,        0.4666666667u)--(        0.0000000000u,        0.2000000000u) withpen pencircle scaled 1.0;
draw (        0.1333333333u,        0.2000000000u)--(        0.2333333333u,        0.3333333333u) withpen pencircle scaled 1.0;
draw (        0.0000000000u,        0.2000000000u)--(        0.0000000000u,       -0.2666666667u) withpen pencircle scaled 1.0;
draw (        0.0000000000u,       -0.2666666667u)--(        0.0000000000u,       -0.3333333333u) withpen pencircle scaled 1.0;
draw (        0.0000000000u,       -0.3333333333u)--(        0.1333333333u,       -0.3333333333u) withpen pencircle scaled 1.0;
draw (        0.1333333333u,       -0.3333333333u)--(        0.1333333333u,       -0.2666666667u) withpen pencircle scaled 1.0;
draw (        0.1333333333u,        0.4666666667u)--(        0.4666666667u,        0.4666666667u) withpen pencircle scaled 1.0;
draw (        0.4666666667u,        0.6666666667u)--(        0.1333333333u,        0.6666666667u) withpen pencircle scaled 1.0;
draw (        0.5333333333u,        0.4666666667u)--(        0.6000000000u,        0.6666666667u) withpen pencircle scaled 1.0;
draw (        0.6000000000u,        0.6666666667u)--(        0.4666666667u,        0.6666666667u) withpen pencircle scaled 1.0;
draw (        0.4666666667u,        0.4666666667u)--(        0.3000000000u,        0.3333333333u) withpen pencircle scaled 1.0;
draw (        0.3000000000u,        0.3333333333u)--(        0.4666666667u,        0.2000000000u) withpen pencircle scaled 1.0;
draw (        0.4666666667u,        0.2000000000u)--(        0.1333333333u,        0.2000000000u) withpen pencircle scaled 1.0;
draw (        0.1333333333u,       -0.2666666667u)--(        0.4000000000u,        0.0000000000u) withpen pencircle scaled 1.0;
draw (        0.7333333333u,        0.1333333333u)--(        0.6000000000u,        0.3333333333u) withpen pencircle scaled 1.0;
draw (        0.4000000000u,        0.0000000000u)--(        0.6000000000u,        0.0000000000u) withpen pencircle scaled 1.0;
draw (        0.8666666667u,       -0.2666666667u)--(        0.8666666667u,       -0.3333333333u) withpen pencircle scaled 1.0;
draw (        0.8666666667u,       -0.3333333333u)--(        1.0000000000u,       -0.3333333333u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,       -0.3333333333u)--(        1.0000000000u,       -0.2666666667u) withpen pencircle scaled 1.0;
draw (        0.8666666667u,        0.3333333333u)--(        0.7333333333u,        0.1333333333u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,       -0.2666666667u)--(        1.0000000000u,        0.3333333333u) withpen pencircle scaled 1.0;
draw (        0.8666666667u,        0.6666666667u)--(        0.8000000000u,        0.6666666667u) withpen pencircle scaled 1.0;
draw (        0.7333333333u,        0.5333333333u)--(        0.8666666667u,        0.3333333333u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        0.3333333333u)--(        1.0000000000u,        0.5333333333u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        0.5333333333u)--(        0.8666666667u,        0.6666666667u) withpen pencircle scaled 1.0;
draw (        0.8000000000u,        0.6666666667u)--(        0.5333333333u,        0.4666666667u) withpen pencircle scaled 1.0;
draw (        0.6000000000u,        0.3333333333u)--(        0.7333333333u,        0.5333333333u) withpen pencircle scaled 1.0;
draw (        0.6000000000u,        0.0000000000u)--(        0.8666666667u,       -0.2666666667u) withpen pencircle scaled 1.0;

x1 = .9333239180u;
y1 = -.2697618589u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .8195635260u;
y2 = -.1560441004u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .9359789014u;
y3 = -.1399625099u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .8586989733u;
y4 = .0139072118u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .9906224614u;
y5 = .0178940348u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .7215833170u;
y6 = -.0392887795u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .9313924961u;
y7 = .1709868359u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .9743991409u;
y8 = .3209813535u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .0666796245u;
y9 = -.2697584450u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = .2020772361u;
y10 = -.1652950001u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;

x11 = .4574282374u;
y11 = .0424733375u;
draw fullcircle scaled D shifted z11 withcolor RoyalBlue withpen pencircle scaled 1.0;

x12 = .6092389919u;
y12 = .0598086975u;
draw fullcircle scaled D shifted z12 withcolor RoyalBlue withpen pencircle scaled 1.0;

x13 = .7742760417u;
y13 = .1061442564u;
draw fullcircle scaled D shifted z13 withcolor RoyalBlue withpen pencircle scaled 1.0;

x14 = .8416070986u;
y14 = .2634880536u;
draw fullcircle scaled D shifted z14 withcolor RoyalBlue withpen pencircle scaled 1.0;

x15 = .0577516326u;
y15 = -.1344255279u;
draw fullcircle scaled D shifted z15 withcolor RoyalBlue withpen pencircle scaled 1.0;

x16 = .3023505295u;
y16 = -.0288202827u;
draw fullcircle scaled D shifted z16 withcolor RoyalBlue withpen pencircle scaled 1.0;

x17 = .3435654252u;
y17 = .1278329339u;
draw fullcircle scaled D shifted z17 withcolor RoyalBlue withpen pencircle scaled 1.0;

x18 = .4920438508u;
y18 = .1872956230u;
draw fullcircle scaled D shifted z18 withcolor RoyalBlue withpen pencircle scaled 1.0;

x19 = .6471686451u;
y19 = .2004806763u;
draw fullcircle scaled D shifted z19 withcolor RoyalBlue withpen pencircle scaled 1.0;

x20 = .8520557521u;
y20 = .4304402320u;
draw fullcircle scaled D shifted z20 withcolor RoyalBlue withpen pencircle scaled 1.0;

x21 = .1623713862u;
y21 = -.0062409618u;
draw fullcircle scaled D shifted z21 withcolor RoyalBlue withpen pencircle scaled 1.0;

x22 = .2251209601u;
y22 = .1311375131u;
draw fullcircle scaled D shifted z22 withcolor RoyalBlue withpen pencircle scaled 1.0;

x23 = .5456882242u;
y23 = .3311338903u;
draw fullcircle scaled D shifted z23 withcolor RoyalBlue withpen pencircle scaled 1.0;

x24 = .9666841951u;
y24 = .4953530208u;
draw fullcircle scaled D shifted z24 withcolor RoyalBlue withpen pencircle scaled 1.0;

x25 = .0209332224u;
y25 = .0270359257u;
draw fullcircle scaled D shifted z25 withcolor RoyalBlue withpen pencircle scaled 1.0;

x26 = .0845830960u;
y26 = .1532242846u;
draw fullcircle scaled D shifted z26 withcolor RoyalBlue withpen pencircle scaled 1.0;

x27 = .3916001728u;
y27 = .3234608684u;
draw fullcircle scaled D shifted z27 withcolor RoyalBlue withpen pencircle scaled 1.0;

x28 = .4736106869u;
y28 = .4693904663u;
draw fullcircle scaled D shifted z28 withcolor RoyalBlue withpen pencircle scaled 1.0;

x29 = .6346281184u;
y29 = .4659104732u;
draw fullcircle scaled D shifted z29 withcolor RoyalBlue withpen pencircle scaled 1.0;

x30 = .7445701922u;
y30 = .5636398292u;
draw fullcircle scaled D shifted z30 withcolor RoyalBlue withpen pencircle scaled 1.0;

x31 = .8576592976u;
y31 = .6045019736u;
draw fullcircle scaled D shifted z31 withcolor RoyalBlue withpen pencircle scaled 1.0;

x32 = .0293107762u;
y32 = .2770314033u;
draw fullcircle scaled D shifted z32 withcolor RoyalBlue withpen pencircle scaled 1.0;

x33 = .1459178650u;
y33 = .3119227870u;
draw fullcircle scaled D shifted z33 withcolor RoyalBlue withpen pencircle scaled 1.0;

x34 = .3100245403u;
y34 = .5245546333u;
draw fullcircle scaled D shifted z34 withcolor RoyalBlue withpen pencircle scaled 1.0;

x35 = .5283593199u;
y35 = .6087586894u;
draw fullcircle scaled D shifted z35 withcolor RoyalBlue withpen pencircle scaled 1.0;

x36 = .0384789954u;
y36 = .4480426299u;
draw fullcircle scaled D shifted z36 withcolor RoyalBlue withpen pencircle scaled 1.0;

x37 = .1585769498u;
y37 = .4812356933u;
draw fullcircle scaled D shifted z37 withcolor RoyalBlue withpen pencircle scaled 1.0;

x38 = .2137582182u;
y38 = .6404984441u;
draw fullcircle scaled D shifted z38 withcolor RoyalBlue withpen pencircle scaled 1.0;

x39 = .3794058650u;
y39 = .6166078773u;
draw fullcircle scaled D shifted z39 withcolor RoyalBlue withpen pencircle scaled 1.0;

x40 = .0627186792u;
y40 = .5992038608u;
draw fullcircle scaled D shifted z40 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.9333239180u,       -0.2697618589u)--(        0.8666666667u,       -0.3333333333u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9333239180u,       -0.2697618589u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9333239180u,       -0.2697618589u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8666666667u,       -0.3333333333u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8666666667u,       -0.3333333333u) withcolor Yellow;
draw (        0.9333239180u,       -0.2697618589u)--(        1.0000000000u,       -0.3333333333u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9333239180u,       -0.2697618589u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9333239180u,       -0.2697618589u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,       -0.3333333333u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,       -0.3333333333u) withcolor Yellow;
draw (        0.9333239180u,       -0.2697618589u)--(        0.8447003423u,       -0.2446582601u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8195635260u,       -0.1560441004u)--(        0.8447003423u,       -0.2446582601u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9333239180u,       -0.2697618589u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9333239180u,       -0.2697618589u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8195635260u,       -0.1560441004u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8195635260u,       -0.1560441004u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8447003423u,       -0.2446582601u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8447003423u,       -0.2446582601u) withcolor Green;
draw (        0.9333239180u,       -0.2697618589u)--(        0.9999875623u,       -0.2061986040u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9359789014u,       -0.1399625099u)--(        0.9999875623u,       -0.2061986040u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9333239180u,       -0.2697618589u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9333239180u,       -0.2697618589u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9359789014u,       -0.1399625099u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9359789014u,       -0.1399625099u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9999875623u,       -0.2061986040u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9999875623u,       -0.2061986040u) withcolor Green;
draw (        0.8195635260u,       -0.1560441004u)--(        0.8680646746u,       -0.0777372206u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9359789014u,       -0.1399625099u)--(        0.8680646746u,       -0.0777372206u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8195635260u,       -0.1560441004u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8195635260u,       -0.1560441004u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9359789014u,       -0.1399625099u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9359789014u,       -0.1399625099u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8680646746u,       -0.0777372206u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8680646746u,       -0.0777372206u) withcolor MidnightBlue;
draw (        0.8195635260u,       -0.1560441004u)--(        0.8102136296u,       -0.0644094548u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8586989733u,        0.0139072118u)--(        0.8102136296u,       -0.0644094548u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8195635260u,       -0.1560441004u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8195635260u,       -0.1560441004u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8586989733u,        0.0139072118u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8586989733u,        0.0139072118u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8102136296u,       -0.0644094548u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8102136296u,       -0.0644094548u) withcolor MidnightBlue;
draw (        0.8195635260u,       -0.1560441004u)--(        0.8680488698u,       -0.0777274337u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8586989733u,        0.0139072118u)--(        0.8680488698u,       -0.0777274337u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8195635260u,       -0.1560441004u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8195635260u,       -0.1560441004u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8586989733u,        0.0139072118u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8586989733u,        0.0139072118u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8680488698u,       -0.0777274337u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8680488698u,       -0.0777274337u) withcolor MidnightBlue;
draw (        0.8195635260u,       -0.1560441004u)--(        0.7309451343u,       -0.1309222080u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7215833170u,       -0.0392887795u)--(        0.7309451343u,       -0.1309222080u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8195635260u,       -0.1560441004u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8195635260u,       -0.1560441004u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7215833170u,       -0.0392887795u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7215833170u,       -0.0392887795u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7309451343u,       -0.1309222080u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7309451343u,       -0.1309222080u) withcolor Green;
draw (        0.8195635260u,       -0.1560441004u)--(        0.8102017087u,       -0.0644106719u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7215833170u,       -0.0392887795u)--(        0.8102017087u,       -0.0644106719u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8195635260u,       -0.1560441004u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8195635260u,       -0.1560441004u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7215833170u,       -0.0392887795u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7215833170u,       -0.0392887795u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8102017087u,       -0.0644106719u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8102017087u,       -0.0644106719u) withcolor MidnightBlue;
draw (        0.9359789014u,       -0.1399625099u)--(        0.8680762028u,       -0.0777246407u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8586989733u,        0.0139072118u)--(        0.8680762028u,       -0.0777246407u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9359789014u,       -0.1399625099u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9359789014u,       -0.1399625099u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8586989733u,        0.0139072118u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8586989733u,        0.0139072118u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8680762028u,       -0.0777246407u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8680762028u,       -0.0777246407u) withcolor MidnightBlue;
draw (        0.9359789014u,       -0.1399625099u)--(        0.9266016719u,       -0.0483306574u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8586989733u,        0.0139072118u)--(        0.9266016719u,       -0.0483306574u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9359789014u,       -0.1399625099u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9359789014u,       -0.1399625099u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8586989733u,        0.0139072118u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8586989733u,        0.0139072118u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9266016719u,       -0.0483306574u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9266016719u,       -0.0483306574u) withcolor MidnightBlue;
draw (        0.9359789014u,       -0.1399625099u)--(        0.9266019721u,       -0.0483306267u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9906224614u,        0.0178940348u)--(        0.9266019721u,       -0.0483306267u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9359789014u,       -0.1399625099u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9359789014u,       -0.1399625099u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9906224614u,        0.0178940348u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9906224614u,        0.0178940348u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9266019721u,       -0.0483306267u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9266019721u,       -0.0483306267u) withcolor MidnightBlue;
draw (        0.9359789014u,       -0.1399625099u)--(        0.9999993907u,       -0.0737378485u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9906224614u,        0.0178940348u)--(        0.9999993907u,       -0.0737378485u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9359789014u,       -0.1399625099u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9359789014u,       -0.1399625099u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9906224614u,        0.0178940348u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9906224614u,        0.0178940348u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9999993907u,       -0.0737378485u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9999993907u,       -0.0737378485u) withcolor Green;
draw (        0.8586989733u,        0.0139072118u)--(        0.9266018278u,       -0.0483304872u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9906224614u,        0.0178940348u)--(        0.9266018278u,       -0.0483304872u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8586989733u,        0.0139072118u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8586989733u,        0.0139072118u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9906224614u,        0.0178940348u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9906224614u,        0.0178940348u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9266018278u,       -0.0483304872u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9266018278u,       -0.0483304872u) withcolor MidnightBlue;
draw (        0.8586989733u,        0.0139072118u)--(        0.8102037820u,       -0.0644033574u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7215833170u,       -0.0392887795u)--(        0.8102037820u,       -0.0644033574u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8586989733u,        0.0139072118u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8586989733u,        0.0139072118u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7215833170u,       -0.0392887795u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7215833170u,       -0.0392887795u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8102037820u,       -0.0644033574u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8102037820u,       -0.0644033574u) withcolor MidnightBlue;
draw (        0.8586989733u,        0.0139072118u)--(        0.8664227061u,        0.1056932281u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9313924961u,        0.1709868359u)--(        0.8664227061u,        0.1056932281u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8586989733u,        0.0139072118u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8586989733u,        0.0139072118u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9313924961u,        0.1709868359u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9313924961u,        0.1709868359u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8664227061u,        0.1056932281u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8664227061u,        0.1056932281u) withcolor MidnightBlue;
draw (        0.8586989733u,        0.0139072118u)--(        0.8663853692u,        0.1056963623u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7742760417u,        0.1061442564u)--(        0.8663853692u,        0.1056963623u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8586989733u,        0.0139072118u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8586989733u,        0.0139072118u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7742760417u,        0.1061442564u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7742760417u,        0.1061442564u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8663853692u,        0.1056963623u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8663853692u,        0.1056963623u) withcolor MidnightBlue;
draw (        0.9906224614u,        0.0178940348u)--(        0.9999992466u,        0.1095259328u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9313924961u,        0.1709868359u)--(        0.9999992466u,        0.1095259328u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9906224614u,        0.0178940348u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9906224614u,        0.0178940348u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9313924961u,        0.1709868359u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9313924961u,        0.1709868359u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9999992466u,        0.1095259328u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9999992466u,        0.1095259328u) withcolor Green;
draw (        0.7215833170u,       -0.0392887795u)--(        0.6299482001u,       -0.0299435031u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6092389919u,        0.0598086975u)--(        0.6299482001u,       -0.0299435031u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7215833170u,       -0.0392887795u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7215833170u,       -0.0392887795u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6092389919u,        0.0598086975u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6092389919u,        0.0598086975u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6299482001u,       -0.0299435031u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6299482001u,       -0.0299435031u) withcolor Green;
draw (        0.7215833170u,       -0.0392887795u)--(        0.7008741088u,        0.0504634212u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6092389919u,        0.0598086975u)--(        0.7008741088u,        0.0504634212u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7215833170u,       -0.0392887795u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7215833170u,       -0.0392887795u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6092389919u,        0.0598086975u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6092389919u,        0.0598086975u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7008741088u,        0.0504634212u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7008741088u,        0.0504634212u) withcolor MidnightBlue;
draw (        0.7215833170u,       -0.0392887795u)--(        0.7008965506u,        0.0504685964u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7742760417u,        0.1061442564u)--(        0.7008965506u,        0.0504685964u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7215833170u,       -0.0392887795u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7215833170u,       -0.0392887795u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7742760417u,        0.1061442564u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7742760417u,        0.1061442564u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7008965506u,        0.0504685964u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7008965506u,        0.0504685964u) withcolor MidnightBlue;
draw (        0.9313924961u,        0.1709868359u)--(        0.9999620963u,        0.2324891835u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9743991409u,        0.3209813535u)--(        0.9999620963u,        0.2324891835u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9313924961u,        0.1709868359u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9313924961u,        0.1709868359u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9743991409u,        0.3209813535u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9743991409u,        0.3209813535u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9999620963u,        0.2324891835u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9999620963u,        0.2324891835u) withcolor Green;
draw (        0.9313924961u,        0.1709868359u)--(        0.8663855278u,        0.1057302431u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7742760417u,        0.1061442564u)--(        0.8663855278u,        0.1057302431u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9313924961u,        0.1709868359u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9313924961u,        0.1709868359u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7742760417u,        0.1061442564u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7742760417u,        0.1061442564u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8663855278u,        0.1057302431u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8663855278u,        0.1057302431u) withcolor MidnightBlue;
draw (        0.9313924961u,        0.1709868359u)--(        0.8392830100u,        0.1714008492u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7742760417u,        0.1061442564u)--(        0.8392830100u,        0.1714008492u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9313924961u,        0.1709868359u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9313924961u,        0.1709868359u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7742760417u,        0.1061442564u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7742760417u,        0.1061442564u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8392830100u,        0.1714008492u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8392830100u,        0.1714008492u) withcolor MidnightBlue;
draw (        0.9313924961u,        0.1709868359u)--(        0.8392830377u,        0.1714069612u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8416070986u,        0.2634880536u)--(        0.8392830377u,        0.1714069612u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9313924961u,        0.1709868359u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9313924961u,        0.1709868359u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8416070986u,        0.2634880536u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8416070986u,        0.2634880536u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8392830377u,        0.1714069612u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8392830377u,        0.1714069612u) withcolor MidnightBlue;
draw (        0.9743991409u,        0.3209813535u)--(        0.8853542139u,        0.3445467900u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8416070986u,        0.2634880536u)--(        0.8853542139u,        0.3445467900u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9743991409u,        0.3209813535u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9743991409u,        0.3209813535u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8416070986u,        0.2634880536u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8416070986u,        0.2634880536u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8853542139u,        0.3445467900u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8853542139u,        0.3445467900u) withcolor MidnightBlue;
draw (        0.9743991409u,        0.3209813535u)--(        0.8853578928u,        0.3445606870u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8520557521u,        0.4304402320u)--(        0.8853578928u,        0.3445606870u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9743991409u,        0.3209813535u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9743991409u,        0.3209813535u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8520557521u,        0.4304402320u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8520557521u,        0.4304402320u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8853578928u,        0.3445606870u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8853578928u,        0.3445606870u) withcolor MidnightBlue;
draw (        0.9743991409u,        0.3209813535u)--(        0.9410970002u,        0.4068608985u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8520557521u,        0.4304402320u)--(        0.9410970002u,        0.4068608985u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9743991409u,        0.3209813535u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9743991409u,        0.3209813535u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8520557521u,        0.4304402320u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8520557521u,        0.4304402320u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9410970002u,        0.4068608985u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9410970002u,        0.4068608985u) withcolor MidnightBlue;
draw (        0.9743991409u,        0.3209813535u)--(        0.9411072808u,        0.4068648844u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9666841951u,        0.4953530208u)--(        0.9411072808u,        0.4068648844u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9743991409u,        0.3209813535u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9743991409u,        0.3209813535u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9666841951u,        0.4953530208u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9666841951u,        0.4953530208u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9411072808u,        0.4068648844u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9411072808u,        0.4068648844u) withcolor MidnightBlue;
draw (        0.0666796245u,       -0.2697584450u)--(        0.0000000000u,       -0.3333333333u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0666796245u,       -0.2697584450u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0666796245u,       -0.2697584450u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,       -0.3333333333u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,       -0.3333333333u) withcolor Yellow;
draw (        0.0666796245u,       -0.2697584450u)--(        0.1333333333u,       -0.3333333333u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0666796245u,       -0.2697584450u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0666796245u,       -0.2697584450u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1333333333u,       -0.3333333333u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1333333333u,       -0.3333333333u) withcolor Yellow;
draw (        0.0666796245u,       -0.2697584450u)--(        0.0000158342u,       -0.2061953434u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0577516326u,       -0.1344255279u)--(        0.0000158342u,       -0.2061953434u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0666796245u,       -0.2697584450u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0666796245u,       -0.2697584450u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0577516326u,       -0.1344255279u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0577516326u,       -0.1344255279u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000158342u,       -0.2061953434u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000158342u,       -0.2061953434u) withcolor Green;
draw (        0.2020772361u,       -0.1652950001u)--(        0.1414439501u,       -0.0959557366u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0577516326u,       -0.1344255279u)--(        0.1414439501u,       -0.0959557366u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2020772361u,       -0.1652950001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2020772361u,       -0.1652950001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0577516326u,       -0.1344255279u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0577516326u,       -0.1344255279u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1414439501u,       -0.0959557366u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1414439501u,       -0.0959557366u) withcolor MidnightBlue;
draw (        0.2020772361u,       -0.1652950001u)--(        0.2229985628u,       -0.0755920073u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3023505295u,       -0.0288202827u)--(        0.2229985628u,       -0.0755920073u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2020772361u,       -0.1652950001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2020772361u,       -0.1652950001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3023505295u,       -0.0288202827u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3023505295u,       -0.0288202827u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2229985628u,       -0.0755920073u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2229985628u,       -0.0755920073u) withcolor MidnightBlue;
draw (        0.2020772361u,       -0.1652950001u)--(        0.2814292027u,       -0.1185232755u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3023505295u,       -0.0288202827u)--(        0.2814292027u,       -0.1185232755u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2020772361u,       -0.1652950001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2020772361u,       -0.1652950001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3023505295u,       -0.0288202827u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3023505295u,       -0.0288202827u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2814292027u,       -0.1185232755u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2814292027u,       -0.1185232755u) withcolor Green;
draw (        0.2020772361u,       -0.1652950001u)--(        0.1414558749u,       -0.0959453108u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1623713862u,       -0.0062409618u)--(        0.1414558749u,       -0.0959453108u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2020772361u,       -0.1652950001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2020772361u,       -0.1652950001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1623713862u,       -0.0062409618u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1623713862u,       -0.0062409618u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1414558749u,       -0.0959453108u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1414558749u,       -0.0959453108u) withcolor MidnightBlue;
draw (        0.2020772361u,       -0.1652950001u)--(        0.2229927473u,       -0.0755906512u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1623713862u,       -0.0062409618u)--(        0.2229927473u,       -0.0755906512u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2020772361u,       -0.1652950001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2020772361u,       -0.1652950001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1623713862u,       -0.0062409618u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1623713862u,       -0.0062409618u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2229927473u,       -0.0755906512u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2229927473u,       -0.0755906512u) withcolor MidnightBlue;
draw (        0.4574282374u,        0.0424733375u)--(        0.5274959226u,        0.1022633728u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6092389919u,        0.0598086975u)--(        0.5274959226u,        0.1022633728u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4574282374u,        0.0424733375u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4574282374u,        0.0424733375u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6092389919u,        0.0598086975u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6092389919u,        0.0598086975u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5274959226u,        0.1022633728u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5274959226u,        0.1022633728u) withcolor MidnightBlue;
draw (        0.4574282374u,        0.0424733375u)--(        0.5391713067u,        0.0000186622u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6092389919u,        0.0598086975u)--(        0.5391713067u,        0.0000186622u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4574282374u,        0.0424733375u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4574282374u,        0.0424733375u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6092389919u,        0.0598086975u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6092389919u,        0.0598086975u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5391713067u,        0.0000186622u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5391713067u,        0.0000186622u) withcolor Green;
draw (        0.4574282374u,        0.0424733375u)--(        0.3654115689u,        0.0383186340u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3023505295u,       -0.0288202827u)--(        0.3654115689u,        0.0383186340u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4574282374u,        0.0424733375u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4574282374u,        0.0424733375u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3023505295u,       -0.0288202827u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3023505295u,       -0.0288202827u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3654115689u,        0.0383186340u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3654115689u,        0.0383186340u) withcolor MidnightBlue;
draw (        0.4574282374u,        0.0424733375u)--(        0.3654101435u,        0.0383503259u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3435654252u,        0.1278329339u)--(        0.3654101435u,        0.0383503259u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4574282374u,        0.0424733375u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4574282374u,        0.0424733375u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3435654252u,        0.1278329339u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3435654252u,        0.1278329339u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3654101435u,        0.0383503259u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3654101435u,        0.0383503259u) withcolor MidnightBlue;
draw (        0.4574282374u,        0.0424733375u)--(        0.5274848019u,        0.1022764027u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4920438508u,        0.1872956230u)--(        0.5274848019u,        0.1022764027u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4574282374u,        0.0424733375u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4574282374u,        0.0424733375u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4920438508u,        0.1872956230u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4920438508u,        0.1872956230u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5274848019u,        0.1022764027u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5274848019u,        0.1022764027u) withcolor MidnightBlue;
draw (        0.6092389919u,        0.0598086975u)--(        0.6826378042u,        0.1154588833u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7742760417u,        0.1061442564u)--(        0.6826378042u,        0.1154588833u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6092389919u,        0.0598086975u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6092389919u,        0.0598086975u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7742760417u,        0.1061442564u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7742760417u,        0.1061442564u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6826378042u,        0.1154588833u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6826378042u,        0.1154588833u) withcolor MidnightBlue;
draw (        0.6092389919u,        0.0598086975u)--(        0.7008772293u,        0.0504940706u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7742760417u,        0.1061442564u)--(        0.7008772293u,        0.0504940706u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6092389919u,        0.0598086975u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6092389919u,        0.0598086975u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7742760417u,        0.1061442564u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7742760417u,        0.1061442564u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7008772293u,        0.0504940706u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7008772293u,        0.0504940706u) withcolor MidnightBlue;
draw (        0.6092389919u,        0.0598086975u)--(        0.5275076403u,        0.1022859266u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4920438508u,        0.1872956230u)--(        0.5275076403u,        0.1022859266u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6092389919u,        0.0598086975u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6092389919u,        0.0598086975u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4920438508u,        0.1872956230u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4920438508u,        0.1872956230u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5275076403u,        0.1022859266u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5275076403u,        0.1022859266u) withcolor MidnightBlue;
draw (        0.6092389919u,        0.0598086975u)--(        0.5737752024u,        0.1448183939u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4920438508u,        0.1872956230u)--(        0.5737752024u,        0.1448183939u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6092389919u,        0.0598086975u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6092389919u,        0.0598086975u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4920438508u,        0.1872956230u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4920438508u,        0.1872956230u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5737752024u,        0.1448183939u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5737752024u,        0.1448183939u) withcolor MidnightBlue;
draw (        0.6092389919u,        0.0598086975u)--(        0.5737780918u,        0.1448195992u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6471686451u,        0.2004806763u)--(        0.5737780918u,        0.1448195992u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6092389919u,        0.0598086975u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6092389919u,        0.0598086975u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6471686451u,        0.2004806763u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6471686451u,        0.2004806763u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5737780918u,        0.1448195992u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5737780918u,        0.1448195992u) withcolor MidnightBlue;
draw (        0.6092389919u,        0.0598086975u)--(        0.6826295453u,        0.1154697747u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6471686451u,        0.2004806763u)--(        0.6826295453u,        0.1154697747u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6092389919u,        0.0598086975u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6092389919u,        0.0598086975u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6471686451u,        0.2004806763u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6471686451u,        0.2004806763u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6826295453u,        0.1154697747u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6826295453u,        0.1154697747u) withcolor MidnightBlue;
draw (        0.7742760417u,        0.1061442564u)--(        0.7766064268u,        0.1982251890u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8416070986u,        0.2634880536u)--(        0.7766064268u,        0.1982251890u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7742760417u,        0.1061442564u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7742760417u,        0.1061442564u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8416070986u,        0.2634880536u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8416070986u,        0.2634880536u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7766064268u,        0.1982251890u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7766064268u,        0.1982251890u) withcolor Green;
draw (        0.7742760417u,        0.1061442564u)--(        0.8392767135u,        0.1714071210u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8416070986u,        0.2634880536u)--(        0.8392767135u,        0.1714071210u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7742760417u,        0.1061442564u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7742760417u,        0.1061442564u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8416070986u,        0.2634880536u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8416070986u,        0.2634880536u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8392767135u,        0.1714071210u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8392767135u,        0.1714071210u) withcolor MidnightBlue;
draw (        0.7742760417u,        0.1061442564u)--(        0.6826393274u,        0.1154738558u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6471686451u,        0.2004806763u)--(        0.6826393274u,        0.1154738558u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7742760417u,        0.1061442564u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7742760417u,        0.1061442564u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6471686451u,        0.2004806763u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6471686451u,        0.2004806763u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6826393274u,        0.1154738558u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6826393274u,        0.1154738558u) withcolor MidnightBlue;
draw (        0.8416070986u,        0.2634880536u)--(        0.8853407646u,        0.3445540470u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8520557521u,        0.4304402320u)--(        0.8853407646u,        0.3445540470u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8416070986u,        0.2634880536u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8416070986u,        0.2634880536u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8520557521u,        0.4304402320u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8520557521u,        0.4304402320u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8853407646u,        0.3445540470u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8853407646u,        0.3445540470u) withcolor MidnightBlue;
draw (        0.0577516326u,       -0.1344255279u)--(        0.0786858706u,       -0.0447255473u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1623713862u,       -0.0062409618u)--(        0.0786858706u,       -0.0447255473u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0577516326u,       -0.1344255279u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0577516326u,       -0.1344255279u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1623713862u,       -0.0062409618u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1623713862u,       -0.0062409618u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0786858706u,       -0.0447255473u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0786858706u,       -0.0447255473u) withcolor MidnightBlue;
draw (        0.0577516326u,       -0.1344255279u)--(        0.1414371482u,       -0.0959409424u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1623713862u,       -0.0062409618u)--(        0.1414371482u,       -0.0959409424u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0577516326u,       -0.1344255279u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0577516326u,       -0.1344255279u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1623713862u,       -0.0062409618u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1623713862u,       -0.0062409618u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1414371482u,       -0.0959409424u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1414371482u,       -0.0959409424u) withcolor MidnightBlue;
draw (        0.0577516326u,       -0.1344255279u)--(        0.0000040085u,       -0.0626652273u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0209332224u,        0.0270359257u)--(        0.0000040085u,       -0.0626652273u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0577516326u,       -0.1344255279u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0577516326u,       -0.1344255279u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0209332224u,        0.0270359257u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0209332224u,        0.0270359257u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000040085u,       -0.0626652273u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000040085u,       -0.0626652273u) withcolor Green;
draw (        0.0577516326u,       -0.1344255279u)--(        0.0786808466u,       -0.0447243749u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0209332224u,        0.0270359257u)--(        0.0786808466u,       -0.0447243749u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0577516326u,       -0.1344255279u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0577516326u,       -0.1344255279u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0209332224u,        0.0270359257u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0209332224u,        0.0270359257u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0786808466u,       -0.0447243749u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0786808466u,       -0.0447243749u) withcolor MidnightBlue;
draw (        0.3023505295u,       -0.0288202827u)--(        0.3653844822u,        0.0383440653u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3435654252u,        0.1278329339u)--(        0.3653844822u,        0.0383440653u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3023505295u,       -0.0288202827u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3023505295u,       -0.0288202827u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3435654252u,        0.1278329339u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3435654252u,        0.1278329339u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3653844822u,        0.0383440653u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3653844822u,        0.0383440653u) withcolor MidnightBlue;
draw (        0.3023505295u,       -0.0288202827u)--(        0.2229960579u,       -0.0755877572u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1623713862u,       -0.0062409618u)--(        0.2229960579u,       -0.0755877572u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3023505295u,       -0.0288202827u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3023505295u,       -0.0288202827u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1623713862u,       -0.0062409618u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1623713862u,       -0.0062409618u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2229960579u,       -0.0755877572u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2229960579u,       -0.0755877572u) withcolor MidnightBlue;
draw (        0.3023505295u,       -0.0288202827u)--(        0.2417258578u,        0.0405265127u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1623713862u,       -0.0062409618u)--(        0.2417258578u,        0.0405265127u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3023505295u,       -0.0288202827u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3023505295u,       -0.0288202827u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1623713862u,       -0.0062409618u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1623713862u,       -0.0062409618u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2417258578u,        0.0405265127u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2417258578u,        0.0405265127u) withcolor MidnightBlue;
draw (        0.3023505295u,       -0.0288202827u)--(        0.2417398172u,        0.0405387138u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2251209601u,        0.1311375131u)--(        0.2417398172u,        0.0405387138u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3023505295u,       -0.0288202827u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3023505295u,       -0.0288202827u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2251209601u,        0.1311375131u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2251209601u,        0.1311375131u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2417398172u,        0.0405387138u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2417398172u,        0.0405387138u) withcolor MidnightBlue;
draw (        0.3435654252u,        0.1278329339u)--(        0.4008127888u,        0.1999929519u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4920438508u,        0.1872956230u)--(        0.4008127888u,        0.1999929519u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3435654252u,        0.1278329339u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3435654252u,        0.1278329339u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4920438508u,        0.1872956230u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4920438508u,        0.1872956230u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4008127888u,        0.1999929519u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4008127888u,        0.1999929519u) withcolor Green;
draw (        0.3435654252u,        0.1278329339u)--(        0.2863101709u,        0.1999866911u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2251209601u,        0.1311375131u)--(        0.2863101709u,        0.1999866911u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3435654252u,        0.1278329339u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3435654252u,        0.1278329339u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2251209601u,        0.1311375131u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2251209601u,        0.1311375131u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2863101709u,        0.1999866911u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2863101709u,        0.1999866911u) withcolor Green;
draw (        0.4920438508u,        0.1872956230u)--(        0.5737767490u,        0.1448213697u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6471686451u,        0.2004806763u)--(        0.5737767490u,        0.1448213697u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4920438508u,        0.1872956230u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4920438508u,        0.1872956230u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6471686451u,        0.2004806763u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6471686451u,        0.2004806763u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5737767490u,        0.1448213697u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5737767490u,        0.1448213697u) withcolor MidnightBlue;
draw (        0.4920438508u,        0.1872956230u)--(        0.4711589093u,        0.2770070941u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5456882242u,        0.3311338903u)--(        0.4711589093u,        0.2770070941u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4920438508u,        0.1872956230u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4920438508u,        0.1872956230u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5456882242u,        0.3311338903u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5456882242u,        0.3311338903u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4711589093u,        0.2770070941u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4711589093u,        0.2770070941u) withcolor MidnightBlue;
draw (        0.4920438508u,        0.1872956230u)--(        0.4711355056u,        0.2770016425u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3916001728u,        0.3234608684u)--(        0.4711355056u,        0.2770016425u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4920438508u,        0.1872956230u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4920438508u,        0.1872956230u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3916001728u,        0.3234608684u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3916001728u,        0.3234608684u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4711355056u,        0.2770016425u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4711355056u,        0.2770016425u) withcolor MidnightBlue;
draw (        0.6471686451u,        0.2004806763u)--(        0.6284326045u,        0.2906654290u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5456882242u,        0.3311338903u)--(        0.6284326045u,        0.2906654290u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6471686451u,        0.2004806763u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6471686451u,        0.2004806763u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5456882242u,        0.3311338903u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5456882242u,        0.3311338903u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6284326045u,        0.2906654290u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6284326045u,        0.2906654290u) withcolor Green;
draw (        0.8520557521u,        0.4304402320u)--(        0.9410987115u,        0.4068673618u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9666841951u,        0.4953530208u)--(        0.9410987115u,        0.4068673618u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8520557521u,        0.4304402320u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8520557521u,        0.4304402320u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9666841951u,        0.4953530208u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9666841951u,        0.4953530208u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9410987115u,        0.4068673618u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9410987115u,        0.4068673618u) withcolor MidnightBlue;
draw (        0.8520557521u,        0.4304402320u)--(        0.7718018320u,        0.4756468373u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7445701922u,        0.5636398292u)--(        0.7718018320u,        0.4756468373u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8520557521u,        0.4304402320u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8520557521u,        0.4304402320u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7445701922u,        0.5636398292u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7445701922u,        0.5636398292u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7718018320u,        0.4756468373u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7718018320u,        0.4756468373u) withcolor Green;
draw (        0.8520557521u,        0.4304402320u)--(        0.8248241122u,        0.5184332239u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7445701922u,        0.5636398292u)--(        0.8248241122u,        0.5184332239u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8520557521u,        0.4304402320u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8520557521u,        0.4304402320u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7445701922u,        0.5636398292u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7445701922u,        0.5636398292u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8248241122u,        0.5184332239u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8248241122u,        0.5184332239u) withcolor MidnightBlue;
draw (        0.8520557521u,        0.4304402320u)--(        0.8248379900u,        0.5184375176u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8576592976u,        0.6045019736u)--(        0.8248379900u,        0.5184375176u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8520557521u,        0.4304402320u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8520557521u,        0.4304402320u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8576592976u,        0.6045019736u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8576592976u,        0.6045019736u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8248379900u,        0.5184375176u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8248379900u,        0.5184375176u) withcolor MidnightBlue;
draw (        0.1623713862u,       -0.0062409618u)--(        0.1457716176u,        0.0843613368u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2251209601u,        0.1311375131u)--(        0.1457716176u,        0.0843613368u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1623713862u,       -0.0062409618u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1623713862u,       -0.0062409618u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2251209601u,        0.1311375131u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2251209601u,        0.1311375131u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1457716176u,        0.0843613368u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1457716176u,        0.0843613368u) withcolor MidnightBlue;
draw (        0.1623713862u,       -0.0062409618u)--(        0.2417207288u,        0.0405352145u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2251209601u,        0.1311375131u)--(        0.2417207288u,        0.0405352145u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1623713862u,       -0.0062409618u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1623713862u,       -0.0062409618u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2251209601u,        0.1311375131u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2251209601u,        0.1311375131u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2417207288u,        0.0405352145u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2417207288u,        0.0405352145u) withcolor MidnightBlue;
draw (        0.1623713862u,       -0.0062409618u)--(        0.0786841201u,       -0.0447217405u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0209332224u,        0.0270359257u)--(        0.0786841201u,       -0.0447217405u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1623713862u,       -0.0062409618u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1623713862u,       -0.0062409618u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0209332224u,        0.0270359257u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0209332224u,        0.0270359257u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0786841201u,       -0.0447217405u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0786841201u,       -0.0447217405u) withcolor MidnightBlue;
draw (        0.1623713862u,       -0.0062409618u)--(        0.1457529591u,        0.0843579163u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0845830960u,        0.1532242846u)--(        0.1457529591u,        0.0843579163u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1623713862u,       -0.0062409618u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1623713862u,       -0.0062409618u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0845830960u,        0.1532242846u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0845830960u,        0.1532242846u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1457529591u,        0.0843579163u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1457529591u,        0.0843579163u) withcolor MidnightBlue;
draw (        0.2251209601u,        0.1311375131u)--(        0.1457665282u,        0.0843699713u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0845830960u,        0.1532242846u)--(        0.1457665282u,        0.0843699713u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2251209601u,        0.1311375131u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2251209601u,        0.1311375131u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0845830960u,        0.1532242846u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0845830960u,        0.1532242846u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1457665282u,        0.0843699713u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1457665282u,        0.0843699713u) withcolor MidnightBlue;
draw (        0.2251209601u,        0.1311375131u)--(        0.1639375279u,        0.1999918264u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0845830960u,        0.1532242846u)--(        0.1639375279u,        0.1999918264u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2251209601u,        0.1311375131u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2251209601u,        0.1311375131u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0845830960u,        0.1532242846u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0845830960u,        0.1532242846u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1639375279u,        0.1999918264u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1639375279u,        0.1999918264u) withcolor Green;
draw (        0.5456882242u,        0.3311338903u)--(        0.4711477023u,        0.2770225289u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3916001728u,        0.3234608684u)--(        0.4711477023u,        0.2770225289u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5456882242u,        0.3311338903u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5456882242u,        0.3311338903u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3916001728u,        0.3234608684u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3916001728u,        0.3234608684u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4711477023u,        0.2770225289u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4711477023u,        0.2770225289u) withcolor MidnightBlue;
draw (        0.5456882242u,        0.3311338903u)--(        0.4661406946u,        0.3775722298u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3916001728u,        0.3234608684u)--(        0.4661406946u,        0.3775722298u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5456882242u,        0.3311338903u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5456882242u,        0.3311338903u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3916001728u,        0.3234608684u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3916001728u,        0.3234608684u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4661406946u,        0.3775722298u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4661406946u,        0.3775722298u) withcolor MidnightBlue;
draw (        0.5456882242u,        0.3311338903u)--(        0.4661469494u,        0.3775829424u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4736106869u,        0.4693904663u)--(        0.4661469494u,        0.3775829424u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5456882242u,        0.3311338903u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5456882242u,        0.3311338903u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4736106869u,        0.4693904663u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4736106869u,        0.4693904663u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4661469494u,        0.3775829424u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4661469494u,        0.3775829424u) withcolor MidnightBlue;
draw (        0.5456882242u,        0.3311338903u)--(        0.5531519616u,        0.4229414143u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4736106869u,        0.4693904663u)--(        0.5531519616u,        0.4229414143u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5456882242u,        0.3311338903u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5456882242u,        0.3311338903u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4736106869u,        0.4693904663u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4736106869u,        0.4693904663u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5531519616u,        0.4229414143u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5531519616u,        0.4229414143u) withcolor MidnightBlue;
draw (        0.5456882242u,        0.3311338903u)--(        0.5531543681u,        0.4229412186u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6346281184u,        0.4659104732u)--(        0.5531543681u,        0.4229412186u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5456882242u,        0.3311338903u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5456882242u,        0.3311338903u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6346281184u,        0.4659104732u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6346281184u,        0.4659104732u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5531543681u,        0.4229412186u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5531543681u,        0.4229412186u) withcolor MidnightBlue;
draw (        0.5456882242u,        0.3311338903u)--(        0.6271619745u,        0.3741031449u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6346281184u,        0.4659104732u)--(        0.6271619745u,        0.3741031449u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5456882242u,        0.3311338903u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5456882242u,        0.3311338903u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6346281184u,        0.4659104732u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6346281184u,        0.4659104732u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6271619745u,        0.3741031449u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6271619745u,        0.3741031449u) withcolor Green;
draw (        0.0209332224u,        0.0270359257u)--(        0.0000072817u,        0.1167378422u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0845830960u,        0.1532242846u)--(        0.0000072817u,        0.1167378422u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0209332224u,        0.0270359257u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0209332224u,        0.0270359257u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0845830960u,        0.1532242846u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0845830960u,        0.1532242846u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000072817u,        0.1167378422u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000072817u,        0.1167378422u) withcolor Green;
draw (        0.0845830960u,        0.1532242846u)--(        0.0000056990u,        0.1897070580u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0293107762u,        0.2770314033u)--(        0.0000056990u,        0.1897070580u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0845830960u,        0.1532242846u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0845830960u,        0.1532242846u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0293107762u,        0.2770314033u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0293107762u,        0.2770314033u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000056990u,        0.1897070580u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000056990u,        0.1897070580u) withcolor Green;
draw (        0.0845830960u,        0.1532242846u)--(        0.1481960066u,        0.2198405471u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1459178650u,        0.3119227870u)--(        0.1481960066u,        0.2198405471u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0845830960u,        0.1532242846u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0845830960u,        0.1532242846u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1459178650u,        0.3119227870u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1459178650u,        0.3119227870u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1481960066u,        0.2198405471u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1481960066u,        0.2198405471u) withcolor Green;
draw (        0.3916001728u,        0.3234608684u)--(        0.3000000000u,        0.3333333333u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3916001728u,        0.3234608684u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3916001728u,        0.3234608684u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3000000000u,        0.3333333333u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3000000000u,        0.3333333333u) withcolor Yellow;
draw (        0.3916001728u,        0.3234608684u)--(        0.4661320355u,        0.3775841560u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4736106869u,        0.4693904663u)--(        0.4661320355u,        0.3775841560u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3916001728u,        0.3234608684u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3916001728u,        0.3234608684u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4736106869u,        0.4693904663u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4736106869u,        0.4693904663u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4661320355u,        0.3775841560u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4661320355u,        0.3775841560u) withcolor MidnightBlue;
draw (        0.4736106869u,        0.4693904663u)--(        0.5531531719u,        0.4229434868u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6346281184u,        0.4659104732u)--(        0.5531531719u,        0.4229434868u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4736106869u,        0.4693904663u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4736106869u,        0.4693904663u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6346281184u,        0.4659104732u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6346281184u,        0.4659104732u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5531531719u,        0.4229434868u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5531531719u,        0.4229434868u) withcolor MidnightBlue;
draw (        0.4736106869u,        0.4693904663u)--(        0.4020899258u,        0.5274344977u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3100245403u,        0.5245546333u)--(        0.4020899258u,        0.5274344977u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4736106869u,        0.4693904663u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4736106869u,        0.4693904663u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3100245403u,        0.5245546333u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3100245403u,        0.5245546333u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4020899258u,        0.5274344977u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4020899258u,        0.5274344977u) withcolor MidnightBlue;
draw (        0.4736106869u,        0.4693904663u)--(        0.4510430535u,        0.5586934940u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5283593199u,        0.6087586894u)--(        0.4510430535u,        0.5586934940u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4736106869u,        0.4693904663u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4736106869u,        0.4693904663u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5283593199u,        0.6087586894u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5283593199u,        0.6087586894u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4510430535u,        0.5586934940u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4510430535u,        0.5586934940u) withcolor MidnightBlue;
draw (        0.4736106869u,        0.4693904663u)--(        0.5509269533u,        0.5194556616u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5283593199u,        0.6087586894u)--(        0.5509269533u,        0.5194556616u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4736106869u,        0.4693904663u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4736106869u,        0.4693904663u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5283593199u,        0.6087586894u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5283593199u,        0.6087586894u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5509269533u,        0.5194556616u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5509269533u,        0.5194556616u) withcolor Green;
draw (        0.4736106869u,        0.4693904663u)--(        0.4510287198u,        0.5586898706u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3794058650u,        0.6166078773u)--(        0.4510287198u,        0.5586898706u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4736106869u,        0.4693904663u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4736106869u,        0.4693904663u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3794058650u,        0.6166078773u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3794058650u,        0.6166078773u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4510287198u,        0.5586898706u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4510287198u,        0.5586898706u) withcolor MidnightBlue;
draw (        0.6346281184u,        0.4659104732u)--(        0.6527592068u,        0.5562187922u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7445701922u,        0.5636398292u)--(        0.6527592068u,        0.5562187922u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6346281184u,        0.4659104732u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6346281184u,        0.4659104732u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7445701922u,        0.5636398292u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7445701922u,        0.5636398292u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6527592068u,        0.5562187922u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6527592068u,        0.5562187922u) withcolor Green;
draw (        0.7445701922u,        0.5636398292u)--(        0.7774009328u,        0.6497006873u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8576592976u,        0.6045019736u)--(        0.7774009328u,        0.6497006873u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7445701922u,        0.5636398292u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7445701922u,        0.5636398292u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8576592976u,        0.6045019736u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8576592976u,        0.6045019736u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7774009328u,        0.6497006873u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7774009328u,        0.6497006873u) withcolor Green;
draw (        0.7445701922u,        0.5636398292u)--(        0.8248285570u,        0.5184411155u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8576592976u,        0.6045019736u)--(        0.8248285570u,        0.5184411155u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7445701922u,        0.5636398292u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7445701922u,        0.5636398292u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8576592976u,        0.6045019736u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8576592976u,        0.6045019736u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8248285570u,        0.5184411155u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8248285570u,        0.5184411155u) withcolor MidnightBlue;
draw (        0.0293107762u,        0.2770314033u)--(        0.0677936910u,        0.3607176871u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1459178650u,        0.3119227870u)--(        0.0677936910u,        0.3607176871u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0293107762u,        0.2770314033u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0293107762u,        0.2770314033u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1459178650u,        0.3119227870u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1459178650u,        0.3119227870u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0677936910u,        0.3607176871u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0677936910u,        0.3607176871u) withcolor MidnightBlue;
draw (        0.0293107762u,        0.2770314033u)--(        0.0000008899u,        0.3643541346u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0384789954u,        0.4480426299u)--(        0.0000008899u,        0.3643541346u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0293107762u,        0.2770314033u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0293107762u,        0.2770314033u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0384789954u,        0.4480426299u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0384789954u,        0.4480426299u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000008899u,        0.3643541346u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000008899u,        0.3643541346u) withcolor Green;
draw (        0.0293107762u,        0.2770314033u)--(        0.0677888816u,        0.3607198986u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0384789954u,        0.4480426299u)--(        0.0677888816u,        0.3607198986u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0293107762u,        0.2770314033u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0293107762u,        0.2770314033u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0384789954u,        0.4480426299u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0384789954u,        0.4480426299u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0677888816u,        0.3607198986u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0677888816u,        0.3607198986u) withcolor MidnightBlue;
draw (        0.1459178650u,        0.3119227870u)--(        0.0677967149u,        0.3607225282u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0384789954u,        0.4480426299u)--(        0.0677967149u,        0.3607225282u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1459178650u,        0.3119227870u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1459178650u,        0.3119227870u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0384789954u,        0.4480426299u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0384789954u,        0.4480426299u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0677967149u,        0.3607225282u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0677967149u,        0.3607225282u) withcolor MidnightBlue;
draw (        0.1459178650u,        0.3119227870u)--(        0.1166001455u,        0.3992428887u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0384789954u,        0.4480426299u)--(        0.1166001455u,        0.3992428887u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1459178650u,        0.3119227870u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1459178650u,        0.3119227870u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0384789954u,        0.4480426299u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0384789954u,        0.4480426299u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1166001455u,        0.3992428887u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1166001455u,        0.3992428887u) withcolor MidnightBlue;
draw (        0.1459178650u,        0.3119227870u)--(        0.1166040539u,        0.3992442008u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1585769498u,        0.4812356933u)--(        0.1166040539u,        0.3992442008u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1459178650u,        0.3119227870u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1459178650u,        0.3119227870u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1585769498u,        0.4812356933u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1585769498u,        0.4812356933u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1166040539u,        0.3992442008u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1166040539u,        0.3992442008u) withcolor MidnightBlue;
draw (        0.1459178650u,        0.3119227870u)--(        0.1878907609u,        0.3939142795u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1585769498u,        0.4812356933u)--(        0.1878907609u,        0.3939142795u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1459178650u,        0.3119227870u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1459178650u,        0.3119227870u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1585769498u,        0.4812356933u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1585769498u,        0.4812356933u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1878907609u,        0.3939142795u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1878907609u,        0.3939142795u) withcolor Green;
draw (        0.3100245403u,        0.5245546333u)--(        0.4020929805u,        0.5273351288u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3794058650u,        0.6166078773u)--(        0.4020929805u,        0.5273351288u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3100245403u,        0.5245546333u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3100245403u,        0.5245546333u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3794058650u,        0.6166078773u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3794058650u,        0.6166078773u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4020929805u,        0.5273351288u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4020929805u,        0.5273351288u) withcolor MidnightBlue;
draw (        0.5283593199u,        0.6087586894u)--(        0.6000000000u,        0.6666666667u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5283593199u,        0.6087586894u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5283593199u,        0.6087586894u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6000000000u,        0.6666666667u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6000000000u,        0.6666666667u) withcolor Yellow;
draw (        0.5283593199u,        0.6087586894u)--(        0.4510379822u,        0.5587013266u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3794058650u,        0.6166078773u)--(        0.4510379822u,        0.5587013266u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5283593199u,        0.6087586894u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5283593199u,        0.6087586894u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3794058650u,        0.6166078773u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3794058650u,        0.6166078773u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4510379822u,        0.5587013266u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4510379822u,        0.5587013266u) withcolor MidnightBlue;
draw (        0.5283593199u,        0.6087586894u)--(        0.4567272027u,        0.6666652401u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3794058650u,        0.6166078773u)--(        0.4567272027u,        0.6666652401u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5283593199u,        0.6087586894u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5283593199u,        0.6087586894u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3794058650u,        0.6166078773u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3794058650u,        0.6166078773u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4567272027u,        0.6666652401u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4567272027u,        0.6666652401u) withcolor Green;
draw (        0.0384789954u,        0.4480426299u)--(        0.1166017135u,        0.3992453990u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1585769498u,        0.4812356933u)--(        0.1166017135u,        0.3992453990u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0384789954u,        0.4480426299u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0384789954u,        0.4480426299u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1585769498u,        0.4812356933u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1585769498u,        0.4812356933u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1166017135u,        0.3992453990u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1166017135u,        0.3992453990u) withcolor MidnightBlue;
draw (        0.0384789954u,        0.4480426299u)--(        0.0000102910u,        0.5317354469u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0627186792u,        0.5992038608u)--(        0.0000102910u,        0.5317354469u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0384789954u,        0.4480426299u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0384789954u,        0.4480426299u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0627186792u,        0.5992038608u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0627186792u,        0.5992038608u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000102910u,        0.5317354469u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000102910u,        0.5317354469u) withcolor Green;
draw (        0.1585769498u,        0.4812356933u)--(        0.1510419095u,        0.5730373927u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2137582182u,        0.6404984441u)--(        0.1510419095u,        0.5730373927u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1585769498u,        0.4812356933u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1585769498u,        0.4812356933u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2137582182u,        0.6404984441u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2137582182u,        0.6404984441u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1510419095u,        0.5730373927u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1510419095u,        0.5730373927u) withcolor MidnightBlue;
draw (        0.1585769498u,        0.4812356933u)--(        0.1510340865u,        0.5730367502u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0627186792u,        0.5992038608u)--(        0.1510340865u,        0.5730367502u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1585769498u,        0.4812356933u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1585769498u,        0.4812356933u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0627186792u,        0.5992038608u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0627186792u,        0.5992038608u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1510340865u,        0.5730367502u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1510340865u,        0.5730367502u) withcolor MidnightBlue;
draw (        0.2137582182u,        0.6404984441u)--(        0.3020771121u,        0.6666537843u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3794058650u,        0.6166078773u)--(        0.3020771121u,        0.6666537843u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2137582182u,        0.6404984441u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2137582182u,        0.6404984441u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3794058650u,        0.6166078773u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3794058650u,        0.6166078773u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3020771121u,        0.6666537843u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3020771121u,        0.6666537843u) withcolor Green;
draw (        0.2137582182u,        0.6404984441u)--(        0.1510359251u,        0.5730429567u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0627186792u,        0.5992038608u)--(        0.1510359251u,        0.5730429567u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2137582182u,        0.6404984441u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2137582182u,        0.6404984441u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0627186792u,        0.5992038608u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0627186792u,        0.5992038608u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1510359251u,        0.5730429567u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1510359251u,        0.5730429567u) withcolor MidnightBlue;
draw (        0.2137582182u,        0.6404984441u)--(        0.1254409723u,        0.6666593482u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0627186792u,        0.5992038608u)--(        0.1254409723u,        0.6666593482u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2137582182u,        0.6404984441u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2137582182u,        0.6404984441u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0627186792u,        0.5992038608u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0627186792u,        0.5992038608u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1254409723u,        0.6666593482u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1254409723u,        0.6666593482u) withcolor Green;
draw (        0.0627186792u,        0.5992038608u)--(        0.0000000000u,        0.6666666667u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0627186792u,        0.5992038608u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0627186792u,        0.5992038608u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.6666666667u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.6666666667u) withcolor Yellow;
% r =         0.0921104165
endfig;
end;
