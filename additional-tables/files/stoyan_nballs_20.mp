prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(20);
u := 20cm;
path pols[];
pair p[];
R := .1327772115u;
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

x1 = .0666625429u;
y1 = -.2184996197u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .2805051586u;
y2 = -.0973106696u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .9308855346u;
y3 = -.2199544913u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .0821592729u;
y4 = .0006622743u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .5294258012u;
y5 = .0430088374u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .7332468348u;
y6 = -.0846842105u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .9250164787u;
y7 = .0030079435u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .0968666128u;
y8 = .1957863424u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .3117482354u;
y9 = .0956910691u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = .5195343715u;
y10 = .2430002972u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;

x11 = .7235486737u;
y11 = .1441060105u;
draw fullcircle scaled D shifted z11 withcolor RoyalBlue withpen pencircle scaled 1.0;

x12 = .9279892810u;
y12 = .2241396603u;
draw fullcircle scaled D shifted z12 withcolor RoyalBlue withpen pencircle scaled 1.0;

x13 = .1050889823u;
y13 = .3677977192u;
draw fullcircle scaled D shifted z13 withcolor RoyalBlue withpen pencircle scaled 1.0;

x14 = .3926632979u;
y14 = .3954647424u;
draw fullcircle scaled D shifted z14 withcolor RoyalBlue withpen pencircle scaled 1.0;

x15 = .6110373210u;
y15 = .4759455834u;
draw fullcircle scaled D shifted z15 withcolor RoyalBlue withpen pencircle scaled 1.0;

x16 = .9338902108u;
y16 = .4508519010u;
draw fullcircle scaled D shifted z16 withcolor RoyalBlue withpen pencircle scaled 1.0;

x17 = .0760653247u;
y17 = .5578371882u;
draw fullcircle scaled D shifted z17 withcolor RoyalBlue withpen pencircle scaled 1.0;

x18 = .2560161760u;
y18 = .5839228127u;
draw fullcircle scaled D shifted z18 withcolor RoyalBlue withpen pencircle scaled 1.0;

x19 = .4799239229u;
y19 = .6099583210u;
draw fullcircle scaled D shifted z19 withcolor RoyalBlue withpen pencircle scaled 1.0;

x20 = .8198184288u;
y20 = .5606356413u;
draw fullcircle scaled D shifted z20 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.0666625429u,       -0.2184996197u)--(        0.0000000000u,       -0.3333333333u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0666625429u,       -0.2184996197u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0666625429u,       -0.2184996197u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,       -0.3333333333u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,       -0.3333333333u) withcolor Yellow;
draw (        0.0666625429u,       -0.2184996197u)--(        0.1333333333u,       -0.3333333333u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0666625429u,       -0.2184996197u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0666625429u,       -0.2184996197u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1333333333u,       -0.3333333333u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1333333333u,       -0.3333333333u) withcolor Yellow;
draw (        0.0666625429u,       -0.2184996197u)--(        0.1488039320u,       -0.1141800142u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2805051586u,       -0.0973106696u)--(        0.1488039320u,       -0.1141800142u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0666625429u,       -0.2184996197u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0666625429u,       -0.2184996197u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2805051586u,       -0.0973106696u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2805051586u,       -0.0973106696u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1488039320u,       -0.1141800142u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1488039320u,       -0.1141800142u) withcolor MidnightBlue;
draw (        0.0666625429u,       -0.2184996197u)--(        0.1983637695u,       -0.2016302751u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2805051586u,       -0.0973106696u)--(        0.1983637695u,       -0.2016302751u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0666625429u,       -0.2184996197u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0666625429u,       -0.2184996197u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2805051586u,       -0.0973106696u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2805051586u,       -0.0973106696u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1983637695u,       -0.2016302751u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1983637695u,       -0.2016302751u) withcolor Green;
draw (        0.0666625429u,       -0.2184996197u)--(        0.0000193925u,       -0.1036585191u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0821592729u,        0.0006622743u)--(        0.0000193925u,       -0.1036585191u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0666625429u,       -0.2184996197u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0666625429u,       -0.2184996197u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0821592729u,        0.0006622743u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0821592729u,        0.0006622743u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000193925u,       -0.1036585191u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000193925u,       -0.1036585191u) withcolor Green;
draw (        0.0666625429u,       -0.2184996197u)--(        0.1488024233u,       -0.1141788263u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0821592729u,        0.0006622743u)--(        0.1488024233u,       -0.1141788263u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0666625429u,       -0.2184996197u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0666625429u,       -0.2184996197u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0821592729u,        0.0006622743u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0821592729u,        0.0006622743u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1488024233u,       -0.1141788263u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1488024233u,       -0.1141788263u) withcolor MidnightBlue;
draw (        0.2805051586u,       -0.0973106696u)--(        0.1488036860u,       -0.1141780936u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0821592729u,        0.0006622743u)--(        0.1488036860u,       -0.1141780936u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2805051586u,       -0.0973106696u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2805051586u,       -0.0973106696u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0821592729u,        0.0006622743u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0821592729u,        0.0006622743u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1488036860u,       -0.1141780936u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1488036860u,       -0.1141780936u) withcolor MidnightBlue;
draw (        0.2805051586u,       -0.0973106696u)--(        0.3848238096u,       -0.0151680683u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3117482354u,        0.0956910691u)--(        0.3848238096u,       -0.0151680683u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2805051586u,       -0.0973106696u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2805051586u,       -0.0973106696u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3117482354u,        0.0956910691u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3117482354u,        0.0956910691u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3848238096u,       -0.0151680683u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3848238096u,       -0.0151680683u) withcolor Green;
draw (        0.9308855346u,       -0.2199544913u)--(        1.0000000000u,       -0.3333333333u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9308855346u,       -0.2199544913u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9308855346u,       -0.2199544913u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,       -0.3333333333u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,       -0.3333333333u) withcolor Yellow;
draw (        0.9308855346u,       -0.2199544913u)--(        0.7996695473u,       -0.1996529503u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7332468348u,       -0.0846842105u)--(        0.7996695473u,       -0.1996529503u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9308855346u,       -0.2199544913u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9308855346u,       -0.2199544913u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7332468348u,       -0.0846842105u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7332468348u,       -0.0846842105u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7996695473u,       -0.1996529503u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7996695473u,       -0.1996529503u) withcolor Green;
draw (        0.9308855346u,       -0.2199544913u)--(        0.9999893211u,       -0.1065770043u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9250164787u,        0.0030079435u)--(        0.9999893211u,       -0.1065770043u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9308855346u,       -0.2199544913u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9308855346u,       -0.2199544913u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9250164787u,        0.0030079435u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9250164787u,        0.0030079435u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9999893211u,       -0.1065770043u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9999893211u,       -0.1065770043u) withcolor Green;
draw (        0.0821592729u,        0.0006622743u)--(        0.0000040417u,        0.1049709790u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0968666128u,        0.1957863424u)--(        0.0000040417u,        0.1049709790u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0821592729u,        0.0006622743u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0821592729u,        0.0006622743u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0968666128u,        0.1957863424u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0968666128u,        0.1957863424u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000040417u,        0.1049709790u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000040417u,        0.1049709790u) withcolor Green;
draw (        0.0821592729u,        0.0006622743u)--(        0.1790218440u,        0.0914776377u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0968666128u,        0.1957863424u)--(        0.1790218440u,        0.0914776377u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0821592729u,        0.0006622743u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0821592729u,        0.0006622743u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0968666128u,        0.1957863424u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0968666128u,        0.1957863424u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1790218440u,        0.0914776377u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1790218440u,        0.0914776377u) withcolor MidnightBlue;
draw (        0.0821592729u,        0.0006622743u)--(        0.1790384572u,        0.0914599149u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3117482354u,        0.0956910691u)--(        0.1790384572u,        0.0914599149u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0821592729u,        0.0006622743u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0821592729u,        0.0006622743u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3117482354u,        0.0956910691u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3117482354u,        0.0956910691u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1790384572u,        0.0914599149u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1790384572u,        0.0914599149u) withcolor MidnightBlue;
draw (        0.5294258012u,        0.0430088374u)--(        0.6612169714u,        0.0268572019u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7332468348u,       -0.0846842105u)--(        0.6612169714u,        0.0268572019u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5294258012u,        0.0430088374u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5294258012u,        0.0430088374u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7332468348u,       -0.0846842105u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7332468348u,       -0.0846842105u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6612169714u,        0.0268572019u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6612169714u,        0.0268572019u) withcolor MidnightBlue;
draw (        0.5294258012u,        0.0430088374u)--(        0.4038047317u,        0.0000072702u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3117482354u,        0.0956910691u)--(        0.4038047317u,        0.0000072702u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5294258012u,        0.0430088374u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5294258012u,        0.0430088374u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3117482354u,        0.0956910691u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3117482354u,        0.0956910691u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4038047317u,        0.0000072702u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4038047317u,        0.0000072702u) withcolor Green;
draw (        0.5294258012u,        0.0430088374u)--(        0.4373693049u,        0.1386926363u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3117482354u,        0.0956910691u)--(        0.4373693049u,        0.1386926363u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5294258012u,        0.0430088374u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5294258012u,        0.0430088374u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3117482354u,        0.0956910691u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3117482354u,        0.0956910691u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4373693049u,        0.1386926363u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4373693049u,        0.1386926363u) withcolor MidnightBlue;
draw (        0.5294258012u,        0.0430088374u)--(        0.4373731359u,        0.1386963220u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5195343715u,        0.2430002972u)--(        0.4373731359u,        0.1386963220u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5294258012u,        0.0430088374u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5294258012u,        0.0430088374u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5195343715u,        0.2430002972u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5195343715u,        0.2430002972u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4373731359u,        0.1386963220u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4373731359u,        0.1386963220u) withcolor MidnightBlue;
draw (        0.5294258012u,        0.0430088374u)--(        0.6612183137u,        0.0268681578u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7235486737u,        0.1441060105u)--(        0.6612183137u,        0.0268681578u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5294258012u,        0.0430088374u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5294258012u,        0.0430088374u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7235486737u,        0.1441060105u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7235486737u,        0.1441060105u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6612183137u,        0.0268681578u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6612183137u,        0.0268681578u) withcolor MidnightBlue;
draw (        0.7332468348u,       -0.0846842105u)--(        0.7955692329u,        0.0325578748u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9250164787u,        0.0030079435u)--(        0.7955692329u,        0.0325578748u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7332468348u,       -0.0846842105u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7332468348u,       -0.0846842105u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9250164787u,        0.0030079435u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9250164787u,        0.0030079435u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7955692329u,        0.0325578748u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7955692329u,        0.0325578748u) withcolor MidnightBlue;
draw (        0.7332468348u,       -0.0846842105u)--(        0.6612268834u,        0.0268636021u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7235486737u,        0.1441060105u)--(        0.6612268834u,        0.0268636021u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7332468348u,       -0.0846842105u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7332468348u,       -0.0846842105u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7235486737u,        0.1441060105u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7235486737u,        0.1441060105u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6612268834u,        0.0268636021u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6612268834u,        0.0268636021u) withcolor MidnightBlue;
draw (        0.7332468348u,       -0.0846842105u)--(        0.7955686251u,        0.0325581979u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7235486737u,        0.1441060105u)--(        0.7955686251u,        0.0325581979u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7332468348u,       -0.0846842105u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7332468348u,       -0.0846842105u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7235486737u,        0.1441060105u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7235486737u,        0.1441060105u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7955686251u,        0.0325581979u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7955686251u,        0.0325581979u) withcolor MidnightBlue;
draw (        0.9250164787u,        0.0030079435u)--(        0.7955694245u,        0.0325587140u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7235486737u,        0.1441060105u)--(        0.7955694245u,        0.0325587140u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9250164787u,        0.0030079435u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9250164787u,        0.0030079435u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7235486737u,        0.1441060105u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7235486737u,        0.1441060105u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7955694245u,        0.0325587140u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7955694245u,        0.0325587140u) withcolor MidnightBlue;
draw (        0.9250164787u,        0.0030079435u)--(        0.8529957280u,        0.1145552400u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7235486737u,        0.1441060105u)--(        0.8529957280u,        0.1145552400u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9250164787u,        0.0030079435u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9250164787u,        0.0030079435u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7235486737u,        0.1441060105u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7235486737u,        0.1441060105u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8529957280u,        0.1145552400u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8529957280u,        0.1145552400u) withcolor MidnightBlue;
draw (        0.9250164787u,        0.0030079435u)--(        0.8530059875u,        0.1145618634u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9279892810u,        0.2241396603u)--(        0.8530059875u,        0.1145618634u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9250164787u,        0.0030079435u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9250164787u,        0.0030079435u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9279892810u,        0.2241396603u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9279892810u,        0.2241396603u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8530059875u,        0.1145618634u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8530059875u,        0.1145618634u) withcolor MidnightBlue;
draw (        0.9250164787u,        0.0030079435u)--(        0.9999997722u,        0.1125857405u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9279892810u,        0.2241396603u)--(        0.9999997722u,        0.1125857405u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9250164787u,        0.0030079435u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9250164787u,        0.0030079435u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9279892810u,        0.2241396603u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9279892810u,        0.2241396603u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9999997722u,        0.1125857405u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9999997722u,        0.1125857405u) withcolor Green;
draw (        0.0968666128u,        0.1957863424u)--(        0.2295773458u,        0.1999874391u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3117482354u,        0.0956910691u)--(        0.2295773458u,        0.1999874391u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0968666128u,        0.1957863424u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0968666128u,        0.1957863424u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3117482354u,        0.0956910691u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3117482354u,        0.0956910691u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2295773458u,        0.1999874391u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2295773458u,        0.1999874391u) withcolor Green;
draw (        0.0968666128u,        0.1957863424u)--(        0.1790375023u,        0.0914899724u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3117482354u,        0.0956910691u)--(        0.1790375023u,        0.0914899724u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0968666128u,        0.1957863424u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0968666128u,        0.1957863424u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3117482354u,        0.0956910691u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3117482354u,        0.0956910691u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1790375023u,        0.0914899724u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1790375023u,        0.0914899724u) withcolor MidnightBlue;
draw (        0.0968666128u,        0.1957863424u)--(        0.0000193020u,        0.2866179793u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1050889823u,        0.3677977192u)--(        0.0000193020u,        0.2866179793u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0968666128u,        0.1957863424u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0968666128u,        0.1957863424u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1050889823u,        0.3677977192u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1050889823u,        0.3677977192u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000193020u,        0.2866179793u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000193020u,        0.2866179793u) withcolor Green;
draw (        0.3117482354u,        0.0956910691u)--(        0.3939164445u,        0.1999895509u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5195343715u,        0.2430002972u)--(        0.3939164445u,        0.1999895509u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3117482354u,        0.0956910691u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3117482354u,        0.0956910691u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5195343715u,        0.2430002972u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5195343715u,        0.2430002972u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3939164445u,        0.1999895509u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3939164445u,        0.1999895509u) withcolor Green;
draw (        0.3117482354u,        0.0956910691u)--(        0.4373661624u,        0.1387018155u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5195343715u,        0.2430002972u)--(        0.4373661624u,        0.1387018155u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3117482354u,        0.0956910691u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3117482354u,        0.0956910691u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5195343715u,        0.2430002972u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5195343715u,        0.2430002972u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4373661624u,        0.1387018155u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4373661624u,        0.1387018155u) withcolor MidnightBlue;
draw (        0.5195343715u,        0.2430002972u)--(        0.6516968541u,        0.2557621966u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7235486737u,        0.1441060105u)--(        0.6516968541u,        0.2557621966u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5195343715u,        0.2430002972u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5195343715u,        0.2430002972u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7235486737u,        0.1441060105u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7235486737u,        0.1441060105u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6516968541u,        0.2557621966u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6516968541u,        0.2557621966u) withcolor Green;
draw (        0.5195343715u,        0.2430002972u)--(        0.3882359022u,        0.2627613662u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3926632979u,        0.3954647424u)--(        0.3882359022u,        0.2627613662u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5195343715u,        0.2430002972u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5195343715u,        0.2430002972u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3926632979u,        0.3954647424u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3926632979u,        0.3954647424u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3882359022u,        0.2627613662u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3882359022u,        0.2627613662u) withcolor Green;
draw (        0.5195343715u,        0.2430002972u)--(        0.5239617672u,        0.3757036734u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3926632979u,        0.3954647424u)--(        0.5239617672u,        0.3757036734u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5195343715u,        0.2430002972u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5195343715u,        0.2430002972u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3926632979u,        0.3954647424u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3926632979u,        0.3954647424u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5239617672u,        0.3757036734u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5239617672u,        0.3757036734u) withcolor MidnightBlue;
draw (        0.5195343715u,        0.2430002972u)--(        0.5239665728u,        0.3757035130u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6110373210u,        0.4759455834u)--(        0.5239665728u,        0.3757035130u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5195343715u,        0.2430002972u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5195343715u,        0.2430002972u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6110373210u,        0.4759455834u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6110373210u,        0.4759455834u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5239665728u,        0.3757035130u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5239665728u,        0.3757035130u) withcolor MidnightBlue;
draw (        0.5195343715u,        0.2430002972u)--(        0.6066051198u,        0.3432423677u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6110373210u,        0.4759455834u)--(        0.6066051198u,        0.3432423677u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5195343715u,        0.2430002972u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5195343715u,        0.2430002972u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6110373210u,        0.4759455834u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6110373210u,        0.4759455834u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6066051198u,        0.3432423677u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6066051198u,        0.3432423677u) withcolor Green;
draw (        0.7235486737u,        0.1441060105u)--(        0.8529984214u,        0.1145670412u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9279892810u,        0.2241396603u)--(        0.8529984214u,        0.1145670412u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7235486737u,        0.1441060105u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7235486737u,        0.1441060105u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9279892810u,        0.2241396603u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9279892810u,        0.2241396603u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8529984214u,        0.1145670412u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8529984214u,        0.1145670412u) withcolor MidnightBlue;
draw (        0.9279892810u,        0.2241396603u)--(        0.9999922062u,        0.3356984639u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9338902108u,        0.4508519010u)--(        0.9999922062u,        0.3356984639u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9279892810u,        0.2241396603u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9279892810u,        0.2241396603u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9338902108u,        0.4508519010u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9338902108u,        0.4508519010u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9999922062u,        0.3356984639u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9999922062u,        0.3356984639u) withcolor Green;
draw (        0.1050889823u,        0.3677977192u)--(        0.2333333333u,        0.3333333333u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1050889823u,        0.3677977192u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1050889823u,        0.3677977192u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2333333333u,        0.3333333333u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2333333333u,        0.3333333333u) withcolor Yellow;
draw (        0.1050889823u,        0.3677977192u)--(        0.1937094292u,        0.4666724090u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2560161760u,        0.5839228127u)--(        0.1937094292u,        0.4666724090u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1050889823u,        0.3677977192u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1050889823u,        0.3677977192u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2560161760u,        0.5839228127u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2560161760u,        0.5839228127u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1937094292u,        0.4666724090u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1937094292u,        0.4666724090u) withcolor Green;
draw (        0.3926632979u,        0.3954647424u)--(        0.5239623014u,        0.3757072233u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6110373210u,        0.4759455834u)--(        0.5239623014u,        0.3757072233u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3926632979u,        0.3954647424u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3926632979u,        0.3954647424u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6110373210u,        0.4759455834u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6110373210u,        0.4759455834u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5239623014u,        0.3757072233u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5239623014u,        0.3757072233u) withcolor MidnightBlue;
draw (        0.3926632979u,        0.3954647424u)--(        0.3760678167u,        0.5272007583u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2560161760u,        0.5839228127u)--(        0.3760678167u,        0.5272007583u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3926632979u,        0.3954647424u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3926632979u,        0.3954647424u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2560161760u,        0.5839228127u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2560161760u,        0.5839228127u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3760678167u,        0.5272007583u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3760678167u,        0.5272007583u) withcolor MidnightBlue;
draw (        0.3926632979u,        0.3954647424u)--(        0.3760902084u,        0.5272035772u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4799239229u,        0.6099583210u)--(        0.3760902084u,        0.5272035772u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3926632979u,        0.3954647424u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3926632979u,        0.3954647424u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4799239229u,        0.6099583210u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4799239229u,        0.6099583210u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3760902084u,        0.5272035772u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3760902084u,        0.5272035772u) withcolor MidnightBlue;
draw (        0.6110373210u,        0.4759455834u)--(        0.6890101206u,        0.5834166512u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8198184288u,        0.5606356413u)--(        0.6890101206u,        0.5834166512u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6110373210u,        0.4759455834u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6110373210u,        0.4759455834u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8198184288u,        0.5606356413u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8198184288u,        0.5606356413u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6890101206u,        0.5834166512u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6890101206u,        0.5834166512u) withcolor Green;
draw (        0.9338902108u,        0.4508519010u)--(        0.8029341277u,        0.4289363312u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8198184288u,        0.5606356413u)--(        0.8029341277u,        0.4289363312u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9338902108u,        0.4508519010u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9338902108u,        0.4508519010u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8198184288u,        0.5606356413u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8198184288u,        0.5606356413u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8029341277u,        0.4289363312u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8029341277u,        0.4289363312u) withcolor Green;
draw (        0.9338902108u,        0.4508519010u)--(        0.9507745118u,        0.5825512110u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8198184288u,        0.5606356413u)--(        0.9507745118u,        0.5825512110u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9338902108u,        0.4508519010u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9338902108u,        0.4508519010u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8198184288u,        0.5606356413u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8198184288u,        0.5606356413u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9507745118u,        0.5825512110u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9507745118u,        0.5825512110u) withcolor Green;
draw (        0.0760653247u,        0.5578371882u)--(        0.0000000000u,        0.6666666667u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0760653247u,        0.5578371882u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0760653247u,        0.5578371882u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.6666666667u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.6666666667u) withcolor Yellow;
draw (        0.0760653247u,        0.5578371882u)--(        0.1521583504u,        0.6666472949u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2560161760u,        0.5839228127u)--(        0.1521583504u,        0.6666472949u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0760653247u,        0.5578371882u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0760653247u,        0.5578371882u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2560161760u,        0.5839228127u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2560161760u,        0.5839228127u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1521583504u,        0.6666472949u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1521583504u,        0.6666472949u) withcolor Green;
draw (        0.2560161760u,        0.5839228127u)--(        0.3598631126u,        0.6666609638u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4799239229u,        0.6099583210u)--(        0.3598631126u,        0.6666609638u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2560161760u,        0.5839228127u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2560161760u,        0.5839228127u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4799239229u,        0.6099583210u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4799239229u,        0.6099583210u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3598631126u,        0.6666609638u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3598631126u,        0.6666609638u) withcolor Green;
draw (        0.2560161760u,        0.5839228127u)--(        0.3760769863u,        0.5272201699u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4799239229u,        0.6099583210u)--(        0.3760769863u,        0.5272201699u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2560161760u,        0.5839228127u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2560161760u,        0.5839228127u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4799239229u,        0.6099583210u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4799239229u,        0.6099583210u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3760769863u,        0.5272201699u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3760769863u,        0.5272201699u) withcolor MidnightBlue;
draw (        0.4799239229u,        0.6099583210u)--(        0.6000000000u,        0.6666666667u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4799239229u,        0.6099583210u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4799239229u,        0.6099583210u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6000000000u,        0.6666666667u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6000000000u,        0.6666666667u) withcolor Yellow;
% r =         0.1327772115
endfig;
end;
