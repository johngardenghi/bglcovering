prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(10);
u := 20cm;
path pols[];
pair p[];
R := .1954663097u;
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

x1 = .9246607385u;
y1 = -.1529622055u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .9484053378u;
y2 = .2159378393u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .6788859427u;
y3 = .0580457175u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .8185544722u;
y4 = .4772116251u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .4951408689u;
y5 = .3447351419u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .4487582297u;
y6 = .5793483967u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .0666504649u;
y7 = -.1495734935u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .3453193974u;
y8 = .0385130353u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .0835924097u;
y9 = .2109245927u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = .1369184425u;
y10 = .5271597957u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.9246607385u,       -0.1529622055u)--(        1.0000000000u,       -0.3333333333u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9246607385u,       -0.1529622055u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9246607385u,       -0.1529622055u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,       -0.3333333333u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,       -0.3333333333u) withcolor Yellow;
draw (        0.9246607385u,       -0.1529622055u)--(        0.8730694296u,        0.0355727119u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9484053378u,        0.2159378393u)--(        0.8730694296u,        0.0355727119u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9246607385u,       -0.1529622055u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9246607385u,       -0.1529622055u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9484053378u,        0.2159378393u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9484053378u,        0.2159378393u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8730694296u,        0.0355727119u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8730694296u,        0.0355727119u) withcolor MidnightBlue;
draw (        0.9246607385u,       -0.1529622055u)--(        0.9999966468u,        0.0274029219u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9484053378u,        0.2159378393u)--(        0.9999966468u,        0.0274029219u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9246607385u,       -0.1529622055u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9246607385u,       -0.1529622055u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9484053378u,        0.2159378393u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9484053378u,        0.2159378393u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9999966468u,        0.0274029219u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9999966468u,        0.0274029219u) withcolor Green;
draw (        0.9246607385u,       -0.1529622055u)--(        0.7304910378u,       -0.1304854268u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6788859427u,        0.0580457175u)--(        0.7304910378u,       -0.1304854268u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9246607385u,       -0.1529622055u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9246607385u,       -0.1529622055u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6788859427u,        0.0580457175u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6788859427u,        0.0580457175u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7304910378u,       -0.1304854268u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7304910378u,       -0.1304854268u) withcolor Green;
draw (        0.9246607385u,       -0.1529622055u)--(        0.8730556433u,        0.0355689389u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6788859427u,        0.0580457175u)--(        0.8730556433u,        0.0355689389u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9246607385u,       -0.1529622055u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9246607385u,       -0.1529622055u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6788859427u,        0.0580457175u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6788859427u,        0.0580457175u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8730556433u,        0.0355689389u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8730556433u,        0.0355689389u) withcolor MidnightBlue;
draw (        0.9484053378u,        0.2159378393u)--(        0.8730566956u,        0.0355780312u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6788859427u,        0.0580457175u)--(        0.8730566956u,        0.0355780312u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9484053378u,        0.2159378393u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9484053378u,        0.2159378393u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6788859427u,        0.0580457175u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6788859427u,        0.0580457175u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8730566956u,        0.0355780312u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8730566956u,        0.0355780312u) withcolor MidnightBlue;
draw (        0.9484053378u,        0.2159378393u)--(        0.9999839132u,        0.4044762407u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8185544722u,        0.4772116251u)--(        0.9999839132u,        0.4044762407u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9484053378u,        0.2159378393u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9484053378u,        0.2159378393u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8185544722u,        0.4772116251u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8185544722u,        0.4772116251u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9999839132u,        0.4044762407u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9999839132u,        0.4044762407u) withcolor Green;
draw (        0.6788859427u,        0.0580457175u)--(        0.5061768255u,        0.1495806239u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4951408689u,        0.3447351419u)--(        0.5061768255u,        0.1495806239u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6788859427u,        0.0580457175u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6788859427u,        0.0580457175u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4951408689u,        0.3447351419u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4951408689u,        0.3447351419u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5061768255u,        0.1495806239u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5061768255u,        0.1495806239u) withcolor MidnightBlue;
draw (        0.6788859427u,        0.0580457175u)--(        0.5061713680u,        0.1495703260u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3453193974u,        0.0385130353u)--(        0.5061713680u,        0.1495703260u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6788859427u,        0.0580457175u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6788859427u,        0.0580457175u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3453193974u,        0.0385130353u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3453193974u,        0.0385130353u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5061713680u,        0.1495703260u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5061713680u,        0.1495703260u) withcolor MidnightBlue;
draw (        0.8185544722u,        0.4772116251u)--(        0.8666666667u,        0.6666666667u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8185544722u,        0.4772116251u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8185544722u,        0.4772116251u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8666666667u,        0.6666666667u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8666666667u,        0.6666666667u) withcolor Yellow;
draw (        0.8185544722u,        0.4772116251u)--(        0.6236498268u,        0.4920189699u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4951408689u,        0.3447351419u)--(        0.6236498268u,        0.4920189699u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8185544722u,        0.4772116251u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8185544722u,        0.4772116251u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4951408689u,        0.3447351419u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4951408689u,        0.3447351419u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6236498268u,        0.4920189699u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6236498268u,        0.4920189699u) withcolor MidnightBlue;
draw (        0.8185544722u,        0.4772116251u)--(        0.6236531278u,        0.4920623563u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4487582297u,        0.5793483967u)--(        0.6236531278u,        0.4920623563u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8185544722u,        0.4772116251u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8185544722u,        0.4772116251u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4487582297u,        0.5793483967u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4487582297u,        0.5793483967u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6236531278u,        0.4920623563u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6236531278u,        0.4920623563u) withcolor MidnightBlue;
draw (        0.4951408689u,        0.3447351419u)--(        0.3000000000u,        0.3333333333u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4951408689u,        0.3447351419u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4951408689u,        0.3447351419u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3000000000u,        0.3333333333u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3000000000u,        0.3333333333u) withcolor Yellow;
draw (        0.4951408689u,        0.3447351419u)--(        0.6236370390u,        0.4920301266u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4487582297u,        0.5793483967u)--(        0.6236370390u,        0.4920301266u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4951408689u,        0.3447351419u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4951408689u,        0.3447351419u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4487582297u,        0.5793483967u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4487582297u,        0.5793483967u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6236370390u,        0.4920301266u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6236370390u,        0.4920301266u) withcolor MidnightBlue;
draw (        0.4951408689u,        0.3447351419u)--(        0.5061647296u,        0.1495799403u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3453193974u,        0.0385130353u)--(        0.5061647296u,        0.1495799403u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4951408689u,        0.3447351419u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4951408689u,        0.3447351419u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3453193974u,        0.0385130353u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3453193974u,        0.0385130353u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5061647296u,        0.1495799403u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5061647296u,        0.1495799403u) withcolor MidnightBlue;
draw (        0.4487582297u,        0.5793483967u)--(        0.2738633316u,        0.6666344371u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1369184425u,        0.5271597957u)--(        0.2738633316u,        0.6666344371u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4487582297u,        0.5793483967u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4487582297u,        0.5793483967u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1369184425u,        0.5271597957u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1369184425u,        0.5271597957u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2738633316u,        0.6666344371u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2738633316u,        0.6666344371u) withcolor Green;
draw (        0.0666504649u,       -0.1495734935u)--(        0.0000000000u,       -0.3333333333u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0666504649u,       -0.1495734935u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0666504649u,       -0.1495734935u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,       -0.3333333333u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,       -0.3333333333u) withcolor Yellow;
draw (        0.0666504649u,       -0.1495734935u)--(        0.1333333333u,       -0.3333333333u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0666504649u,       -0.1495734935u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0666504649u,       -0.1495734935u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1333333333u,       -0.3333333333u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1333333333u,       -0.3333333333u) withcolor Yellow;
draw (        0.0666504649u,       -0.1495734935u)--(        0.1501839740u,        0.0271444487u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3453193974u,        0.0385130353u)--(        0.1501839740u,        0.0271444487u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0666504649u,       -0.1495734935u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0666504649u,       -0.1495734935u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3453193974u,        0.0385130353u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3453193974u,        0.0385130353u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1501839740u,        0.0271444487u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1501839740u,        0.0271444487u) withcolor MidnightBlue;
draw (        0.0666504649u,       -0.1495734935u)--(        0.2617858883u,       -0.1382049069u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3453193974u,        0.0385130353u)--(        0.2617858883u,       -0.1382049069u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0666504649u,       -0.1495734935u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0666504649u,       -0.1495734935u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3453193974u,        0.0385130353u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3453193974u,        0.0385130353u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2617858883u,       -0.1382049069u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2617858883u,       -0.1382049069u) withcolor Green;
draw (        0.0666504649u,       -0.1495734935u)--(        0.1501758298u,        0.0271482982u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0835924097u,        0.2109245927u)--(        0.1501758298u,        0.0271482982u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0666504649u,       -0.1495734935u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0666504649u,       -0.1495734935u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0835924097u,        0.2109245927u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0835924097u,        0.2109245927u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1501758298u,        0.0271482982u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1501758298u,        0.0271482982u) withcolor MidnightBlue;
draw (        0.3453193974u,        0.0385130353u)--(        0.1501835861u,        0.0271511086u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0835924097u,        0.2109245927u)--(        0.1501835861u,        0.0271511086u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3453193974u,        0.0385130353u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3453193974u,        0.0385130353u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0835924097u,        0.2109245927u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0835924097u,        0.2109245927u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1501835861u,        0.0271511086u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1501835861u,        0.0271511086u) withcolor MidnightBlue;
draw (        0.0835924097u,        0.2109245927u)--(        0.2204802449u,        0.3504552303u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1369184425u,        0.5271597957u)--(        0.2204802449u,        0.3504552303u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0835924097u,        0.2109245927u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0835924097u,        0.2109245927u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1369184425u,        0.5271597957u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1369184425u,        0.5271597957u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2204802449u,        0.3504552303u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2204802449u,        0.3504552303u) withcolor Green;
draw (        0.1369184425u,        0.5271597957u)--(        0.0000000000u,        0.6666666667u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1369184425u,        0.5271597957u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1369184425u,        0.5271597957u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.6666666667u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.6666666667u) withcolor Yellow;
% r =         0.1954663097
endfig;
end;
