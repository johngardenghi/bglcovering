prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(30);
u := 20cm;
path pols[];
pair p[];
R := .1097578738u;
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

x1 = .9698343386u;
y1 = .2036968162u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .9056059394u;
y2 = .3652557745u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .9181280619u;
y3 = .4943778321u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .9754085858u;
y4 = -.1009716655u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .9003760147u;
y5 = .0520620791u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .7937904263u;
y6 = .1863718461u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .7274112963u;
y7 = .4909871412u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .8200734908u;
y8 = .6014648575u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .9099050159u;
y9 = -.2706277291u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = .7906017260u;
y10 = -.1081946508u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;

x11 = .6963168780u;
y11 = .0254028677u;
draw fullcircle scaled D shifted z11 withcolor RoyalBlue withpen pencircle scaled 1.0;

x12 = .6137589764u;
y12 = .2070517684u;
draw fullcircle scaled D shifted z12 withcolor RoyalBlue withpen pencircle scaled 1.0;

x13 = .5576463832u;
y13 = .4056269100u;
draw fullcircle scaled D shifted z13 withcolor RoyalBlue withpen pencircle scaled 1.0;

x14 = .5399619167u;
y14 = .5747744183u;
draw fullcircle scaled D shifted z14 withcolor RoyalBlue withpen pencircle scaled 1.0;

x15 = .5065914555u;
y15 = .0719135323u;
draw fullcircle scaled D shifted z15 withcolor RoyalBlue withpen pencircle scaled 1.0;

x16 = .4436694178u;
y16 = .2676294293u;
draw fullcircle scaled D shifted z16 withcolor RoyalBlue withpen pencircle scaled 1.0;

x17 = .3737924749u;
y17 = .4544520445u;
draw fullcircle scaled D shifted z17 withcolor RoyalBlue withpen pencircle scaled 1.0;

x18 = .3739378977u;
y18 = .6381290719u;
draw fullcircle scaled D shifted z18 withcolor RoyalBlue withpen pencircle scaled 1.0;

x19 = .3197331378u;
y19 = -.0030480941u;
draw fullcircle scaled D shifted z19 withcolor RoyalBlue withpen pencircle scaled 1.0;

x20 = .3287603085u;
y20 = .1611041873u;
draw fullcircle scaled D shifted z20 withcolor RoyalBlue withpen pencircle scaled 1.0;

x21 = .2567148418u;
y21 = .3571795641u;
draw fullcircle scaled D shifted z21 withcolor RoyalBlue withpen pencircle scaled 1.0;

x22 = .2210886587u;
y22 = .5674005019u;
draw fullcircle scaled D shifted z22 withcolor RoyalBlue withpen pencircle scaled 1.0;

x23 = .1814367963u;
y23 = -.1153681548u;
draw fullcircle scaled D shifted z23 withcolor RoyalBlue withpen pencircle scaled 1.0;

x24 = .1862819111u;
y24 = .0976996048u;
draw fullcircle scaled D shifted z24 withcolor RoyalBlue withpen pencircle scaled 1.0;

x25 = .0785464688u;
y25 = .2861832645u;
draw fullcircle scaled D shifted z25 withcolor RoyalBlue withpen pencircle scaled 1.0;

x26 = .0693241769u;
y26 = .4479599365u;
draw fullcircle scaled D shifted z26 withcolor RoyalBlue withpen pencircle scaled 1.0;

x27 = .0870979687u;
y27 = .5998779292u;
draw fullcircle scaled D shifted z27 withcolor RoyalBlue withpen pencircle scaled 1.0;

x28 = .0674193742u;
y28 = -.2466921418u;
draw fullcircle scaled D shifted z28 withcolor RoyalBlue withpen pencircle scaled 1.0;

x29 = .0724245090u;
y29 = -.0775630941u;
draw fullcircle scaled D shifted z29 withcolor RoyalBlue withpen pencircle scaled 1.0;

x30 = .0397641253u;
y30 = .1072155353u;
draw fullcircle scaled D shifted z30 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.9698343386u,        0.2036968162u)--(        0.8754518861u,        0.2597213065u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9056059394u,        0.3652557745u)--(        0.8754518861u,        0.2597213065u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9698343386u,        0.2036968162u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9698343386u,        0.2036968162u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9056059394u,        0.3652557745u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9056059394u,        0.3652557745u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8754518861u,        0.2597213065u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8754518861u,        0.2597213065u) withcolor MidnightBlue;
draw (        0.9698343386u,        0.2036968162u)--(        0.9999883919u,        0.3092312841u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9056059394u,        0.3652557745u)--(        0.9999883919u,        0.3092312841u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9698343386u,        0.2036968162u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9698343386u,        0.2036968162u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9056059394u,        0.3652557745u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9056059394u,        0.3652557745u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9999883919u,        0.3092312841u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9999883919u,        0.3092312841u) withcolor Green;
draw (        0.9698343386u,        0.2036968162u)--(        0.9999836498u,        0.0981609933u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9003760147u,        0.0520620791u)--(        0.9999836498u,        0.0981609933u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9698343386u,        0.2036968162u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9698343386u,        0.2036968162u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9003760147u,        0.0520620791u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9003760147u,        0.0520620791u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9999836498u,        0.0981609933u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9999836498u,        0.0981609933u) withcolor Green;
draw (        0.9698343386u,        0.2036968162u)--(        0.8754471439u,        0.2597133168u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7937904263u,        0.1863718461u)--(        0.8754471439u,        0.2597133168u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9698343386u,        0.2036968162u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9698343386u,        0.2036968162u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7937904263u,        0.1863718461u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7937904263u,        0.1863718461u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8754471439u,        0.2597133168u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8754471439u,        0.2597133168u) withcolor MidnightBlue;
draw (        0.9056059394u,        0.3652557745u)--(        0.8237397008u,        0.4383632970u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9181280619u,        0.4943778321u)--(        0.8237397008u,        0.4383632970u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9056059394u,        0.3652557745u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9056059394u,        0.3652557745u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9181280619u,        0.4943778321u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9181280619u,        0.4943778321u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8237397008u,        0.4383632970u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8237397008u,        0.4383632970u) withcolor MidnightBlue;
draw (        0.9056059394u,        0.3652557745u)--(        0.9999943006u,        0.4212703097u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9181280619u,        0.4943778321u)--(        0.9999943006u,        0.4212703097u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9056059394u,        0.3652557745u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9056059394u,        0.3652557745u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9181280619u,        0.4943778321u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9181280619u,        0.4943778321u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9999943006u,        0.4212703097u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9999943006u,        0.4212703097u) withcolor Green;
draw (        0.9056059394u,        0.3652557745u)--(        0.8754358494u,        0.2597258900u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7937904263u,        0.1863718461u)--(        0.8754358494u,        0.2597258900u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9056059394u,        0.3652557745u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9056059394u,        0.3652557745u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7937904263u,        0.1863718461u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7937904263u,        0.1863718461u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8754358494u,        0.2597258900u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8754358494u,        0.2597258900u) withcolor MidnightBlue;
draw (        0.9056059394u,        0.3652557745u)--(        0.8237177557u,        0.4383387154u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7274112963u,        0.4909871412u)--(        0.8237177557u,        0.4383387154u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9056059394u,        0.3652557745u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9056059394u,        0.3652557745u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7274112963u,        0.4909871412u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7274112963u,        0.4909871412u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8237177557u,        0.4383387154u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8237177557u,        0.4383387154u) withcolor MidnightBlue;
draw (        0.9181280619u,        0.4943778321u)--(        0.8237352694u,        0.4383707648u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7274112963u,        0.4909871412u)--(        0.8237352694u,        0.4383707648u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9181280619u,        0.4943778321u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9181280619u,        0.4943778321u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7274112963u,        0.4909871412u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7274112963u,        0.4909871412u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8237352694u,        0.4383707648u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8237352694u,        0.4383707648u) withcolor MidnightBlue;
draw (        0.9181280619u,        0.4943778321u)--(        0.9298122716u,        0.6035120161u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8200734908u,        0.6014648575u)--(        0.9298122716u,        0.6035120161u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9181280619u,        0.4943778321u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9181280619u,        0.4943778321u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8200734908u,        0.6014648575u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8200734908u,        0.6014648575u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9298122716u,        0.6035120161u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9298122716u,        0.6035120161u) withcolor Green;
draw (        0.9754085858u,       -0.1009716655u)--(        0.9999988958u,        0.0059961218u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9003760147u,        0.0520620791u)--(        0.9999988958u,        0.0059961218u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9754085858u,       -0.1009716655u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9754085858u,       -0.1009716655u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9003760147u,        0.0520620791u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9003760147u,        0.0520620791u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9999988958u,        0.0059961218u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9999988958u,        0.0059961218u) withcolor Green;
draw (        0.9754085858u,       -0.1009716655u)--(        0.9999992005u,       -0.2079393828u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9099050159u,       -0.2706277291u)--(        0.9999992005u,       -0.2079393828u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9754085858u,       -0.1009716655u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9754085858u,       -0.1009716655u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9099050159u,       -0.2706277291u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9099050159u,       -0.2706277291u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9999992005u,       -0.2079393828u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9999992005u,       -0.2079393828u) withcolor Green;
draw (        0.9754085858u,       -0.1009716655u)--(        0.8853144012u,       -0.1636600118u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9099050159u,       -0.2706277291u)--(        0.8853144012u,       -0.1636600118u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9754085858u,       -0.1009716655u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9754085858u,       -0.1009716655u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9099050159u,       -0.2706277291u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9099050159u,       -0.2706277291u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8853144012u,       -0.1636600118u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8853144012u,       -0.1636600118u) withcolor MidnightBlue;
draw (        0.9754085858u,       -0.1009716655u)--(        0.8853140966u,       -0.1636595741u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7906017260u,       -0.1081946508u)--(        0.8853140966u,       -0.1636595741u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9754085858u,       -0.1009716655u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9754085858u,       -0.1009716655u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7906017260u,       -0.1081946508u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7906017260u,       -0.1081946508u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8853140966u,       -0.1636595741u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8853140966u,       -0.1636595741u) withcolor MidnightBlue;
draw (        0.9003760147u,        0.0520620791u)--(        0.7933995629u,        0.0766146682u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7937904263u,        0.1863718461u)--(        0.7933995629u,        0.0766146682u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9003760147u,        0.0520620791u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9003760147u,        0.0520620791u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7937904263u,        0.1863718461u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7937904263u,        0.1863718461u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7933995629u,        0.0766146682u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7933995629u,        0.0766146682u) withcolor MidnightBlue;
draw (        0.9003760147u,        0.0520620791u)--(        0.8033114031u,        0.0008248657u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7906017260u,       -0.1081946508u)--(        0.8033114031u,        0.0008248657u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9003760147u,        0.0520620791u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9003760147u,        0.0520620791u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7906017260u,       -0.1081946508u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7906017260u,       -0.1081946508u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8033114031u,        0.0008248657u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8033114031u,        0.0008248657u) withcolor MidnightBlue;
draw (        0.9003760147u,        0.0520620791u)--(        0.8032947434u,        0.0008564385u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6963168780u,        0.0254028677u)--(        0.8032947434u,        0.0008564385u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9003760147u,        0.0520620791u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9003760147u,        0.0520620791u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6963168780u,        0.0254028677u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6963168780u,        0.0254028677u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8032947434u,        0.0008564385u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8032947434u,        0.0008564385u) withcolor MidnightBlue;
draw (        0.9003760147u,        0.0520620791u)--(        0.7933981493u,        0.0766085083u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6963168780u,        0.0254028677u)--(        0.7933981493u,        0.0766085083u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9003760147u,        0.0520620791u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9003760147u,        0.0520620791u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6963168780u,        0.0254028677u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6963168780u,        0.0254028677u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7933981493u,        0.0766085083u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7933981493u,        0.0766085083u) withcolor MidnightBlue;
draw (        0.7937904263u,        0.1863718461u)--(        0.7933948912u,        0.0766146849u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6963168780u,        0.0254028677u)--(        0.7933948912u,        0.0766146849u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7937904263u,        0.1863718461u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7937904263u,        0.1863718461u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6963168780u,        0.0254028677u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6963168780u,        0.0254028677u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7933948912u,        0.0766146849u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7933948912u,        0.0766146849u) withcolor MidnightBlue;
draw (        0.7937904263u,        0.1863718461u)--(        0.6967124131u,        0.1351600288u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6963168780u,        0.0254028677u)--(        0.6967124131u,        0.1351600288u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7937904263u,        0.1863718461u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7937904263u,        0.1863718461u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6963168780u,        0.0254028677u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6963168780u,        0.0254028677u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6967124131u,        0.1351600288u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6967124131u,        0.1351600288u) withcolor MidnightBlue;
draw (        0.7937904263u,        0.1863718461u)--(        0.6967057902u,        0.1351725854u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6137589764u,        0.2070517684u)--(        0.6967057902u,        0.1351725854u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7937904263u,        0.1863718461u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7937904263u,        0.1863718461u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6137589764u,        0.2070517684u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6137589764u,        0.2070517684u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6967057902u,        0.1351725854u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6967057902u,        0.1351725854u) withcolor MidnightBlue;
draw (        0.7274112963u,        0.4909871412u)--(        0.7103348824u,        0.5994084811u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8200734908u,        0.6014648575u)--(        0.7103348824u,        0.5994084811u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7274112963u,        0.4909871412u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7274112963u,        0.4909871412u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8200734908u,        0.6014648575u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8200734908u,        0.6014648575u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7103348824u,        0.5994084811u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7103348824u,        0.5994084811u) withcolor Green;
draw (        0.7274112963u,        0.4909871412u)--(        0.6178412660u,        0.4974058125u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5576463832u,        0.4056269100u)--(        0.6178412660u,        0.4974058125u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7274112963u,        0.4909871412u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7274112963u,        0.4909871412u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5576463832u,        0.4056269100u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5576463832u,        0.4056269100u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6178412660u,        0.4974058125u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6178412660u,        0.4974058125u) withcolor MidnightBlue;
draw (        0.7274112963u,        0.4909871412u)--(        0.6178430019u,        0.4974353761u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5399619167u,        0.5747744183u)--(        0.6178430019u,        0.4974353761u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7274112963u,        0.4909871412u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7274112963u,        0.4909871412u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5399619167u,        0.5747744183u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5399619167u,        0.5747744183u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6178430019u,        0.4974353761u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6178430019u,        0.4974353761u) withcolor MidnightBlue;
draw (        0.9099050159u,       -0.2706277291u)--(        1.0000000000u,       -0.3333333333u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9099050159u,       -0.2706277291u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9099050159u,       -0.2706277291u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,       -0.3333333333u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,       -0.3333333333u) withcolor Yellow;
draw (        0.9099050159u,       -0.2706277291u)--(        0.8151929889u,       -0.2151622190u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7906017260u,       -0.1081946508u)--(        0.8151929889u,       -0.2151622190u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9099050159u,       -0.2706277291u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9099050159u,       -0.2706277291u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7906017260u,       -0.1081946508u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7906017260u,       -0.1081946508u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8151929889u,       -0.2151622190u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8151929889u,       -0.2151622190u) withcolor Green;
draw (        0.9099050159u,       -0.2706277291u)--(        0.8853137530u,       -0.1636601608u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7906017260u,       -0.1081946508u)--(        0.8853137530u,       -0.1636601608u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9099050159u,       -0.2706277291u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9099050159u,       -0.2706277291u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7906017260u,       -0.1081946508u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7906017260u,       -0.1081946508u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8853137530u,       -0.1636601608u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8853137530u,       -0.1636601608u) withcolor MidnightBlue;
draw (        0.7906017260u,       -0.1081946508u)--(        0.6836304869u,       -0.0836193610u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6963168780u,        0.0254028677u)--(        0.6836304869u,       -0.0836193610u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7906017260u,       -0.1081946508u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7906017260u,       -0.1081946508u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6963168780u,        0.0254028677u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6963168780u,        0.0254028677u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6836304869u,       -0.0836193610u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6836304869u,       -0.0836193610u) withcolor Green;
draw (        0.7906017260u,       -0.1081946508u)--(        0.8032881171u,        0.0008275779u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6963168780u,        0.0254028677u)--(        0.8032881171u,        0.0008275779u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7906017260u,       -0.1081946508u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7906017260u,       -0.1081946508u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6963168780u,        0.0254028677u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6963168780u,        0.0254028677u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8032881171u,        0.0008275779u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8032881171u,        0.0008275779u) withcolor MidnightBlue;
draw (        0.6963168780u,        0.0254028677u)--(        0.6133808937u,        0.0972945458u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6137589764u,        0.2070517684u)--(        0.6133808937u,        0.0972945458u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6963168780u,        0.0254028677u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6963168780u,        0.0254028677u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6137589764u,        0.2070517684u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6137589764u,        0.2070517684u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6133808937u,        0.0972945458u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6133808937u,        0.0972945458u) withcolor MidnightBlue;
draw (        0.6963168780u,        0.0254028677u)--(        0.6966949607u,        0.1351600903u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6137589764u,        0.2070517684u)--(        0.6966949607u,        0.1351600903u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6963168780u,        0.0254028677u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6963168780u,        0.0254028677u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6137589764u,        0.2070517684u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6137589764u,        0.2070517684u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6966949607u,        0.1351600903u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6966949607u,        0.1351600903u) withcolor MidnightBlue;
draw (        0.6963168780u,        0.0254028677u)--(        0.6133757892u,        0.0972886566u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5065914555u,        0.0719135323u)--(        0.6133757892u,        0.0972886566u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6963168780u,        0.0254028677u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6963168780u,        0.0254028677u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5065914555u,        0.0719135323u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5065914555u,        0.0719135323u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6133757892u,        0.0972886566u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6133757892u,        0.0972886566u) withcolor MidnightBlue;
draw (        0.6137589764u,        0.2070517684u)--(        0.5496753108u,        0.2961588650u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5576463832u,        0.4056269100u)--(        0.5496753108u,        0.2961588650u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6137589764u,        0.2070517684u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6137589764u,        0.2070517684u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5576463832u,        0.4056269100u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5576463832u,        0.4056269100u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5496753108u,        0.2961588650u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5496753108u,        0.2961588650u) withcolor MidnightBlue;
draw (        0.6137589764u,        0.2070517684u)--(        0.6133743842u,        0.0972945684u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5065914555u,        0.0719135323u)--(        0.6133743842u,        0.0972945684u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6137589764u,        0.2070517684u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6137589764u,        0.2070517684u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5065914555u,        0.0719135323u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5065914555u,        0.0719135323u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6133743842u,        0.0972945684u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6133743842u,        0.0972945684u) withcolor MidnightBlue;
draw (        0.6137589764u,        0.2070517684u)--(        0.5496579721u,        0.2961463929u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4436694178u,        0.2676294293u)--(        0.5496579721u,        0.2961463929u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6137589764u,        0.2070517684u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6137589764u,        0.2070517684u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4436694178u,        0.2676294293u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4436694178u,        0.2676294293u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5496579721u,        0.2961463929u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5496579721u,        0.2961463929u) withcolor MidnightBlue;
draw (        0.5576463832u,        0.4056269100u)--(        0.4797837652u,        0.4829845443u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5399619167u,        0.5747744183u)--(        0.4797837652u,        0.4829845443u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5576463832u,        0.4056269100u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5576463832u,        0.4056269100u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5399619167u,        0.5747744183u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5399619167u,        0.5747744183u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4797837652u,        0.4829845443u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4797837652u,        0.4829845443u) withcolor MidnightBlue;
draw (        0.5576463832u,        0.4056269100u)--(        0.6178245347u,        0.4974167839u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5399619167u,        0.5747744183u)--(        0.6178245347u,        0.4974167839u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5576463832u,        0.4056269100u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5576463832u,        0.4056269100u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5399619167u,        0.5747744183u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5399619167u,        0.5747744183u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6178245347u,        0.4974167839u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6178245347u,        0.4974167839u) withcolor MidnightBlue;
draw (        0.5576463832u,        0.4056269100u)--(        0.5496542013u,        0.2961604042u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4436694178u,        0.2676294293u)--(        0.5496542013u,        0.2961604042u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5576463832u,        0.4056269100u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5576463832u,        0.4056269100u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4436694178u,        0.2676294293u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4436694178u,        0.2676294293u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5496542013u,        0.2961604042u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5496542013u,        0.2961604042u) withcolor MidnightBlue;
draw (        0.5576463832u,        0.4056269100u)--(        0.4516615997u,        0.3770959351u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4436694178u,        0.2676294293u)--(        0.4516615997u,        0.3770959351u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5576463832u,        0.4056269100u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5576463832u,        0.4056269100u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4436694178u,        0.2676294293u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4436694178u,        0.2676294293u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4516615997u,        0.3770959351u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4516615997u,        0.3770959351u) withcolor MidnightBlue;
draw (        0.5576463832u,        0.4056269100u)--(        0.4516605355u,        0.3770998887u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3737924749u,        0.4544520445u)--(        0.4516605355u,        0.3770998887u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5576463832u,        0.4056269100u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5576463832u,        0.4056269100u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3737924749u,        0.4544520445u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3737924749u,        0.4544520445u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4516605355u,        0.3770998887u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4516605355u,        0.3770998887u) withcolor MidnightBlue;
draw (        0.5576463832u,        0.4056269100u)--(        0.4797783227u,        0.4829790659u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3737924749u,        0.4544520445u)--(        0.4797783227u,        0.4829790659u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5576463832u,        0.4056269100u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5576463832u,        0.4056269100u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3737924749u,        0.4544520445u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3737924749u,        0.4544520445u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4797783227u,        0.4829790659u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4797783227u,        0.4829790659u) withcolor MidnightBlue;
draw (        0.5399619167u,        0.5747744183u)--(        0.6000000000u,        0.6666666667u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5399619167u,        0.5747744183u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5399619167u,        0.5747744183u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6000000000u,        0.6666666667u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6000000000u,        0.6666666667u) withcolor Yellow;
draw (        0.5399619167u,        0.5747744183u)--(        0.4797753650u,        0.4829900521u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3737924749u,        0.4544520445u)--(        0.4797753650u,        0.4829900521u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5399619167u,        0.5747744183u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5399619167u,        0.5747744183u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3737924749u,        0.4544520445u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3737924749u,        0.4544520445u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4797753650u,        0.4829900521u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4797753650u,        0.4829900521u) withcolor MidnightBlue;
draw (        0.5399619167u,        0.5747744183u)--(        0.4339790266u,        0.5462364107u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3737924749u,        0.4544520445u)--(        0.4339790266u,        0.5462364107u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5399619167u,        0.5747744183u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5399619167u,        0.5747744183u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3737924749u,        0.4544520445u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3737924749u,        0.4544520445u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4339790266u,        0.5462364107u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4339790266u,        0.5462364107u) withcolor MidnightBlue;
draw (        0.5399619167u,        0.5747744183u)--(        0.4339760274u,        0.5462475512u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3739378977u,        0.6381290719u)--(        0.4339760274u,        0.5462475512u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5399619167u,        0.5747744183u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5399619167u,        0.5747744183u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3739378977u,        0.6381290719u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3739378977u,        0.6381290719u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4339760274u,        0.5462475512u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4339760274u,        0.5462475512u) withcolor MidnightBlue;
draw (        0.5399619167u,        0.5747744183u)--(        0.4799237869u,        0.6666559389u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3739378977u,        0.6381290719u)--(        0.4799237869u,        0.6666559389u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5399619167u,        0.5747744183u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5399619167u,        0.5747744183u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3739378977u,        0.6381290719u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3739378977u,        0.6381290719u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4799237869u,        0.6666559389u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4799237869u,        0.6666559389u) withcolor Green;
draw (        0.5065914555u,        0.0719135323u)--(        0.4384957864u,        0.1579935574u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4436694178u,        0.2676294293u)--(        0.4384957864u,        0.1579935574u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5065914555u,        0.0719135323u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5065914555u,        0.0719135323u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4436694178u,        0.2676294293u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4436694178u,        0.2676294293u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4384957864u,        0.1579935574u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4384957864u,        0.1579935574u) withcolor MidnightBlue;
draw (        0.5065914555u,        0.0719135323u)--(        0.3968777109u,        0.0750256256u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3197331378u,       -0.0030480941u)--(        0.3968777109u,        0.0750256256u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5065914555u,        0.0719135323u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5065914555u,        0.0719135323u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3197331378u,       -0.0030480941u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3197331378u,       -0.0030480941u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3968777109u,        0.0750256256u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3968777109u,        0.0750256256u) withcolor MidnightBlue;
draw (        0.5065914555u,        0.0719135323u)--(        0.3968781687u,        0.0750417217u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3287603085u,        0.1611041873u)--(        0.3968781687u,        0.0750417217u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5065914555u,        0.0719135323u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5065914555u,        0.0719135323u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3287603085u,        0.1611041873u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3287603085u,        0.1611041873u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3968781687u,        0.0750417217u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3968781687u,        0.0750417217u) withcolor MidnightBlue;
draw (        0.5065914555u,        0.0719135323u)--(        0.4384735953u,        0.1579759979u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3287603085u,        0.1611041873u)--(        0.4384735953u,        0.1579759979u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5065914555u,        0.0719135323u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5065914555u,        0.0719135323u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3287603085u,        0.1611041873u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3287603085u,        0.1611041873u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4384735953u,        0.1579759979u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4384735953u,        0.1579759979u) withcolor MidnightBlue;
draw (        0.4436694178u,        0.2676294293u)--(        0.3658049473u,        0.3449851990u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3737924749u,        0.4544520445u)--(        0.3658049473u,        0.3449851990u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4436694178u,        0.2676294293u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4436694178u,        0.2676294293u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3737924749u,        0.4544520445u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3737924749u,        0.4544520445u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3658049473u,        0.3449851990u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3658049473u,        0.3449851990u) withcolor MidnightBlue;
draw (        0.4436694178u,        0.2676294293u)--(        0.4516569454u,        0.3770962748u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3737924749u,        0.4544520445u)--(        0.4516569454u,        0.3770962748u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4436694178u,        0.2676294293u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4436694178u,        0.2676294293u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3737924749u,        0.4544520445u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3737924749u,        0.4544520445u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4516569454u,        0.3770962748u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4516569454u,        0.3770962748u) withcolor MidnightBlue;
draw (        0.4436694178u,        0.2676294293u)--(        0.4384741236u,        0.1579945818u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3287603085u,        0.1611041873u)--(        0.4384741236u,        0.1579945818u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4436694178u,        0.2676294293u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4436694178u,        0.2676294293u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3287603085u,        0.1611041873u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3287603085u,        0.1611041873u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4384741236u,        0.1579945818u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4384741236u,        0.1579945818u) withcolor MidnightBlue;
draw (        0.4436694178u,        0.2676294293u)--(        0.3657917097u,        0.3449718721u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2567148418u,        0.3571795641u)--(        0.3657917097u,        0.3449718721u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4436694178u,        0.2676294293u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4436694178u,        0.2676294293u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2567148418u,        0.3571795641u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2567148418u,        0.3571795641u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3657917097u,        0.3449718721u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3657917097u,        0.3449718721u) withcolor MidnightBlue;
draw (        0.3737924749u,        0.4544520445u)--(        0.4339690191u,        0.5462429722u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3739378977u,        0.6381290719u)--(        0.4339690191u,        0.5462429722u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3737924749u,        0.4544520445u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3737924749u,        0.4544520445u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3739378977u,        0.6381290719u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3739378977u,        0.6381290719u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4339690191u,        0.5462429722u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4339690191u,        0.5462429722u) withcolor MidnightBlue;
draw (        0.3737924749u,        0.4544520445u)--(        0.3657932955u,        0.3449860498u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2567148418u,        0.3571795641u)--(        0.3657932955u,        0.3449860498u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3737924749u,        0.4544520445u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3737924749u,        0.4544520445u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2567148418u,        0.3571795641u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2567148418u,        0.3571795641u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3657932955u,        0.3449860498u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3657932955u,        0.3449860498u) withcolor MidnightBlue;
draw (        0.3737924749u,        0.4544520445u)--(        0.2647140212u,        0.4666455588u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2567148418u,        0.3571795641u)--(        0.2647140212u,        0.4666455588u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3737924749u,        0.4544520445u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3737924749u,        0.4544520445u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2567148418u,        0.3571795641u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2567148418u,        0.3571795641u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2647140212u,        0.4666455588u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2647140212u,        0.4666455588u) withcolor MidnightBlue;
draw (        0.3737924749u,        0.4544520445u)--(        0.2647186592u,        0.4666869778u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2210886587u,        0.5674005019u)--(        0.2647186592u,        0.4666869778u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3737924749u,        0.4544520445u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3737924749u,        0.4544520445u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2210886587u,        0.5674005019u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2210886587u,        0.5674005019u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2647186592u,        0.4666869778u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2647186592u,        0.4666869778u) withcolor MidnightBlue;
draw (        0.3739378977u,        0.6381290719u)--(        0.2679507760u,        0.6666513599u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2210886587u,        0.5674005019u)--(        0.2679507760u,        0.6666513599u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3739378977u,        0.6381290719u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3739378977u,        0.6381290719u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2210886587u,        0.5674005019u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2210886587u,        0.5674005019u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2679507760u,        0.6666513599u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2679507760u,        0.6666513599u) withcolor Green;
draw (        0.3197331378u,       -0.0030480941u)--(        0.3968688687u,        0.0750343616u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3287603085u,        0.1611041873u)--(        0.3968688687u,        0.0750343616u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3197331378u,       -0.0030480941u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3197331378u,       -0.0030480941u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3287603085u,        0.1611041873u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3287603085u,        0.1611041873u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3968688687u,        0.0750343616u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3968688687u,        0.0750343616u) withcolor MidnightBlue;
draw (        0.3197331378u,       -0.0030480941u)--(        0.2910086445u,       -0.1089805931u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1814367963u,       -0.1153681548u)--(        0.2910086445u,       -0.1089805931u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3197331378u,       -0.0030480941u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3197331378u,       -0.0030480941u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1814367963u,       -0.1153681548u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1814367963u,       -0.1153681548u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2910086445u,       -0.1089805931u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2910086445u,       -0.1089805931u) withcolor Green;
draw (        0.3197331378u,       -0.0030480941u)--(        0.2101612895u,       -0.0094356557u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1814367963u,       -0.1153681548u)--(        0.2101612895u,       -0.0094356557u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3197331378u,       -0.0030480941u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3197331378u,       -0.0030480941u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1814367963u,       -0.1153681548u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1814367963u,       -0.1153681548u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2101612895u,       -0.0094356557u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2101612895u,       -0.0094356557u) withcolor MidnightBlue;
draw (        0.3197331378u,       -0.0030480941u)--(        0.2101609143u,       -0.0094292152u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1862819111u,        0.0976996048u)--(        0.2101609143u,       -0.0094292152u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3197331378u,       -0.0030480941u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3197331378u,       -0.0030480941u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1862819111u,        0.0976996048u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1862819111u,        0.0976996048u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2101609143u,       -0.0094292152u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2101609143u,       -0.0094292152u) withcolor MidnightBlue;
draw (        0.3287603085u,        0.1611041873u)--(        0.2261156586u,        0.1999740563u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1862819111u,        0.0976996048u)--(        0.2261156586u,        0.1999740563u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3287603085u,        0.1611041873u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3287603085u,        0.1611041873u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1862819111u,        0.0976996048u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1862819111u,        0.0976996048u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2261156586u,        0.1999740563u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2261156586u,        0.1999740563u) withcolor Green;
draw (        0.2567148418u,        0.3571795641u)--(        0.2646360453u,        0.4666512290u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2210886587u,        0.5674005019u)--(        0.2646360453u,        0.4666512290u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2567148418u,        0.3571795641u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2567148418u,        0.3571795641u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2210886587u,        0.5674005019u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2210886587u,        0.5674005019u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2646360453u,        0.4666512290u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2646360453u,        0.4666512290u) withcolor MidnightBlue;
draw (        0.2567148418u,        0.3571795641u)--(        0.1873947919u,        0.2720824340u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0785464688u,        0.2861832645u)--(        0.1873947919u,        0.2720824340u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2567148418u,        0.3571795641u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2567148418u,        0.3571795641u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0785464688u,        0.2861832645u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0785464688u,        0.2861832645u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1873947919u,        0.2720824340u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1873947919u,        0.2720824340u) withcolor Green;
draw (        0.2567148418u,        0.3571795641u)--(        0.1478665186u,        0.3712803946u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0785464688u,        0.2861832645u)--(        0.1478665186u,        0.3712803946u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2567148418u,        0.3571795641u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2567148418u,        0.3571795641u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0785464688u,        0.2861832645u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0785464688u,        0.2861832645u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1478665186u,        0.3712803946u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1478665186u,        0.3712803946u) withcolor MidnightBlue;
draw (        0.2567148418u,        0.3571795641u)--(        0.1478683321u,        0.3712943865u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0693241769u,        0.4479599365u)--(        0.1478683321u,        0.3712943865u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2567148418u,        0.3571795641u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2567148418u,        0.3571795641u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0693241769u,        0.4479599365u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0693241769u,        0.4479599365u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1478683321u,        0.3712943865u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1478683321u,        0.3712943865u) withcolor MidnightBlue;
draw (        0.2210886587u,        0.5674005019u)--(        0.1774737670u,        0.4666804339u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0693241769u,        0.4479599365u)--(        0.1774737670u,        0.4666804339u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2210886587u,        0.5674005019u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2210886587u,        0.5674005019u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0693241769u,        0.4479599365u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0693241769u,        0.4479599365u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1774737670u,        0.4666804339u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1774737670u,        0.4666804339u) withcolor Green;
draw (        0.2210886587u,        0.5674005019u)--(        0.1742126844u,        0.6666448160u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0870979687u,        0.5998779292u)--(        0.1742126844u,        0.6666448160u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2210886587u,        0.5674005019u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2210886587u,        0.5674005019u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0870979687u,        0.5998779292u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0870979687u,        0.5998779292u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1742126844u,        0.6666448160u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1742126844u,        0.6666448160u) withcolor Green;
draw (        0.1814367963u,       -0.1153681548u)--(        0.1575706264u,       -0.0082364750u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1862819111u,        0.0976996048u)--(        0.1575706264u,       -0.0082364750u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1814367963u,       -0.1153681548u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1814367963u,       -0.1153681548u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1862819111u,        0.0976996048u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1862819111u,        0.0976996048u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1575706264u,       -0.0082364750u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1575706264u,       -0.0082364750u) withcolor MidnightBlue;
draw (        0.1814367963u,       -0.1153681548u)--(        0.2101480809u,       -0.0094320750u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1862819111u,        0.0976996048u)--(        0.2101480809u,       -0.0094320750u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1814367963u,       -0.1153681548u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1814367963u,       -0.1153681548u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1862819111u,        0.0976996048u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1862819111u,        0.0976996048u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2101480809u,       -0.0094320750u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2101480809u,       -0.0094320750u) withcolor MidnightBlue;
draw (        0.1814367963u,       -0.1153681548u)--(        0.1575246411u,       -0.0082467297u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0724245090u,       -0.0775630941u)--(        0.1575246411u,       -0.0082467297u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1814367963u,       -0.1153681548u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1814367963u,       -0.1153681548u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0724245090u,       -0.0775630941u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0724245090u,       -0.0775630941u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1575246411u,       -0.0082467297u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1575246411u,       -0.0082467297u) withcolor MidnightBlue;
draw (        0.1862819111u,        0.0976996048u)--(        0.1183202673u,        0.1838854844u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0785464688u,        0.2861832645u)--(        0.1183202673u,        0.1838854844u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1862819111u,        0.0976996048u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1862819111u,        0.0976996048u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0785464688u,        0.2861832645u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0785464688u,        0.2861832645u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1183202673u,        0.1838854844u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1183202673u,        0.1838854844u) withcolor MidnightBlue;
draw (        0.1862819111u,        0.0976996048u)--(        0.1465081125u,        0.1999973849u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0785464688u,        0.2861832645u)--(        0.1465081125u,        0.1999973849u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1862819111u,        0.0976996048u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1862819111u,        0.0976996048u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0785464688u,        0.2861832645u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0785464688u,        0.2861832645u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1465081125u,        0.1999973849u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1465081125u,        0.1999973849u) withcolor Green;
draw (        0.1862819111u,        0.0976996048u)--(        0.1575008582u,       -0.0082175414u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0724245090u,       -0.0775630941u)--(        0.1575008582u,       -0.0082175414u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1862819111u,        0.0976996048u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1862819111u,        0.0976996048u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0724245090u,       -0.0775630941u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0724245090u,       -0.0775630941u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1575008582u,       -0.0082175414u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1575008582u,       -0.0082175414u) withcolor MidnightBlue;
draw (        0.1862819111u,        0.0976996048u)--(        0.1183110679u,        0.1838782295u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0397641253u,        0.1072155353u)--(        0.1183110679u,        0.1838782295u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1862819111u,        0.0976996048u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1862819111u,        0.0976996048u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0397641253u,        0.1072155353u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0397641253u,        0.1072155353u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1183110679u,        0.1838782295u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1183110679u,        0.1838782295u) withcolor MidnightBlue;
draw (        0.0785464688u,        0.2861832645u)--(        0.0000107255u,        0.3628574316u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0693241769u,        0.4479599365u)--(        0.0000107255u,        0.3628574316u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0785464688u,        0.2861832645u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0785464688u,        0.2861832645u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0693241769u,        0.4479599365u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0693241769u,        0.4479599365u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000107255u,        0.3628574316u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000107255u,        0.3628574316u) withcolor Green;
draw (        0.0785464688u,        0.2861832645u)--(        0.1478599201u,        0.3712857693u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0693241769u,        0.4479599365u)--(        0.1478599201u,        0.3712857693u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0785464688u,        0.2861832645u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0785464688u,        0.2861832645u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0693241769u,        0.4479599365u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0693241769u,        0.4479599365u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1478599201u,        0.3712857693u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1478599201u,        0.3712857693u) withcolor MidnightBlue;
draw (        0.0785464688u,        0.2861832645u)--(        0.1183084649u,        0.1838808964u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0397641253u,        0.1072155353u)--(        0.1183084649u,        0.1838808964u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0785464688u,        0.2861832645u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0785464688u,        0.2861832645u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0397641253u,        0.1072155353u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0397641253u,        0.1072155353u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1183084649u,        0.1838808964u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1183084649u,        0.1838808964u) withcolor MidnightBlue;
draw (        0.0785464688u,        0.2861832645u)--(        0.0000021292u,        0.2095179034u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0397641253u,        0.1072155353u)--(        0.0000021292u,        0.2095179034u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0785464688u,        0.2861832645u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0785464688u,        0.2861832645u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0397641253u,        0.1072155353u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0397641253u,        0.1072155353u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000021292u,        0.2095179034u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000021292u,        0.2095179034u) withcolor Green;
draw (        0.0693241769u,        0.4479599365u)--(        0.0000167556u,        0.5330673522u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0870979687u,        0.5998779292u)--(        0.0000167556u,        0.5330673522u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0693241769u,        0.4479599365u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0693241769u,        0.4479599365u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0870979687u,        0.5998779292u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0870979687u,        0.5998779292u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000167556u,        0.5330673522u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000167556u,        0.5330673522u) withcolor Green;
draw (        0.0870979687u,        0.5998779292u)--(        0.0000000000u,        0.6666666667u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0870979687u,        0.5998779292u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0870979687u,        0.5998779292u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.6666666667u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.6666666667u) withcolor Yellow;
draw (        0.0674193742u,       -0.2466921418u)--(        0.0000000000u,       -0.3333333333u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0674193742u,       -0.2466921418u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0674193742u,       -0.2466921418u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,       -0.3333333333u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,       -0.3333333333u) withcolor Yellow;
draw (        0.0674193742u,       -0.2466921418u)--(        0.0000285151u,       -0.1600592209u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0724245090u,       -0.0775630941u)--(        0.0000285151u,       -0.1600592209u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0674193742u,       -0.2466921418u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0674193742u,       -0.2466921418u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0724245090u,       -0.0775630941u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0724245090u,       -0.0775630941u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000285151u,       -0.1600592209u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000285151u,       -0.1600592209u) withcolor Green;
draw (        0.0724245090u,       -0.0775630941u)--(        0.0000047322u,        0.0049121554u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0397641253u,        0.1072155353u)--(        0.0000047322u,        0.0049121554u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0724245090u,       -0.0775630941u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0724245090u,       -0.0775630941u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0397641253u,        0.1072155353u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0397641253u,        0.1072155353u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000047322u,        0.0049121554u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000047322u,        0.0049121554u) withcolor Green;
% r =         0.1097578738
endfig;
end;
