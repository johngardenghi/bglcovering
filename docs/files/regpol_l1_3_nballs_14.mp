prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(14);
u := 20cm;
path pols[];
pair p[];
R := .1275047926u;
D := 2*R;

pols[1] = 
(        0.0000000000u,        0.0000000000u)--
(        1.0000000000u,        0.0000000000u)--
(        0.5000000000u,        0.8660254038u)--
(        0.0000000000u,        0.0000000000u)--
cycle;
fill pols[1] withcolor 0.9Dandelion;
draw (        0.0000000000u,        0.0000000000u)--(        1.0000000000u,        0.0000000000u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        0.0000000000u)--(        0.5000000000u,        0.8660254038u) withpen pencircle scaled 1.0;
draw (        0.5000000000u,        0.8660254038u)--(        0.0000000000u,        0.0000000000u) withpen pencircle scaled 1.0;

x1 = .4322103128u;
y1 = .2495898893u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .4938260362u;
y2 = .4136828564u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .6358910517u;
y3 = .2932093651u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .3369448721u;
y4 = .0501314683u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .3123147804u;
y5 = .4489767206u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .4558041487u;
y6 = .6053246854u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .6617029829u;
y7 = .5253006634u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .2282321533u;
y8 = .2423782350u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .7506338330u;
y9 = .0937508610u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = .5461897674u;
y10 = .7471313741u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;

x11 = .7878944039u;
y11 = .3114798894u;
draw fullcircle scaled D shifted z11 withcolor RoyalBlue withpen pencircle scaled 1.0;

x12 = .5591936015u;
y12 = .0723349192u;
draw fullcircle scaled D shifted z12 withcolor RoyalBlue withpen pencircle scaled 1.0;

x13 = .1098597302u;
y13 = .0647428228u;
draw fullcircle scaled D shifted z13 withcolor RoyalBlue withpen pencircle scaled 1.0;

x14 = .9184965843u;
y14 = .0981093038u;
draw fullcircle scaled D shifted z14 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.4322103128u,        0.2495898893u)--(        0.5494536944u,        0.1994754043u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.6358910517u,        0.2932093651u)--(        0.5494536944u,        0.1994754043u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4322103128u,        0.2495898893u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4322103128u,        0.2495898893u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6358910517u,        0.2932093651u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6358910517u,        0.2932093651u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5494536944u,        0.1994754043u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5494536944u,        0.1994754043u) withcolor MidnightBlue;
draw (        0.4322103128u,        0.2495898893u)--(        0.4419514527u,        0.1224577442u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.3369448721u,        0.0501314683u)--(        0.4419514527u,        0.1224577442u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4322103128u,        0.2495898893u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4322103128u,        0.2495898893u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3369448721u,        0.0501314683u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3369448721u,        0.0501314683u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4419514527u,        0.1224577442u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4419514527u,        0.1224577442u) withcolor MidnightBlue;
draw (        0.4322103128u,        0.2495898893u)--(        0.3275265247u,        0.3223825858u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.3123147804u,        0.4489767206u)--(        0.3275265247u,        0.3223825858u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4322103128u,        0.2495898893u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4322103128u,        0.2495898893u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3123147804u,        0.4489767206u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3123147804u,        0.4489767206u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3275265247u,        0.3223825858u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3275265247u,        0.3223825858u) withcolor MidnightBlue;
draw (        0.4322103128u,        0.2495898893u)--(        0.3275204687u,        0.3223738759u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2282321533u,        0.2423782350u)--(        0.3275204687u,        0.3223738759u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4322103128u,        0.2495898893u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4322103128u,        0.2495898893u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2282321533u,        0.2423782350u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2282321533u,        0.2423782350u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3275204687u,        0.3223738759u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3275204687u,        0.3223738759u) withcolor MidnightBlue;
draw (        0.4322103128u,        0.2495898893u)--(        0.5494500502u,        0.1994668796u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5591936015u,        0.0723349192u)--(        0.5494500502u,        0.1994668796u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4322103128u,        0.2495898893u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4322103128u,        0.2495898893u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5591936015u,        0.0723349192u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5591936015u,        0.0723349192u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5494500502u,        0.1994668796u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5494500502u,        0.1994668796u) withcolor MidnightBlue;
draw (        0.4322103128u,        0.2495898893u)--(        0.4419538640u,        0.1224579290u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5591936015u,        0.0723349192u)--(        0.4419538640u,        0.1224579290u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4322103128u,        0.2495898893u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4322103128u,        0.2495898893u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5591936015u,        0.0723349192u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5591936015u,        0.0723349192u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4419538640u,        0.1224579290u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4419538640u,        0.1224579290u) withcolor MidnightBlue;
draw (        0.6358910517u,        0.2932093651u)--(        0.6997111879u,        0.4035926173u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.6617029829u,        0.5253006634u)--(        0.6997111879u,        0.4035926173u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6358910517u,        0.2932093651u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6358910517u,        0.2932093651u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6617029829u,        0.5253006634u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6617029829u,        0.5253006634u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6997111879u,        0.4035926173u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6997111879u,        0.4035926173u) withcolor MidnightBlue;
draw (        0.6358910517u,        0.2932093651u)--(        0.6456270377u,        0.1660768251u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.7506338330u,        0.0937508610u)--(        0.6456270377u,        0.1660768251u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6358910517u,        0.2932093651u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6358910517u,        0.2932093651u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7506338330u,        0.0937508610u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7506338330u,        0.0937508610u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6456270377u,        0.1660768251u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6456270377u,        0.1660768251u) withcolor MidnightBlue;
draw (        0.6358910517u,        0.2932093651u)--(        0.6997237802u,        0.4035853358u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.7878944039u,        0.3114798894u)--(        0.6997237802u,        0.4035853358u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6358910517u,        0.2932093651u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6358910517u,        0.2932093651u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7878944039u,        0.3114798894u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7878944039u,        0.3114798894u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6997237802u,        0.4035853358u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6997237802u,        0.4035853358u) withcolor MidnightBlue;
draw (        0.6358910517u,        0.2932093651u)--(        0.6456227022u,        0.1660764932u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5591936015u,        0.0723349192u)--(        0.6456227022u,        0.1660764932u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6358910517u,        0.2932093651u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6358910517u,        0.2932093651u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5591936015u,        0.0723349192u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5591936015u,        0.0723349192u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6456227022u,        0.1660764932u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6456227022u,        0.1660764932u) withcolor MidnightBlue;
draw (        0.6358910517u,        0.2932093651u)--(        0.5494619510u,        0.1994677911u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5591936015u,        0.0723349192u)--(        0.5494619510u,        0.1994677911u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6358910517u,        0.2932093651u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6358910517u,        0.2932093651u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5591936015u,        0.0723349192u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5591936015u,        0.0723349192u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5494619510u,        0.1994677911u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5494619510u,        0.1994677911u) withcolor MidnightBlue;
draw (        0.3369448721u,        0.0501314683u)--(        0.2271027073u,        0.1148784448u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2282321533u,        0.2423782350u)--(        0.2271027073u,        0.1148784448u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3369448721u,        0.0501314683u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3369448721u,        0.0501314683u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2282321533u,        0.2423782350u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2282321533u,        0.2423782350u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2271027073u,        0.1148784448u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2271027073u,        0.1148784448u) withcolor MidnightBlue;
draw (        0.3369448721u,        0.0501314683u)--(        0.4419528918u,        0.1224556548u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5591936015u,        0.0723349192u)--(        0.4419528918u,        0.1224556548u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3369448721u,        0.0501314683u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3369448721u,        0.0501314683u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5591936015u,        0.0723349192u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5591936015u,        0.0723349192u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4419528918u,        0.1224556548u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4419528918u,        0.1224556548u) withcolor MidnightBlue;
draw (        0.3369448721u,        0.0501314683u)--(        0.4541855818u,        0.0000107327u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5591936015u,        0.0723349192u)--(        0.4541855818u,        0.0000107327u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3369448721u,        0.0501314683u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3369448721u,        0.0501314683u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5591936015u,        0.0723349192u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5591936015u,        0.0723349192u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4541855818u,        0.0000107327u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4541855818u,        0.0000107327u) withcolor Green;
draw (        0.3369448721u,        0.0501314683u)--(        0.2197068980u,        0.0000043343u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.1098597302u,        0.0647428228u)--(        0.2197068980u,        0.0000043343u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3369448721u,        0.0501314683u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3369448721u,        0.0501314683u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1098597302u,        0.0647428228u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1098597302u,        0.0647428228u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2197068980u,        0.0000043343u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2197068980u,        0.0000043343u) withcolor Green;
draw (        0.3369448721u,        0.0501314683u)--(        0.2270977044u,        0.1148699568u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.1098597302u,        0.0647428228u)--(        0.2270977044u,        0.1148699568u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3369448721u,        0.0501314683u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3369448721u,        0.0501314683u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1098597302u,        0.0647428228u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1098597302u,        0.0647428228u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2270977044u,        0.1148699568u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2270977044u,        0.1148699568u) withcolor MidnightBlue;
draw (        0.3123147804u,        0.4489767206u)--(        0.3319682920u,        0.5749577177u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.4558041487u,        0.6053246854u)--(        0.3319682920u,        0.5749577177u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3123147804u,        0.4489767206u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3123147804u,        0.4489767206u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4558041487u,        0.6053246854u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4558041487u,        0.6053246854u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3319682920u,        0.5749577177u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3319682920u,        0.5749577177u) withcolor Green;
draw (        0.3123147804u,        0.4489767206u)--(        0.3275146044u,        0.3223811540u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2282321533u,        0.2423782350u)--(        0.3275146044u,        0.3223811540u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3123147804u,        0.4489767206u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3123147804u,        0.4489767206u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2282321533u,        0.2423782350u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2282321533u,        0.2423782350u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3275146044u,        0.3223811540u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3275146044u,        0.3223811540u) withcolor MidnightBlue;
draw (        0.3123147804u,        0.4489767206u)--(        0.2130323293u,        0.3689738017u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2282321533u,        0.2423782350u)--(        0.2130323293u,        0.3689738017u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3123147804u,        0.4489767206u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3123147804u,        0.4489767206u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2282321533u,        0.2423782350u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2282321533u,        0.2423782350u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2130323293u,        0.3689738017u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2130323293u,        0.3689738017u) withcolor Green;
draw (        0.4558041487u,        0.6053246854u)--(        0.5818301111u,        0.6246877677u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.6617029829u,        0.5253006634u)--(        0.5818301111u,        0.6246877677u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4558041487u,        0.6053246854u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4558041487u,        0.6053246854u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6617029829u,        0.5253006634u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6617029829u,        0.5253006634u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5818301111u,        0.6246877677u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5818301111u,        0.6246877677u) withcolor MidnightBlue;
draw (        0.4558041487u,        0.6053246854u)--(        0.4201669289u,        0.7277479704u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5461897674u,        0.7471313741u)--(        0.4201669289u,        0.7277479704u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4558041487u,        0.6053246854u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4558041487u,        0.6053246854u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5461897674u,        0.7471313741u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5461897674u,        0.7471313741u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4201669289u,        0.7277479704u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4201669289u,        0.7277479704u) withcolor Green;
draw (        0.4558041487u,        0.6053246854u)--(        0.5818269872u,        0.6247080891u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5461897674u,        0.7471313741u)--(        0.5818269872u,        0.6247080891u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4558041487u,        0.6053246854u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4558041487u,        0.6053246854u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5461897674u,        0.7471313741u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5461897674u,        0.7471313741u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5818269872u,        0.6247080891u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5818269872u,        0.6247080891u) withcolor MidnightBlue;
draw (        0.6617029829u,        0.5253006634u)--(        0.5818715775u,        0.6247210781u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5461897674u,        0.7471313741u)--(        0.5818715775u,        0.6247210781u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6617029829u,        0.5253006634u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6617029829u,        0.5253006634u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5461897674u,        0.7471313741u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5461897674u,        0.7471313741u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5818715775u,        0.6247210781u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5818715775u,        0.6247210781u) withcolor MidnightBlue;
draw (        0.6617029829u,        0.5253006634u)--(        0.6260211728u,        0.6477109595u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5461897674u,        0.7471313741u)--(        0.6260211728u,        0.6477109595u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6617029829u,        0.5253006634u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6617029829u,        0.5253006634u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5461897674u,        0.7471313741u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5461897674u,        0.7471313741u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6260211728u,        0.6477109595u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6260211728u,        0.6477109595u) withcolor Green;
draw (        0.6617029829u,        0.5253006634u)--(        0.7498562092u,        0.4331785660u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.7878944039u,        0.3114798894u)--(        0.7498562092u,        0.4331785660u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6617029829u,        0.5253006634u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6617029829u,        0.5253006634u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7878944039u,        0.3114798894u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7878944039u,        0.3114798894u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7498562092u,        0.4331785660u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7498562092u,        0.4331785660u) withcolor Green;
draw (        0.6617029829u,        0.5253006634u)--(        0.6997411776u,        0.4036019868u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.7878944039u,        0.3114798894u)--(        0.6997411776u,        0.4036019868u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6617029829u,        0.5253006634u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6617029829u,        0.5253006634u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7878944039u,        0.3114798894u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7878944039u,        0.3114798894u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6997411776u,        0.4036019868u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6997411776u,        0.4036019868u) withcolor MidnightBlue;
draw (        0.2282321533u,        0.2423782350u)--(        0.2270940419u,        0.1148785219u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.1098597302u,        0.0647428228u)--(        0.2270940419u,        0.1148785219u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2282321533u,        0.2423782350u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2282321533u,        0.2423782350u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1098597302u,        0.0647428228u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1098597302u,        0.0647428228u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2270940419u,        0.1148785219u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2270940419u,        0.1148785219u) withcolor MidnightBlue;
draw (        0.2282321533u,        0.2423782350u)--(        0.1109978416u,        0.1922425360u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.1098597302u,        0.0647428228u)--(        0.1109978416u,        0.1922425360u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2282321533u,        0.2423782350u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2282321533u,        0.2423782350u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1098597302u,        0.0647428228u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1098597302u,        0.0647428228u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1109978416u,        0.1922425360u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1109978416u,        0.1922425360u) withcolor Green;
draw (        0.7506338330u,        0.0937508610u)--(        0.8320607241u,        0.1918688184u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.7878944039u,        0.3114798894u)--(        0.8320607241u,        0.1918688184u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7506338330u,        0.0937508610u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7506338330u,        0.0937508610u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7878944039u,        0.3114798894u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7878944039u,        0.3114798894u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8320607241u,        0.1918688184u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8320607241u,        0.1918688184u) withcolor MidnightBlue;
draw (        0.7506338330u,        0.0937508610u)--(        0.6642022205u,        0.0000116029u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5591936015u,        0.0723349192u)--(        0.6642022205u,        0.0000116029u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7506338330u,        0.0937508610u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7506338330u,        0.0937508610u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5591936015u,        0.0723349192u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5591936015u,        0.0723349192u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6642022205u,        0.0000116029u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6642022205u,        0.0000116029u) withcolor Green;
draw (        0.7506338330u,        0.0937508610u)--(        0.6456252140u,        0.1660741773u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5591936015u,        0.0723349192u)--(        0.6456252140u,        0.1660741773u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7506338330u,        0.0937508610u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7506338330u,        0.0937508610u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5591936015u,        0.0723349192u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5591936015u,        0.0723349192u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6456252140u,        0.1660741773u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6456252140u,        0.1660741773u) withcolor MidnightBlue;
draw (        0.7506338330u,        0.0937508610u)--(        0.8320745216u,        0.1918573663u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.9184965843u,        0.0981093038u)--(        0.8320745216u,        0.1918573663u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7506338330u,        0.0937508610u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7506338330u,        0.0937508610u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9184965843u,        0.0981093038u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9184965843u,        0.0981093038u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8320745216u,        0.1918573663u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8320745216u,        0.1918573663u) withcolor MidnightBlue;
draw (        0.7506338330u,        0.0937508610u)--(        0.8370558957u,        0.0000027985u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.9184965843u,        0.0981093038u)--(        0.8370558957u,        0.0000027985u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7506338330u,        0.0937508610u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7506338330u,        0.0937508610u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9184965843u,        0.0981093038u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9184965843u,        0.0981093038u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8370558957u,        0.0000027985u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8370558957u,        0.0000027985u) withcolor Green;
draw (        0.5461897674u,        0.7471313741u)--(        0.5000000000u,        0.8660254038u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5461897674u,        0.7471313741u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5461897674u,        0.7471313741u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.8660254038u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.8660254038u) withcolor Yellow;
draw (        0.7878944039u,        0.3114798894u)--(        0.8742864908u,        0.2177042024u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.9184965843u,        0.0981093038u)--(        0.8742864908u,        0.2177042024u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7878944039u,        0.3114798894u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7878944039u,        0.3114798894u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9184965843u,        0.0981093038u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9184965843u,        0.0981093038u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8742864908u,        0.2177042024u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8742864908u,        0.2177042024u) withcolor Green;
draw (        0.7878944039u,        0.3114798894u)--(        0.8321044974u,        0.1918849908u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.9184965843u,        0.0981093038u)--(        0.8321044974u,        0.1918849908u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7878944039u,        0.3114798894u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7878944039u,        0.3114798894u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9184965843u,        0.0981093038u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9184965843u,        0.0981093038u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8321044974u,        0.1918849908u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8321044974u,        0.1918849908u) withcolor MidnightBlue;
draw (        0.1098597302u,        0.0647428228u)--(        0.0000000000u,        0.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1098597302u,        0.0647428228u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1098597302u,        0.0647428228u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.9184965843u,        0.0981093038u)--(        1.0000000000u,        0.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9184965843u,        0.0981093038u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9184965843u,        0.0981093038u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor Yellow;
% r =         0.1275047926
endfig;
end;
