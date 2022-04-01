prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(20);
u := 20cm;
path pols[];
pair p[];
R := .1332687803u;
D := 2*R;

pols[1] = 
(        0.0000000000u,        0.0000000000u)--
(        0.1111111111u,        0.0000000000u)--
(        0.1666666667u,        0.0962250449u)--
(        0.0962250449u,        0.1666666667u)--
(        0.0000000000u,        0.1111111111u)--
(        0.0000000000u,        0.0000000000u)--
cycle;
pols[2] = 
(        0.1666666667u,        0.0962250449u)--
(        0.2222222222u,        0.0000000000u)--
(        0.3333333333u,        0.0000000000u)--
(        0.3888888889u,        0.0962250449u)--
(        0.3333333333u,        0.1924500897u)--
(        0.1666666667u,        0.0962250449u)--
cycle;
pols[3] = 
(        0.1924500897u,        0.3333333333u)--
(        0.0962250449u,        0.3888888889u)--
(        0.0000000000u,        0.3333333333u)--
(        0.0000000000u,        0.2222222222u)--
(        0.0962250449u,        0.1666666667u)--
(        0.1924500897u,        0.3333333333u)--
cycle;
pols[4] = 
(        0.0962250449u,        0.1666666667u)--
(        0.1666666667u,        0.0962250449u)--
(        0.3333333333u,        0.1924500897u)--
(        0.1924500897u,        0.3333333333u)--
(        0.0962250449u,        0.1666666667u)--
cycle;
pols[5] = 
(        0.1924500897u,        0.3333333333u)--
(        0.3333333333u,        0.1924500897u)--
(        0.4444444444u,        0.1924500897u)--
(        0.5000000000u,        0.2886751346u)--
(        0.2886751346u,        0.5000000000u)--
(        0.1924500897u,        0.4444444444u)--
(        0.1924500897u,        0.3333333333u)--
cycle;
pols[6] = 
(        0.6666666667u,        0.1924500897u)--
(        0.6111111111u,        0.0962250449u)--
(        0.6666666667u,        0.0000000000u)--
(        0.7777777778u,        0.0000000000u)--
(        0.8333333333u,        0.0962250449u)--
(        0.6666666667u,        0.1924500897u)--
cycle;
pols[7] = 
(        0.8333333333u,        0.0962250449u)--
(        0.8888888889u,        0.0000000000u)--
(        1.0000000000u,        0.0000000000u)--
(        1.0000000000u,        0.1111111111u)--
(        0.9037749551u,        0.1666666667u)--
(        0.8333333333u,        0.0962250449u)--
cycle;
pols[8] = 
(        0.9037749551u,        0.1666666667u)--
(        1.0000000000u,        0.2222222222u)--
(        1.0000000000u,        0.3333333333u)--
(        0.9037749551u,        0.3888888889u)--
(        0.8075499103u,        0.3333333333u)--
(        0.9037749551u,        0.1666666667u)--
cycle;
pols[9] = 
(        0.9037749551u,        0.1666666667u)--
(        0.8075499103u,        0.3333333333u)--
(        0.6666666667u,        0.1924500897u)--
(        0.8333333333u,        0.0962250449u)--
(        0.9037749551u,        0.1666666667u)--
cycle;
pols[10] = 
(        0.5000000000u,        0.2886751346u)--
(        0.5555555556u,        0.1924500897u)--
(        0.6666666667u,        0.1924500897u)--
(        0.8075499103u,        0.3333333333u)--
(        0.8075499103u,        0.4444444444u)--
(        0.7113248654u,        0.5000000000u)--
(        0.5000000000u,        0.2886751346u)--
cycle;
pols[11] = 
(        0.8075499103u,        0.6666666667u)--
(        0.9037749551u,        0.6111111111u)--
(        1.0000000000u,        0.6666666667u)--
(        1.0000000000u,        0.7777777778u)--
(        0.9037749551u,        0.8333333333u)--
(        0.8075499103u,        0.6666666667u)--
cycle;
pols[12] = 
(        0.9037749551u,        0.8333333333u)--
(        1.0000000000u,        0.8888888889u)--
(        1.0000000000u,        1.0000000000u)--
(        0.8888888889u,        1.0000000000u)--
(        0.8333333333u,        0.9037749551u)--
(        0.9037749551u,        0.8333333333u)--
cycle;
pols[13] = 
(        0.8333333333u,        0.9037749551u)--
(        0.7777777778u,        1.0000000000u)--
(        0.6666666667u,        1.0000000000u)--
(        0.6111111111u,        0.9037749551u)--
(        0.6666666667u,        0.8075499103u)--
(        0.8333333333u,        0.9037749551u)--
cycle;
pols[14] = 
(        0.8333333333u,        0.9037749551u)--
(        0.6666666667u,        0.8075499103u)--
(        0.8075499103u,        0.6666666667u)--
(        0.9037749551u,        0.8333333333u)--
(        0.8333333333u,        0.9037749551u)--
cycle;
pols[15] = 
(        0.8075499103u,        0.6666666667u)--
(        0.6666666667u,        0.8075499103u)--
(        0.5555555556u,        0.8075499103u)--
(        0.5000000000u,        0.7113248654u)--
(        0.7113248654u,        0.5000000000u)--
(        0.8075499103u,        0.5555555556u)--
(        0.8075499103u,        0.6666666667u)--
cycle;
pols[16] = 
(        0.3333333333u,        0.8075499103u)--
(        0.3888888889u,        0.9037749551u)--
(        0.3333333333u,        1.0000000000u)--
(        0.2222222222u,        1.0000000000u)--
(        0.1666666667u,        0.9037749551u)--
(        0.3333333333u,        0.8075499103u)--
cycle;
pols[17] = 
(        0.1666666667u,        0.9037749551u)--
(        0.1111111111u,        1.0000000000u)--
(        0.0000000000u,        1.0000000000u)--
(        0.0000000000u,        0.8888888889u)--
(        0.0962250449u,        0.8333333333u)--
(        0.1666666667u,        0.9037749551u)--
cycle;
pols[18] = 
(        0.0962250449u,        0.8333333333u)--
(        0.0000000000u,        0.7777777778u)--
(        0.0000000000u,        0.6666666667u)--
(        0.0962250449u,        0.6111111111u)--
(        0.1924500897u,        0.6666666667u)--
(        0.0962250449u,        0.8333333333u)--
cycle;
pols[19] = 
(        0.0962250449u,        0.8333333333u)--
(        0.1924500897u,        0.6666666667u)--
(        0.3333333333u,        0.8075499103u)--
(        0.1666666667u,        0.9037749551u)--
(        0.0962250449u,        0.8333333333u)--
cycle;
pols[20] = 
(        0.3333333333u,        0.8075499103u)--
(        0.1924500897u,        0.6666666667u)--
(        0.1924500897u,        0.5555555556u)--
(        0.2886751346u,        0.5000000000u)--
(        0.5000000000u,        0.7113248654u)--
(        0.4444444444u,        0.8075499103u)--
(        0.3333333333u,        0.8075499103u)--
cycle;
pols[21] = 
(        0.5000000000u,        0.7113248654u)--
(        0.2886751346u,        0.5000000000u)--
(        0.5000000000u,        0.2886751346u)--
(        0.7113248654u,        0.5000000000u)--
(        0.5000000000u,        0.7113248654u)--
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
draw (        0.1666666667u,        0.0962250449u)--(        0.0962250449u,        0.1666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3333333333u,        0.1924500897u)--(        0.1666666667u,        0.0962250449u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.0962250449u,        0.1666666667u)--(        0.1924500897u,        0.3333333333u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.0962250449u,        0.1666666667u)--(        0.1666666667u,        0.0962250449u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1666666667u,        0.0962250449u)--(        0.3333333333u,        0.1924500897u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3333333333u,        0.1924500897u)--(        0.1924500897u,        0.3333333333u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1924500897u,        0.3333333333u)--(        0.0962250449u,        0.1666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1924500897u,        0.3333333333u)--(        0.3333333333u,        0.1924500897u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5000000000u,        0.2886751346u)--(        0.2886751346u,        0.5000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.8333333333u,        0.0962250449u)--(        0.6666666667u,        0.1924500897u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.9037749551u,        0.1666666667u)--(        0.8333333333u,        0.0962250449u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.8075499103u,        0.3333333333u)--(        0.9037749551u,        0.1666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.9037749551u,        0.1666666667u)--(        0.8075499103u,        0.3333333333u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.8075499103u,        0.3333333333u)--(        0.6666666667u,        0.1924500897u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.6666666667u,        0.1924500897u)--(        0.8333333333u,        0.0962250449u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.8333333333u,        0.0962250449u)--(        0.9037749551u,        0.1666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.6666666667u,        0.1924500897u)--(        0.8075499103u,        0.3333333333u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.7113248654u,        0.5000000000u)--(        0.5000000000u,        0.2886751346u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.9037749551u,        0.8333333333u)--(        0.8075499103u,        0.6666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.8333333333u,        0.9037749551u)--(        0.9037749551u,        0.8333333333u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.6666666667u,        0.8075499103u)--(        0.8333333333u,        0.9037749551u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.8333333333u,        0.9037749551u)--(        0.6666666667u,        0.8075499103u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.6666666667u,        0.8075499103u)--(        0.8075499103u,        0.6666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.8075499103u,        0.6666666667u)--(        0.9037749551u,        0.8333333333u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.9037749551u,        0.8333333333u)--(        0.8333333333u,        0.9037749551u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.8075499103u,        0.6666666667u)--(        0.6666666667u,        0.8075499103u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5000000000u,        0.7113248654u)--(        0.7113248654u,        0.5000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1666666667u,        0.9037749551u)--(        0.3333333333u,        0.8075499103u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.0962250449u,        0.8333333333u)--(        0.1666666667u,        0.9037749551u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1924500897u,        0.6666666667u)--(        0.0962250449u,        0.8333333333u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.0962250449u,        0.8333333333u)--(        0.1924500897u,        0.6666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1924500897u,        0.6666666667u)--(        0.3333333333u,        0.8075499103u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3333333333u,        0.8075499103u)--(        0.1666666667u,        0.9037749551u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.1666666667u,        0.9037749551u)--(        0.0962250449u,        0.8333333333u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.3333333333u,        0.8075499103u)--(        0.1924500897u,        0.6666666667u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2886751346u,        0.5000000000u)--(        0.5000000000u,        0.7113248654u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5000000000u,        0.7113248654u)--(        0.2886751346u,        0.5000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.2886751346u,        0.5000000000u)--(        0.5000000000u,        0.2886751346u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.5000000000u,        0.2886751346u)--(        0.7113248654u,        0.5000000000u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.7113248654u,        0.5000000000u)--(        0.5000000000u,        0.7113248654u) withpen pencircle scaled 1.0 withcolor 0.9Dandelion;
draw (        0.0000000000u,        0.0000000000u)--(        0.1111111111u,        0.0000000000u) withpen pencircle scaled 1.0;
draw (        0.1111111111u,        0.0000000000u)--(        0.1666666667u,        0.0962250449u) withpen pencircle scaled 1.0;
draw (        0.0962250449u,        0.1666666667u)--(        0.0000000000u,        0.1111111111u) withpen pencircle scaled 1.0;
draw (        0.0000000000u,        0.1111111111u)--(        0.0000000000u,        0.0000000000u) withpen pencircle scaled 1.0;
draw (        0.1666666667u,        0.0962250449u)--(        0.2222222222u,        0.0000000000u) withpen pencircle scaled 1.0;
draw (        0.2222222222u,        0.0000000000u)--(        0.3333333333u,        0.0000000000u) withpen pencircle scaled 1.0;
draw (        0.3333333333u,        0.0000000000u)--(        0.3888888889u,        0.0962250449u) withpen pencircle scaled 1.0;
draw (        0.3888888889u,        0.0962250449u)--(        0.3333333333u,        0.1924500897u) withpen pencircle scaled 1.0;
draw (        0.1924500897u,        0.3333333333u)--(        0.0962250449u,        0.3888888889u) withpen pencircle scaled 1.0;
draw (        0.0962250449u,        0.3888888889u)--(        0.0000000000u,        0.3333333333u) withpen pencircle scaled 1.0;
draw (        0.0000000000u,        0.3333333333u)--(        0.0000000000u,        0.2222222222u) withpen pencircle scaled 1.0;
draw (        0.0000000000u,        0.2222222222u)--(        0.0962250449u,        0.1666666667u) withpen pencircle scaled 1.0;
draw (        0.3333333333u,        0.1924500897u)--(        0.4444444444u,        0.1924500897u) withpen pencircle scaled 1.0;
draw (        0.4444444444u,        0.1924500897u)--(        0.5000000000u,        0.2886751346u) withpen pencircle scaled 1.0;
draw (        0.2886751346u,        0.5000000000u)--(        0.1924500897u,        0.4444444444u) withpen pencircle scaled 1.0;
draw (        0.1924500897u,        0.4444444444u)--(        0.1924500897u,        0.3333333333u) withpen pencircle scaled 1.0;
draw (        0.6666666667u,        0.1924500897u)--(        0.6111111111u,        0.0962250449u) withpen pencircle scaled 1.0;
draw (        0.6111111111u,        0.0962250449u)--(        0.6666666667u,        0.0000000000u) withpen pencircle scaled 1.0;
draw (        0.6666666667u,        0.0000000000u)--(        0.7777777778u,        0.0000000000u) withpen pencircle scaled 1.0;
draw (        0.7777777778u,        0.0000000000u)--(        0.8333333333u,        0.0962250449u) withpen pencircle scaled 1.0;
draw (        0.8333333333u,        0.0962250449u)--(        0.8888888889u,        0.0000000000u) withpen pencircle scaled 1.0;
draw (        0.8888888889u,        0.0000000000u)--(        1.0000000000u,        0.0000000000u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        0.0000000000u)--(        1.0000000000u,        0.1111111111u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        0.1111111111u)--(        0.9037749551u,        0.1666666667u) withpen pencircle scaled 1.0;
draw (        0.9037749551u,        0.1666666667u)--(        1.0000000000u,        0.2222222222u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        0.2222222222u)--(        1.0000000000u,        0.3333333333u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        0.3333333333u)--(        0.9037749551u,        0.3888888889u) withpen pencircle scaled 1.0;
draw (        0.9037749551u,        0.3888888889u)--(        0.8075499103u,        0.3333333333u) withpen pencircle scaled 1.0;
draw (        0.5000000000u,        0.2886751346u)--(        0.5555555556u,        0.1924500897u) withpen pencircle scaled 1.0;
draw (        0.5555555556u,        0.1924500897u)--(        0.6666666667u,        0.1924500897u) withpen pencircle scaled 1.0;
draw (        0.8075499103u,        0.3333333333u)--(        0.8075499103u,        0.4444444444u) withpen pencircle scaled 1.0;
draw (        0.8075499103u,        0.4444444444u)--(        0.7113248654u,        0.5000000000u) withpen pencircle scaled 1.0;
draw (        0.8075499103u,        0.6666666667u)--(        0.9037749551u,        0.6111111111u) withpen pencircle scaled 1.0;
draw (        0.9037749551u,        0.6111111111u)--(        1.0000000000u,        0.6666666667u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        0.6666666667u)--(        1.0000000000u,        0.7777777778u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        0.7777777778u)--(        0.9037749551u,        0.8333333333u) withpen pencircle scaled 1.0;
draw (        0.9037749551u,        0.8333333333u)--(        1.0000000000u,        0.8888888889u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        0.8888888889u)--(        1.0000000000u,        1.0000000000u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        1.0000000000u)--(        0.8888888889u,        1.0000000000u) withpen pencircle scaled 1.0;
draw (        0.8888888889u,        1.0000000000u)--(        0.8333333333u,        0.9037749551u) withpen pencircle scaled 1.0;
draw (        0.8333333333u,        0.9037749551u)--(        0.7777777778u,        1.0000000000u) withpen pencircle scaled 1.0;
draw (        0.7777777778u,        1.0000000000u)--(        0.6666666667u,        1.0000000000u) withpen pencircle scaled 1.0;
draw (        0.6666666667u,        1.0000000000u)--(        0.6111111111u,        0.9037749551u) withpen pencircle scaled 1.0;
draw (        0.6111111111u,        0.9037749551u)--(        0.6666666667u,        0.8075499103u) withpen pencircle scaled 1.0;
draw (        0.6666666667u,        0.8075499103u)--(        0.5555555556u,        0.8075499103u) withpen pencircle scaled 1.0;
draw (        0.5555555556u,        0.8075499103u)--(        0.5000000000u,        0.7113248654u) withpen pencircle scaled 1.0;
draw (        0.7113248654u,        0.5000000000u)--(        0.8075499103u,        0.5555555556u) withpen pencircle scaled 1.0;
draw (        0.8075499103u,        0.5555555556u)--(        0.8075499103u,        0.6666666667u) withpen pencircle scaled 1.0;
draw (        0.3333333333u,        0.8075499103u)--(        0.3888888889u,        0.9037749551u) withpen pencircle scaled 1.0;
draw (        0.3888888889u,        0.9037749551u)--(        0.3333333333u,        1.0000000000u) withpen pencircle scaled 1.0;
draw (        0.3333333333u,        1.0000000000u)--(        0.2222222222u,        1.0000000000u) withpen pencircle scaled 1.0;
draw (        0.2222222222u,        1.0000000000u)--(        0.1666666667u,        0.9037749551u) withpen pencircle scaled 1.0;
draw (        0.1666666667u,        0.9037749551u)--(        0.1111111111u,        1.0000000000u) withpen pencircle scaled 1.0;
draw (        0.1111111111u,        1.0000000000u)--(        0.0000000000u,        1.0000000000u) withpen pencircle scaled 1.0;
draw (        0.0000000000u,        1.0000000000u)--(        0.0000000000u,        0.8888888889u) withpen pencircle scaled 1.0;
draw (        0.0000000000u,        0.8888888889u)--(        0.0962250449u,        0.8333333333u) withpen pencircle scaled 1.0;
draw (        0.0962250449u,        0.8333333333u)--(        0.0000000000u,        0.7777777778u) withpen pencircle scaled 1.0;
draw (        0.0000000000u,        0.7777777778u)--(        0.0000000000u,        0.6666666667u) withpen pencircle scaled 1.0;
draw (        0.0000000000u,        0.6666666667u)--(        0.0962250449u,        0.6111111111u) withpen pencircle scaled 1.0;
draw (        0.0962250449u,        0.6111111111u)--(        0.1924500897u,        0.6666666667u) withpen pencircle scaled 1.0;
draw (        0.1924500897u,        0.6666666667u)--(        0.1924500897u,        0.5555555556u) withpen pencircle scaled 1.0;
draw (        0.1924500897u,        0.5555555556u)--(        0.2886751346u,        0.5000000000u) withpen pencircle scaled 1.0;
draw (        0.5000000000u,        0.7113248654u)--(        0.4444444444u,        0.8075499103u) withpen pencircle scaled 1.0;
draw (        0.4444444444u,        0.8075499103u)--(        0.3333333333u,        0.8075499103u) withpen pencircle scaled 1.0;

x1 = .9444503584u;
y1 = .8788603887u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .8836989872u;
y2 = .6926105155u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .7222223225u;
y3 = .8788630393u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .9444503584u;
y4 = .1211396113u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .8836989871u;
y5 = .3073894845u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .7252162935u;
y6 = .5000000001u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .6491291128u;
y7 = .6861018454u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .7222223226u;
y8 = .1211369606u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .6491291128u;
y9 = .3138981547u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = .4752048465u;
y10 = .2523993539u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;

x11 = .4871511679u;
y11 = .4999999999u;
draw fullcircle scaled D shifted z11 withcolor RoyalBlue withpen pencircle scaled 1.0;

x12 = .4752048465u;
y12 = .7476006462u;
draw fullcircle scaled D shifted z12 withcolor RoyalBlue withpen pencircle scaled 1.0;

x13 = .3064643645u;
y13 = .3754369250u;
draw fullcircle scaled D shifted z13 withcolor RoyalBlue withpen pencircle scaled 1.0;

x14 = .3064643646u;
y14 = .6245630749u;
draw fullcircle scaled D shifted z14 withcolor RoyalBlue withpen pencircle scaled 1.0;

x15 = .2878089772u;
y15 = .8747404399u;
draw fullcircle scaled D shifted z15 withcolor RoyalBlue withpen pencircle scaled 1.0;

x16 = .1211383817u;
y16 = .9444476789u;
draw fullcircle scaled D shifted z16 withcolor RoyalBlue withpen pencircle scaled 1.0;

x17 = .2878089771u;
y17 = .1252595602u;
draw fullcircle scaled D shifted z17 withcolor RoyalBlue withpen pencircle scaled 1.0;

x18 = .1211404610u;
y18 = .7222253340u;
draw fullcircle scaled D shifted z18 withcolor RoyalBlue withpen pencircle scaled 1.0;

x19 = .1211383816u;
y19 = .0555523212u;
draw fullcircle scaled D shifted z19 withcolor RoyalBlue withpen pencircle scaled 1.0;

x20 = .1211404611u;
y20 = .2777746659u;
draw fullcircle scaled D shifted z20 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.9444503584u,        0.8788603887u)--(        1.0000000000u,        1.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9444503584u,        0.8788603887u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9444503584u,        0.8788603887u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        1.0000000000u) withcolor Yellow;
draw (        0.9444503584u,        0.8788603887u)--(        0.8888888889u,        1.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9444503584u,        0.8788603887u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9444503584u,        0.8788603887u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8888888889u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8888888889u,        1.0000000000u) withcolor Yellow;
draw (        0.8836989872u,        0.6926105155u)--(        0.7646613099u,        0.7525321240u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7222223225u,        0.8788630393u)--(        0.7646613099u,        0.7525321240u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8836989872u,        0.6926105155u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8836989872u,        0.6926105155u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7222223225u,        0.8788630393u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7222223225u,        0.8788630393u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7646613099u,        0.7525321240u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7646613099u,        0.7525321240u) withcolor MidnightBlue;
draw (        0.8836989872u,        0.6926105155u)--(        0.7681832994u,        0.6261523056u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7252162935u,        0.5000000001u)--(        0.7681832994u,        0.6261523056u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8836989872u,        0.6926105155u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8836989872u,        0.6926105155u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7252162935u,        0.5000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7252162935u,        0.5000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7681832994u,        0.6261523056u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7681832994u,        0.6261523056u) withcolor MidnightBlue;
draw (        0.8836989872u,        0.6926105155u)--(        0.7681669976u,        0.6261806490u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6491291128u,        0.6861018454u)--(        0.7681669976u,        0.6261806490u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8836989872u,        0.6926105155u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8836989872u,        0.6926105155u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6491291128u,        0.6861018454u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6491291128u,        0.6861018454u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7681669976u,        0.6261806490u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7681669976u,        0.6261806490u) withcolor MidnightBlue;
draw (        0.8836989872u,        0.6926105155u)--(        0.7646611024u,        0.7525317120u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6491291128u,        0.6861018454u)--(        0.7646611024u,        0.7525317120u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8836989872u,        0.6926105155u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8836989872u,        0.6926105155u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6491291128u,        0.6861018454u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6491291128u,        0.6861018454u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7646611024u,        0.7525317120u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7646611024u,        0.7525317120u) withcolor MidnightBlue;
draw (        0.7222223225u,        0.8788630393u)--(        0.7777777778u,        1.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7222223225u,        0.8788630393u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7222223225u,        0.8788630393u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7777777778u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7777777778u,        1.0000000000u) withcolor Yellow;
draw (        0.7222223225u,        0.8788630393u)--(        0.6666666667u,        1.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7222223225u,        0.8788630393u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7222223225u,        0.8788630393u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6666666667u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6666666667u,        1.0000000000u) withcolor Yellow;
draw (        0.7222223225u,        0.8788630393u)--(        0.7646609374u,        0.7525319989u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6491291128u,        0.6861018454u)--(        0.7646609374u,        0.7525319989u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7222223225u,        0.8788630393u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7222223225u,        0.8788630393u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6491291128u,        0.6861018454u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6491291128u,        0.6861018454u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7646609374u,        0.7525319989u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7646609374u,        0.7525319989u) withcolor MidnightBlue;
draw (        0.9444503584u,        0.1211396113u)--(        0.8888888889u,        0.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9444503584u,        0.1211396113u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9444503584u,        0.1211396113u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8888888889u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8888888889u,        0.0000000000u) withcolor Yellow;
draw (        0.9444503584u,        0.1211396113u)--(        1.0000000000u,        0.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9444503584u,        0.1211396113u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9444503584u,        0.1211396113u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.8836989871u,        0.3073894845u)--(        0.7681832994u,        0.3738476945u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7252162935u,        0.5000000001u)--(        0.7681832994u,        0.3738476945u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8836989871u,        0.3073894845u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8836989871u,        0.3073894845u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7252162935u,        0.5000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7252162935u,        0.5000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7681832994u,        0.3738476945u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7681832994u,        0.3738476945u) withcolor MidnightBlue;
draw (        0.8836989871u,        0.3073894845u)--(        0.7646613098u,        0.2474678759u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7222223226u,        0.1211369606u)--(        0.7646613098u,        0.2474678759u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8836989871u,        0.3073894845u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8836989871u,        0.3073894845u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7222223226u,        0.1211369606u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7222223226u,        0.1211369606u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7646613098u,        0.2474678759u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7646613098u,        0.2474678759u) withcolor MidnightBlue;
draw (        0.8836989871u,        0.3073894845u)--(        0.7646611024u,        0.2474682881u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6491291128u,        0.3138981547u)--(        0.7646611024u,        0.2474682881u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8836989871u,        0.3073894845u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8836989871u,        0.3073894845u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6491291128u,        0.3138981547u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6491291128u,        0.3138981547u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7646611024u,        0.2474682881u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7646611024u,        0.2474682881u) withcolor MidnightBlue;
draw (        0.8836989871u,        0.3073894845u)--(        0.7681669975u,        0.3738193512u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6491291128u,        0.3138981547u)--(        0.7681669975u,        0.3738193512u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8836989871u,        0.3073894845u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8836989871u,        0.3073894845u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6491291128u,        0.3138981547u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6491291128u,        0.3138981547u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7681669975u,        0.3738193512u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7681669975u,        0.3738193512u) withcolor MidnightBlue;
draw (        0.7252162935u,        0.5000000001u)--(        0.6061882051u,        0.5599406539u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6491291128u,        0.6861018454u)--(        0.6061882051u,        0.5599406539u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7252162935u,        0.5000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7252162935u,        0.5000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6491291128u,        0.6861018454u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6491291128u,        0.6861018454u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6061882051u,        0.5599406539u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6061882051u,        0.5599406539u) withcolor MidnightBlue;
draw (        0.7252162935u,        0.5000000001u)--(        0.7681572011u,        0.6261611916u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6491291128u,        0.6861018454u)--(        0.7681572011u,        0.6261611916u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7252162935u,        0.5000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7252162935u,        0.5000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6491291128u,        0.6861018454u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6491291128u,        0.6861018454u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7681572011u,        0.6261611916u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7681572011u,        0.6261611916u) withcolor MidnightBlue;
draw (        0.7252162935u,        0.5000000001u)--(        0.7681572011u,        0.3738388086u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6491291128u,        0.3138981547u)--(        0.7681572011u,        0.3738388086u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7252162935u,        0.5000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7252162935u,        0.5000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6491291128u,        0.3138981547u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6491291128u,        0.3138981547u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7681572011u,        0.3738388086u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7681572011u,        0.3738388086u) withcolor MidnightBlue;
draw (        0.7252162935u,        0.5000000001u)--(        0.6061882052u,        0.4400593462u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6491291128u,        0.3138981547u)--(        0.6061882052u,        0.4400593462u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7252162935u,        0.5000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7252162935u,        0.5000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6491291128u,        0.3138981547u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6491291128u,        0.3138981547u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6061882052u,        0.4400593462u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6061882052u,        0.4400593462u) withcolor MidnightBlue;
draw (        0.7252162935u,        0.5000000001u)--(        0.6061837307u,        0.4400682322u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4871511679u,        0.4999999999u)--(        0.6061837307u,        0.4400682322u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7252162935u,        0.5000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7252162935u,        0.5000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4871511679u,        0.4999999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4871511679u,        0.4999999999u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6061837307u,        0.4400682322u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6061837307u,        0.4400682322u) withcolor MidnightBlue;
draw (        0.7252162935u,        0.5000000001u)--(        0.6061837306u,        0.5599317678u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4871511679u,        0.4999999999u)--(        0.6061837306u,        0.5599317678u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7252162935u,        0.5000000001u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7252162935u,        0.5000000001u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4871511679u,        0.4999999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4871511679u,        0.4999999999u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6061837306u,        0.5599317678u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6061837306u,        0.5599317678u) withcolor MidnightBlue;
draw (        0.6491291128u,        0.6861018454u)--(        0.6061774046u,        0.5599443306u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4871511679u,        0.4999999999u)--(        0.6061774046u,        0.5599443306u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6491291128u,        0.6861018454u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6491291128u,        0.6861018454u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4871511679u,        0.4999999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4871511679u,        0.4999999999u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6061774046u,        0.5599443306u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6061774046u,        0.5599443306u) withcolor MidnightBlue;
draw (        0.6491291128u,        0.6861018454u)--(        0.5301028762u,        0.6261575148u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4871511679u,        0.4999999999u)--(        0.5301028762u,        0.6261575148u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6491291128u,        0.6861018454u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6491291128u,        0.6861018454u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4871511679u,        0.4999999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4871511679u,        0.4999999999u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5301028762u,        0.6261575148u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5301028762u,        0.6261575148u) withcolor MidnightBlue;
draw (        0.6491291128u,        0.6861018454u)--(        0.5301000724u,        0.6261630822u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4752048465u,        0.7476006462u)--(        0.5301000724u,        0.6261630822u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6491291128u,        0.6861018454u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6491291128u,        0.6861018454u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4752048465u,        0.7476006462u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4752048465u,        0.7476006462u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5301000724u,        0.6261630822u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5301000724u,        0.6261630822u) withcolor MidnightBlue;
draw (        0.6491291128u,        0.6861018454u)--(        0.5942338869u,        0.8075394095u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4752048465u,        0.7476006462u)--(        0.5942338869u,        0.8075394095u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6491291128u,        0.6861018454u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6491291128u,        0.6861018454u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4752048465u,        0.7476006462u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4752048465u,        0.7476006462u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5942338869u,        0.8075394095u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5942338869u,        0.8075394095u) withcolor Green;
draw (        0.7222223226u,        0.1211369606u)--(        0.6666666667u,        0.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7222223226u,        0.1211369606u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7222223226u,        0.1211369606u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6666666667u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6666666667u,        0.0000000000u) withcolor Yellow;
draw (        0.7222223226u,        0.1211369606u)--(        0.7777777778u,        0.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7222223226u,        0.1211369606u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7222223226u,        0.1211369606u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7777777778u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7777777778u,        0.0000000000u) withcolor Yellow;
draw (        0.7222223226u,        0.1211369606u)--(        0.7646609373u,        0.2474680011u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6491291128u,        0.3138981547u)--(        0.7646609373u,        0.2474680011u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7222223226u,        0.1211369606u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7222223226u,        0.1211369606u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6491291128u,        0.3138981547u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6491291128u,        0.3138981547u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7646609373u,        0.2474680011u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7646609373u,        0.2474680011u) withcolor MidnightBlue;
draw (        0.6491291128u,        0.3138981547u)--(        0.5942338869u,        0.1924605907u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4752048465u,        0.2523993539u)--(        0.5942338869u,        0.1924605907u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6491291128u,        0.3138981547u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6491291128u,        0.3138981547u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4752048465u,        0.2523993539u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4752048465u,        0.2523993539u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5942338869u,        0.1924605907u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5942338869u,        0.1924605907u) withcolor Green;
draw (        0.6491291128u,        0.3138981547u)--(        0.5301000724u,        0.3738369179u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4752048465u,        0.2523993539u)--(        0.5301000724u,        0.3738369179u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6491291128u,        0.3138981547u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6491291128u,        0.3138981547u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4752048465u,        0.2523993539u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4752048465u,        0.2523993539u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5301000724u,        0.3738369179u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5301000724u,        0.3738369179u) withcolor MidnightBlue;
draw (        0.6491291128u,        0.3138981547u)--(        0.5301028760u,        0.3738424850u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4871511679u,        0.4999999999u)--(        0.5301028760u,        0.3738424850u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6491291128u,        0.3138981547u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6491291128u,        0.3138981547u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4871511679u,        0.4999999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4871511679u,        0.4999999999u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5301028760u,        0.3738424850u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5301028760u,        0.3738424850u) withcolor MidnightBlue;
draw (        0.6491291128u,        0.3138981547u)--(        0.6061774047u,        0.4400556696u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4871511679u,        0.4999999999u)--(        0.6061774047u,        0.4400556696u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6491291128u,        0.3138981547u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6491291128u,        0.3138981547u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4871511679u,        0.4999999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4871511679u,        0.4999999999u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6061774047u,        0.4400556696u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6061774047u,        0.4400556696u) withcolor MidnightBlue;
draw (        0.4752048465u,        0.2523993539u)--(        0.5300942515u,        0.3738395491u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4871511679u,        0.4999999999u)--(        0.5300942515u,        0.3738395491u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4752048465u,        0.2523993539u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4752048465u,        0.2523993539u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4871511679u,        0.4999999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4871511679u,        0.4999999999u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5300942515u,        0.3738395491u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5300942515u,        0.3738395491u) withcolor MidnightBlue;
draw (        0.4752048465u,        0.2523993539u)--(        0.3420452551u,        0.2470057449u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3064643645u,        0.3754369250u)--(        0.3420452551u,        0.2470057449u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4752048465u,        0.2523993539u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4752048465u,        0.2523993539u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3064643645u,        0.3754369250u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3064643645u,        0.3754369250u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3420452551u,        0.2470057449u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3420452551u,        0.2470057449u) withcolor MidnightBlue;
draw (        0.4752048465u,        0.2523993539u)--(        0.3420457876u,        0.2469926131u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2878089771u,        0.1252595602u)--(        0.3420457876u,        0.2469926131u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4752048465u,        0.2523993539u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4752048465u,        0.2523993539u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2878089771u,        0.1252595602u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2878089771u,        0.1252595602u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3420457876u,        0.2469926131u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3420457876u,        0.2469926131u) withcolor MidnightBlue;
draw (        0.4871511679u,        0.4999999999u)--(        0.5300942512u,        0.6261604509u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4752048465u,        0.7476006462u)--(        0.5300942512u,        0.6261604509u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4871511679u,        0.4999999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4871511679u,        0.4999999999u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4752048465u,        0.7476006462u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4752048465u,        0.7476006462u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5300942512u,        0.6261604509u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5300942512u,        0.6261604509u) withcolor MidnightBlue;
draw (        0.4871511679u,        0.4999999999u)--(        0.3538823885u,        0.4999845014u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3064643645u,        0.3754369250u)--(        0.3538823885u,        0.4999845014u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4871511679u,        0.4999999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4871511679u,        0.4999999999u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3064643645u,        0.3754369250u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3064643645u,        0.3754369250u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3538823885u,        0.4999845014u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3538823885u,        0.4999845014u) withcolor MidnightBlue;
draw (        0.4871511679u,        0.4999999999u)--(        0.3538823885u,        0.5000154985u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3064643646u,        0.6245630749u)--(        0.3538823885u,        0.5000154985u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4871511679u,        0.4999999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4871511679u,        0.4999999999u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3064643646u,        0.6245630749u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3064643646u,        0.6245630749u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3538823885u,        0.5000154985u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3538823885u,        0.5000154985u) withcolor MidnightBlue;
draw (        0.4752048465u,        0.7476006462u)--(        0.3420452551u,        0.7529942550u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3064643646u,        0.6245630749u)--(        0.3420452551u,        0.7529942550u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4752048465u,        0.7476006462u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4752048465u,        0.7476006462u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3064643646u,        0.6245630749u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3064643646u,        0.6245630749u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3420452551u,        0.7529942550u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3420452551u,        0.7529942550u) withcolor MidnightBlue;
draw (        0.4752048465u,        0.7476006462u)--(        0.3420457877u,        0.7530073869u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2878089772u,        0.8747404399u)--(        0.3420457877u,        0.7530073869u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4752048465u,        0.7476006462u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4752048465u,        0.7476006462u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2878089772u,        0.8747404399u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2878089772u,        0.8747404399u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3420457877u,        0.7530073869u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3420457877u,        0.7530073869u) withcolor MidnightBlue;
draw (        0.3064643645u,        0.3754369250u)--(        0.1924500897u,        0.4444444444u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3064643645u,        0.3754369250u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3064643645u,        0.3754369250u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1924500897u,        0.4444444444u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1924500897u,        0.4444444444u) withcolor Yellow;
draw (        0.3064643645u,        0.3754369250u)--(        0.3538416602u,        0.4999999999u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3064643646u,        0.6245630749u)--(        0.3538416602u,        0.4999999999u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3064643645u,        0.3754369250u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3064643645u,        0.3754369250u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3064643646u,        0.6245630749u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3064643646u,        0.6245630749u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3538416602u,        0.4999999999u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3538416602u,        0.4999999999u) withcolor MidnightBlue;
draw (        0.3064643645u,        0.3754369250u)--(        0.3420274059u,        0.2470008012u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2878089771u,        0.1252595602u)--(        0.3420274059u,        0.2470008012u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3064643645u,        0.3754369250u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3064643645u,        0.3754369250u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2878089771u,        0.1252595602u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2878089771u,        0.1252595602u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3420274059u,        0.2470008012u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3420274059u,        0.2470008012u) withcolor MidnightBlue;
draw (        0.3064643645u,        0.3754369250u)--(        0.2522459358u,        0.2536956839u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2878089771u,        0.1252595602u)--(        0.2522459358u,        0.2536956839u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3064643645u,        0.3754369250u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3064643645u,        0.3754369250u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2878089771u,        0.1252595602u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2878089771u,        0.1252595602u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2522459358u,        0.2536956839u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2522459358u,        0.2536956839u) withcolor MidnightBlue;
draw (        0.3064643645u,        0.3754369250u)--(        0.2522181707u,        0.2537080531u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1211404611u,        0.2777746659u)--(        0.2522181707u,        0.2537080531u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3064643645u,        0.3754369250u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3064643645u,        0.3754369250u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1211404611u,        0.2777746659u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1211404611u,        0.2777746659u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2522181707u,        0.2537080531u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2522181707u,        0.2537080531u) withcolor MidnightBlue;
draw (        0.3064643646u,        0.6245630749u)--(        0.1924500897u,        0.5555555556u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3064643646u,        0.6245630749u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3064643646u,        0.6245630749u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1924500897u,        0.5555555556u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1924500897u,        0.5555555556u) withcolor Yellow;
draw (        0.3064643646u,        0.6245630749u)--(        0.2522459360u,        0.7463043161u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2878089772u,        0.8747404399u)--(        0.2522459360u,        0.7463043161u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3064643646u,        0.6245630749u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3064643646u,        0.6245630749u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2878089772u,        0.8747404399u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2878089772u,        0.8747404399u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2522459360u,        0.7463043161u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2522459360u,        0.7463043161u) withcolor MidnightBlue;
draw (        0.3064643646u,        0.6245630749u)--(        0.3420274057u,        0.7529991987u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2878089772u,        0.8747404399u)--(        0.3420274057u,        0.7529991987u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3064643646u,        0.6245630749u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3064643646u,        0.6245630749u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2878089772u,        0.8747404399u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2878089772u,        0.8747404399u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3420274057u,        0.7529991987u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3420274057u,        0.7529991987u) withcolor MidnightBlue;
draw (        0.3064643646u,        0.6245630749u)--(        0.2522181706u,        0.7462919467u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1211404610u,        0.7222253340u)--(        0.2522181706u,        0.7462919467u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3064643646u,        0.6245630749u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3064643646u,        0.6245630749u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1211404610u,        0.7222253340u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1211404610u,        0.7222253340u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2522181706u,        0.7462919467u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2522181706u,        0.7462919467u) withcolor MidnightBlue;
draw (        0.2878089772u,        0.8747404399u)--(        0.3333333333u,        1.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2878089772u,        0.8747404399u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2878089772u,        0.8747404399u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3333333333u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3333333333u,        1.0000000000u) withcolor Yellow;
draw (        0.2878089772u,        0.8747404399u)--(        0.2422807985u,        0.9999911993u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1211383817u,        0.9444476789u)--(        0.2422807985u,        0.9999911993u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2878089772u,        0.8747404399u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2878089772u,        0.8747404399u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1211383817u,        0.9444476789u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1211383817u,        0.9444476789u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2422807985u,        0.9999911993u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2422807985u,        0.9999911993u) withcolor Green;
draw (        0.2878089772u,        0.8747404399u)--(        0.2522142880u,        0.7463130834u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1211404610u,        0.7222253340u)--(        0.2522142880u,        0.7463130834u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2878089772u,        0.8747404399u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2878089772u,        0.8747404399u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1211404610u,        0.7222253340u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1211404610u,        0.7222253340u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2522142880u,        0.7463130834u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2522142880u,        0.7463130834u) withcolor MidnightBlue;
draw (        0.1211383817u,        0.9444476789u)--(        0.0000000000u,        1.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1211383817u,        0.9444476789u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1211383817u,        0.9444476789u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        1.0000000000u) withcolor Yellow;
draw (        0.1211383817u,        0.9444476789u)--(        0.0000000000u,        0.8888888889u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1211383817u,        0.9444476789u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1211383817u,        0.9444476789u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.8888888889u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.8888888889u) withcolor Yellow;
draw (        0.2878089771u,        0.1252595602u)--(        0.3333333333u,        0.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2878089771u,        0.1252595602u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2878089771u,        0.1252595602u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3333333333u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3333333333u,        0.0000000000u) withcolor Yellow;
draw (        0.2878089771u,        0.1252595602u)--(        0.2422807983u,        0.0000088008u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1211383816u,        0.0555523212u)--(        0.2422807983u,        0.0000088008u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2878089771u,        0.1252595602u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2878089771u,        0.1252595602u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1211383816u,        0.0555523212u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1211383816u,        0.0555523212u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2422807983u,        0.0000088008u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2422807983u,        0.0000088008u) withcolor Green;
draw (        0.2878089771u,        0.1252595602u)--(        0.2522142881u,        0.2536869167u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1211404611u,        0.2777746659u)--(        0.2522142881u,        0.2536869167u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2878089771u,        0.1252595602u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2878089771u,        0.1252595602u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1211404611u,        0.2777746659u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1211404611u,        0.2777746659u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2522142881u,        0.2536869167u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2522142881u,        0.2536869167u) withcolor MidnightBlue;
draw (        0.1211404610u,        0.7222253340u)--(        0.0000000000u,        0.7777777778u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1211404610u,        0.7222253340u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1211404610u,        0.7222253340u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.7777777778u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.7777777778u) withcolor Yellow;
draw (        0.1211404610u,        0.7222253340u)--(        0.0000000000u,        0.6666666667u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1211404610u,        0.7222253340u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1211404610u,        0.7222253340u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.6666666667u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.6666666667u) withcolor Yellow;
draw (        0.1211383816u,        0.0555523212u)--(        0.0000000000u,        0.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1211383816u,        0.0555523212u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1211383816u,        0.0555523212u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.1211383816u,        0.0555523212u)--(        0.0000000000u,        0.1111111111u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1211383816u,        0.0555523212u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1211383816u,        0.0555523212u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.1111111111u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.1111111111u) withcolor Yellow;
draw (        0.1211404611u,        0.2777746659u)--(        0.0000000000u,        0.3333333333u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1211404611u,        0.2777746659u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1211404611u,        0.2777746659u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.3333333333u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.3333333333u) withcolor Yellow;
draw (        0.1211404611u,        0.2777746659u)--(        0.0000000000u,        0.2222222222u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1211404611u,        0.2777746659u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1211404611u,        0.2777746659u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.2222222222u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.2222222222u) withcolor Yellow;
% r =         0.1332687803
endfig;
end;
