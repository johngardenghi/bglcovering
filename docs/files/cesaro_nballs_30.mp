prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(30);
u := 20cm;
path pols[];
pair p[];
R := .1052319205u;
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

x1 = .9149954260u;
y1 = .7157450788u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .9118422444u;
y2 = .9445752229u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .8226753758u;
y3 = .8015505302u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .7157439321u;
y4 = .9149951365u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .7356969796u;
y5 = .6238861432u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .6428870097u;
y6 = .7358405984u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .9149789919u;
y7 = .2842365060u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .7509667924u;
y8 = .3613305434u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .6700907592u;
y9 = .4443257749u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = .5648585032u;
y10 = .5492450609u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;

x11 = .4999999999u;
y11 = .7254820138u;
draw fullcircle scaled D shifted z11 withcolor RoyalBlue withpen pencircle scaled 1.0;

x12 = .2842560679u;
y12 = .9149951366u;
draw fullcircle scaled D shifted z12 withcolor RoyalBlue withpen pencircle scaled 1.0;

x13 = .9256280076u;
y13 = .0745155666u;
draw fullcircle scaled D shifted z13 withcolor RoyalBlue withpen pencircle scaled 1.0;

x14 = .8085724098u;
y14 = .1902558381u;
draw fullcircle scaled D shifted z14 withcolor RoyalBlue withpen pencircle scaled 1.0;

x15 = .6232558851u;
y15 = .2577806083u;
draw fullcircle scaled D shifted z15 withcolor RoyalBlue withpen pencircle scaled 1.0;

x16 = .5000000001u;
y16 = .3611431511u;
draw fullcircle scaled D shifted z16 withcolor RoyalBlue withpen pencircle scaled 1.0;

x17 = .4351414967u;
y17 = .5492450610u;
draw fullcircle scaled D shifted z17 withcolor RoyalBlue withpen pencircle scaled 1.0;

x18 = .3571129903u;
y18 = .7358405983u;
draw fullcircle scaled D shifted z18 withcolor RoyalBlue withpen pencircle scaled 1.0;

x19 = .1773246241u;
y19 = .8015505302u;
draw fullcircle scaled D shifted z19 withcolor RoyalBlue withpen pencircle scaled 1.0;

x20 = .0881577556u;
y20 = .9445752228u;
draw fullcircle scaled D shifted z20 withcolor RoyalBlue withpen pencircle scaled 1.0;

x21 = .7157632804u;
y21 = .0850205270u;
draw fullcircle scaled D shifted z21 withcolor RoyalBlue withpen pencircle scaled 1.0;

x22 = .3767441149u;
y22 = .2577806084u;
draw fullcircle scaled D shifted z22 withcolor RoyalBlue withpen pencircle scaled 1.0;

x23 = .3299092410u;
y23 = .4443257748u;
draw fullcircle scaled D shifted z23 withcolor RoyalBlue withpen pencircle scaled 1.0;

x24 = .2643030202u;
y24 = .6238861431u;
draw fullcircle scaled D shifted z24 withcolor RoyalBlue withpen pencircle scaled 1.0;

x25 = .0850045740u;
y25 = .7157450790u;
draw fullcircle scaled D shifted z25 withcolor RoyalBlue withpen pencircle scaled 1.0;

x26 = .1914275902u;
y26 = .1902558381u;
draw fullcircle scaled D shifted z26 withcolor RoyalBlue withpen pencircle scaled 1.0;

x27 = .2490332076u;
y27 = .3613305435u;
draw fullcircle scaled D shifted z27 withcolor RoyalBlue withpen pencircle scaled 1.0;

x28 = .2842367195u;
y28 = .0850205270u;
draw fullcircle scaled D shifted z28 withcolor RoyalBlue withpen pencircle scaled 1.0;

x29 = .0850210082u;
y29 = .2842365060u;
draw fullcircle scaled D shifted z29 withcolor RoyalBlue withpen pencircle scaled 1.0;

x30 = .0743719922u;
y30 = .0745155666u;
draw fullcircle scaled D shifted z30 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.9149954260u,        0.7157450788u)--(        0.9037749551u,        0.6111111111u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9149954260u,        0.7157450788u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9149954260u,        0.7157450788u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9037749551u,        0.6111111111u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9037749551u,        0.6111111111u) withcolor Yellow;
draw (        0.9149954260u,        0.7157450788u)--(        1.0000000000u,        0.7777777778u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9149954260u,        0.7157450788u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9149954260u,        0.7157450788u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        0.7777777778u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        0.7777777778u) withcolor Yellow;
draw (        0.9149954260u,        0.7157450788u)--(        0.8114614152u,        0.6969178193u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8226753758u,        0.8015505302u)--(        0.8114614152u,        0.6969178193u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9149954260u,        0.7157450788u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9149954260u,        0.7157450788u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8226753758u,        0.8015505302u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8226753758u,        0.8015505302u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8114614152u,        0.6969178193u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8114614152u,        0.6969178193u) withcolor MidnightBlue;
draw (        0.9149954260u,        0.7157450788u)--(        0.9262093866u,        0.8203777897u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8226753758u,        0.8015505302u)--(        0.9262093866u,        0.8203777897u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9149954260u,        0.7157450788u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9149954260u,        0.7157450788u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8226753758u,        0.8015505302u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8226753758u,        0.8015505302u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9262093866u,        0.8203777897u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9262093866u,        0.8203777897u) withcolor Green;
draw (        0.9149954260u,        0.7157450788u)--(        0.8114615967u,        0.6969168210u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7356969796u,        0.6238861432u)--(        0.8114615967u,        0.6969168210u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9149954260u,        0.7157450788u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9149954260u,        0.7157450788u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7356969796u,        0.6238861432u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7356969796u,        0.6238861432u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8114615967u,        0.6969168210u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8114615967u,        0.6969168210u) withcolor MidnightBlue;
draw (        0.8226753758u,        0.8015505302u)--(        0.7177659415u,        0.8097826440u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7157439321u,        0.9149951365u)--(        0.7177659415u,        0.8097826440u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8226753758u,        0.8015505302u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8226753758u,        0.8015505302u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7157439321u,        0.9149951365u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7157439321u,        0.9149951365u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7177659415u,        0.8097826440u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7177659415u,        0.8097826440u) withcolor MidnightBlue;
draw (        0.8226753758u,        0.8015505302u)--(        0.8114605445u,        0.6969179126u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7356969796u,        0.6238861432u)--(        0.8114605445u,        0.6969179126u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8226753758u,        0.8015505302u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8226753758u,        0.8015505302u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7356969796u,        0.6238861432u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7356969796u,        0.6238861432u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8114605445u,        0.6969179126u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8114605445u,        0.6969179126u) withcolor MidnightBlue;
draw (        0.8226753758u,        0.8015505302u)--(        0.7177656843u,        0.8097793648u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6428870097u,        0.7358405984u)--(        0.7177656843u,        0.8097793648u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8226753758u,        0.8015505302u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8226753758u,        0.8015505302u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6428870097u,        0.7358405984u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6428870097u,        0.7358405984u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7177656843u,        0.8097793648u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7177656843u,        0.8097793648u) withcolor MidnightBlue;
draw (        0.7157439321u,        0.9149951365u)--(        0.7777777778u,        1.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7157439321u,        0.9149951365u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7157439321u,        0.9149951365u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7777777778u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7777777778u,        1.0000000000u) withcolor Yellow;
draw (        0.7157439321u,        0.9149951365u)--(        0.6111111111u,        0.9037749551u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7157439321u,        0.9149951365u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7157439321u,        0.9149951365u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6111111111u,        0.9037749551u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6111111111u,        0.9037749551u) withcolor Yellow;
draw (        0.7157439321u,        0.9149951365u)--(        0.7177625112u,        0.8097825781u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6428870097u,        0.7358405984u)--(        0.7177625112u,        0.8097825781u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7157439321u,        0.9149951365u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7157439321u,        0.9149951365u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6428870097u,        0.7358405984u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6428870097u,        0.7358405984u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7177625112u,        0.8097825781u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7177625112u,        0.8097825781u) withcolor MidnightBlue;
draw (        0.7356969796u,        0.6238861432u)--(        0.6307275045u,        0.6313135515u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6428870097u,        0.7358405984u)--(        0.6307275045u,        0.6313135515u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7356969796u,        0.6238861432u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7356969796u,        0.6238861432u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6428870097u,        0.7358405984u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6428870097u,        0.7358405984u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6307275045u,        0.6313135515u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6307275045u,        0.6313135515u) withcolor MidnightBlue;
draw (        0.7356969796u,        0.6238861432u)--(        0.7442348456u,        0.5190011493u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6700907592u,        0.4443257749u)--(        0.7442348456u,        0.5190011493u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7356969796u,        0.6238861432u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7356969796u,        0.6238861432u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6700907592u,        0.4443257749u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6700907592u,        0.4443257749u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7442348456u,        0.5190011493u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7442348456u,        0.5190011493u) withcolor Green;
draw (        0.7356969796u,        0.6238861432u)--(        0.6307274208u,        0.6313123674u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5648585032u,        0.5492450609u)--(        0.6307274208u,        0.6313123674u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7356969796u,        0.6238861432u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7356969796u,        0.6238861432u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5648585032u,        0.5492450609u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5648585032u,        0.5492450609u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6307274208u,        0.6313123674u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6307274208u,        0.6313123674u) withcolor MidnightBlue;
draw (        0.6428870097u,        0.7358405984u)--(        0.6307256812u,        0.6313137636u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5648585032u,        0.5492450609u)--(        0.6307256812u,        0.6313137636u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6428870097u,        0.7358405984u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6428870097u,        0.7358405984u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5648585032u,        0.5492450609u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5648585032u,        0.5492450609u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6307256812u,        0.6313137636u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6307256812u,        0.6313137636u) withcolor MidnightBlue;
draw (        0.6428870097u,        0.7358405984u)--(        0.5770198316u,        0.6537718957u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5648585032u,        0.5492450609u)--(        0.5770198316u,        0.6537718957u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6428870097u,        0.7358405984u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6428870097u,        0.7358405984u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5648585032u,        0.5492450609u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5648585032u,        0.5492450609u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5770198316u,        0.6537718957u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5770198316u,        0.6537718957u) withcolor MidnightBlue;
draw (        0.6428870097u,        0.7358405984u)--(        0.5770174597u,        0.6537737995u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4999999999u,        0.7254820138u)--(        0.5770174597u,        0.6537737995u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6428870097u,        0.7358405984u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6428870097u,        0.7358405984u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4999999999u,        0.7254820138u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999999u,        0.7254820138u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5770174597u,        0.6537737995u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5770174597u,        0.6537737995u) withcolor MidnightBlue;
draw (        0.6428870097u,        0.7358405984u)--(        0.5658695499u,        0.8075488126u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4999999999u,        0.7254820138u)--(        0.5658695499u,        0.8075488126u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6428870097u,        0.7358405984u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6428870097u,        0.7358405984u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4999999999u,        0.7254820138u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999999u,        0.7254820138u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5658695499u,        0.8075488126u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5658695499u,        0.8075488126u) withcolor Green;
draw (        0.9149789919u,        0.2842365060u)--(        1.0000000000u,        0.2222222222u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9149789919u,        0.2842365060u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9149789919u,        0.2842365060u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        0.2222222222u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        0.2222222222u) withcolor Yellow;
draw (        0.9149789919u,        0.2842365060u)--(        0.9037749551u,        0.3888888889u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9149789919u,        0.2842365060u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9149789919u,        0.2842365060u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9037749551u,        0.3888888889u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9037749551u,        0.3888888889u) withcolor Yellow;
draw (        0.9149789919u,        0.2842365060u)--(        0.9132635546u,        0.1790185685u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9256280076u,        0.0745155666u)--(        0.9132635546u,        0.1790185685u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9149789919u,        0.2842365060u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9149789919u,        0.2842365060u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9256280076u,        0.0745155666u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9256280076u,        0.0745155666u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9132635546u,        0.1790185685u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9132635546u,        0.1790185685u) withcolor MidnightBlue;
draw (        0.9149789919u,        0.2842365060u)--(        0.9132027283u,        0.1790195778u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8085724098u,        0.1902558381u)--(        0.9132027283u,        0.1790195778u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9149789919u,        0.2842365060u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9149789919u,        0.2842365060u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8085724098u,        0.1902558381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8085724098u,        0.1902558381u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9132027283u,        0.1790195778u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9132027283u,        0.1790195778u) withcolor MidnightBlue;
draw (        0.7509667924u,        0.3613305434u)--(        0.7734414878u,        0.4641344595u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6700907592u,        0.4443257749u)--(        0.7734414878u,        0.4641344595u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7509667924u,        0.3613305434u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7509667924u,        0.3613305434u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6700907592u,        0.4443257749u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6700907592u,        0.4443257749u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7734414878u,        0.4641344595u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7734414878u,        0.4641344595u) withcolor Green;
draw (        0.7509667924u,        0.3613305434u)--(        0.7284920973u,        0.2585266272u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8085724098u,        0.1902558381u)--(        0.7284920973u,        0.2585266272u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7509667924u,        0.3613305434u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7509667924u,        0.3613305434u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8085724098u,        0.1902558381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8085724098u,        0.1902558381u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7284920973u,        0.2585266272u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7284920973u,        0.2585266272u) withcolor MidnightBlue;
draw (        0.7509667924u,        0.3613305434u)--(        0.7284851505u,        0.2585281462u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6232558851u,        0.2577806083u)--(        0.7284851505u,        0.2585281462u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7509667924u,        0.3613305434u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7509667924u,        0.3613305434u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6232558851u,        0.2577806083u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6232558851u,        0.2577806083u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7284851505u,        0.2585281462u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7284851505u,        0.2585281462u) withcolor MidnightBlue;
draw (        0.6700907592u,        0.4443257749u)--(        0.5648593030u,        0.4440131404u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5648585032u,        0.5492450609u)--(        0.5648593030u,        0.4440131404u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6700907592u,        0.4443257749u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6700907592u,        0.4443257749u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5648585032u,        0.5492450609u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5648585032u,        0.5492450609u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5648593030u,        0.4440131404u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5648593030u,        0.4440131404u) withcolor MidnightBlue;
draw (        0.6700907592u,        0.4443257749u)--(        0.6052323375u,        0.3614575528u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6232558851u,        0.2577806083u)--(        0.6052323375u,        0.3614575528u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6700907592u,        0.4443257749u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6700907592u,        0.4443257749u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6232558851u,        0.2577806083u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6232558851u,        0.2577806083u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6052323375u,        0.3614575528u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6052323375u,        0.3614575528u) withcolor MidnightBlue;
draw (        0.6700907592u,        0.4443257749u)--(        0.6052314489u,        0.3614582483u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5000000001u,        0.3611431511u)--(        0.6052314489u,        0.3614582483u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6700907592u,        0.4443257749u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6700907592u,        0.4443257749u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6052314489u,        0.3614582483u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6052314489u,        0.3614582483u) withcolor MidnightBlue;
draw (        0.6700907592u,        0.4443257749u)--(        0.5648593104u,        0.4440106777u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5000000001u,        0.3611431511u)--(        0.5648593104u,        0.4440106777u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6700907592u,        0.4443257749u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6700907592u,        0.4443257749u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5648593104u,        0.4440106777u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5648593104u,        0.4440106777u) withcolor MidnightBlue;
draw (        0.5648585032u,        0.5492450609u)--(        0.5770160921u,        0.6537723307u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4999999999u,        0.7254820138u)--(        0.5770160921u,        0.6537723307u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5648585032u,        0.5492450609u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5648585032u,        0.5492450609u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4999999999u,        0.7254820138u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999999u,        0.7254820138u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5770160921u,        0.6537723307u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5770160921u,        0.6537723307u) withcolor MidnightBlue;
draw (        0.5648585032u,        0.5492450609u)--(        0.5648561638u,        0.4440131404u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5000000001u,        0.3611431511u)--(        0.5648561638u,        0.4440131404u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5648585032u,        0.5492450609u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5648585032u,        0.5492450609u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5648561638u,        0.4440131404u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5648561638u,        0.4440131404u) withcolor MidnightBlue;
draw (        0.5648585032u,        0.5492450609u)--(        0.5000023395u,        0.4663750716u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5000000001u,        0.3611431511u)--(        0.5000023395u,        0.4663750716u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5648585032u,        0.5492450609u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5648585032u,        0.5492450609u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000023395u,        0.4663750716u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000023395u,        0.4663750716u) withcolor MidnightBlue;
draw (        0.5648585032u,        0.5492450609u)--(        0.4999999999u,        0.4663769027u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4351414967u,        0.5492450610u)--(        0.4999999999u,        0.4663769027u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5648585032u,        0.5492450609u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5648585032u,        0.5492450609u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4351414967u,        0.5492450610u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4351414967u,        0.5492450610u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4999999999u,        0.4663769027u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999999u,        0.4663769027u) withcolor MidnightBlue;
draw (        0.4999999999u,        0.7254820138u)--(        0.4229839077u,        0.6537723308u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4351414967u,        0.5492450610u)--(        0.4229839077u,        0.6537723308u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4999999999u,        0.7254820138u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999999u,        0.7254820138u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4351414967u,        0.5492450610u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4351414967u,        0.5492450610u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4229839077u,        0.6537723308u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4229839077u,        0.6537723308u) withcolor MidnightBlue;
draw (        0.4999999999u,        0.7254820138u)--(        0.4229825403u,        0.6537737994u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3571129903u,        0.7358405983u)--(        0.4229825403u,        0.6537737994u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4999999999u,        0.7254820138u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999999u,        0.7254820138u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3571129903u,        0.7358405983u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3571129903u,        0.7358405983u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4229825403u,        0.6537737994u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4229825403u,        0.6537737994u) withcolor MidnightBlue;
draw (        0.4999999999u,        0.7254820138u)--(        0.4341304499u,        0.8075488126u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3571129903u,        0.7358405983u)--(        0.4341304499u,        0.8075488126u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4999999999u,        0.7254820138u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999999u,        0.7254820138u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3571129903u,        0.7358405983u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3571129903u,        0.7358405983u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4341304499u,        0.8075488126u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4341304499u,        0.8075488126u) withcolor Green;
draw (        0.2842560679u,        0.9149951366u)--(        0.3888888889u,        0.9037749551u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2842560679u,        0.9149951366u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2842560679u,        0.9149951366u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3888888889u,        0.9037749551u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3888888889u,        0.9037749551u) withcolor Yellow;
draw (        0.2842560679u,        0.9149951366u)--(        0.2222222222u,        1.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2842560679u,        0.9149951366u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2842560679u,        0.9149951366u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2222222222u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2222222222u,        1.0000000000u) withcolor Yellow;
draw (        0.2842560679u,        0.9149951366u)--(        0.2822374889u,        0.8097825782u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3571129903u,        0.7358405983u)--(        0.2822374889u,        0.8097825782u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2842560679u,        0.9149951366u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2842560679u,        0.9149951366u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3571129903u,        0.7358405983u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3571129903u,        0.7358405983u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2822374889u,        0.8097825782u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2822374889u,        0.8097825782u) withcolor MidnightBlue;
draw (        0.2842560679u,        0.9149951366u)--(        0.2822340584u,        0.8097826441u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1773246241u,        0.8015505302u)--(        0.2822340584u,        0.8097826441u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2842560679u,        0.9149951366u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2842560679u,        0.9149951366u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1773246241u,        0.8015505302u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1773246241u,        0.8015505302u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2822340584u,        0.8097826441u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2822340584u,        0.8097826441u) withcolor MidnightBlue;
draw (        0.9256280076u,        0.0745155666u)--(        1.0000000000u,        0.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9256280076u,        0.0745155666u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9256280076u,        0.0745155666u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.9256280076u,        0.0745155666u)--(        0.8209985841u,        0.0857601574u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8085724098u,        0.1902558381u)--(        0.8209985841u,        0.0857601574u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9256280076u,        0.0745155666u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9256280076u,        0.0745155666u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8085724098u,        0.1902558381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8085724098u,        0.1902558381u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8209985841u,        0.0857601574u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8209985841u,        0.0857601574u) withcolor MidnightBlue;
draw (        0.9256280076u,        0.0745155666u)--(        0.9132018334u,        0.1790112474u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8085724098u,        0.1902558381u)--(        0.9132018334u,        0.1790112474u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9256280076u,        0.0745155666u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9256280076u,        0.0745155666u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8085724098u,        0.1902558381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8085724098u,        0.1902558381u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9132018334u,        0.1790112474u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9132018334u,        0.1790112474u) withcolor MidnightBlue;
draw (        0.9256280076u,        0.0745155666u)--(        0.8203983306u,        0.0738284156u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7157632804u,        0.0850205270u)--(        0.8203983306u,        0.0738284156u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9256280076u,        0.0745155666u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9256280076u,        0.0745155666u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7157632804u,        0.0850205270u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7157632804u,        0.0850205270u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8203983306u,        0.0738284156u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8203983306u,        0.0738284156u) withcolor Green;
draw (        0.9256280076u,        0.0745155666u)--(        0.8209929574u,        0.0857076780u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7157632804u,        0.0850205270u)--(        0.8209929574u,        0.0857076780u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9256280076u,        0.0745155666u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9256280076u,        0.0745155666u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7157632804u,        0.0850205270u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7157632804u,        0.0850205270u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8209929574u,        0.0857076780u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8209929574u,        0.0857076780u) withcolor MidnightBlue;
draw (        0.8085724098u,        0.1902558381u)--(        0.7033430768u,        0.1895178891u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6232558851u,        0.2577806083u)--(        0.7033430768u,        0.1895178891u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8085724098u,        0.1902558381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8085724098u,        0.1902558381u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6232558851u,        0.2577806083u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6232558851u,        0.2577806083u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7033430768u,        0.1895178891u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7033430768u,        0.1895178891u) withcolor MidnightBlue;
draw (        0.8085724098u,        0.1902558381u)--(        0.7284852182u,        0.2585185573u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6232558851u,        0.2577806083u)--(        0.7284852182u,        0.2585185573u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8085724098u,        0.1902558381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8085724098u,        0.1902558381u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6232558851u,        0.2577806083u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6232558851u,        0.2577806083u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7284852182u,        0.2585185573u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7284852182u,        0.2585185573u) withcolor MidnightBlue;
draw (        0.8085724098u,        0.1902558381u)--(        0.8209926066u,        0.0857594467u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7157632804u,        0.0850205270u)--(        0.8209926066u,        0.0857594467u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8085724098u,        0.1902558381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8085724098u,        0.1902558381u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7157632804u,        0.0850205270u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7157632804u,        0.0850205270u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8209926066u,        0.0857594467u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8209926066u,        0.0857594467u) withcolor MidnightBlue;
draw (        0.8085724098u,        0.1902558381u)--(        0.7033430836u,        0.1895169183u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7157632804u,        0.0850205270u)--(        0.7033430836u,        0.1895169183u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8085724098u,        0.1902558381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8085724098u,        0.1902558381u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7157632804u,        0.0850205270u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7157632804u,        0.0850205270u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7033430836u,        0.1895169183u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7033430836u,        0.1895169183u) withcolor MidnightBlue;
draw (        0.6232558851u,        0.2577806083u)--(        0.5180244338u,        0.2574663606u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5000000001u,        0.3611431511u)--(        0.5180244338u,        0.2574663606u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6232558851u,        0.2577806083u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6232558851u,        0.2577806083u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5180244338u,        0.2574663606u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5180244338u,        0.2574663606u) withcolor Green;
draw (        0.6232558851u,        0.2577806083u)--(        0.6052314514u,        0.3614573988u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5000000001u,        0.3611431511u)--(        0.6052314514u,        0.3614573988u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6232558851u,        0.2577806083u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6232558851u,        0.2577806083u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6052314514u,        0.3614573988u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6052314514u,        0.3614573988u) withcolor MidnightBlue;
draw (        0.6232558851u,        0.2577806083u)--(        0.7033421553u,        0.1895168080u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7157632804u,        0.0850205270u)--(        0.7033421553u,        0.1895168080u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6232558851u,        0.2577806083u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6232558851u,        0.2577806083u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7157632804u,        0.0850205270u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7157632804u,        0.0850205270u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7033421553u,        0.1895168080u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7033421553u,        0.1895168080u) withcolor MidnightBlue;
draw (        0.5000000001u,        0.3611431511u)--(        0.4351438365u,        0.4440131405u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4351414967u,        0.5492450610u)--(        0.4351438365u,        0.4440131405u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4351414967u,        0.5492450610u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4351414967u,        0.5492450610u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4351438365u,        0.4440131405u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4351438365u,        0.4440131405u) withcolor MidnightBlue;
draw (        0.5000000001u,        0.3611431511u)--(        0.4999976603u,        0.4663750716u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4351414967u,        0.5492450610u)--(        0.4999976603u,        0.4663750716u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4351414967u,        0.5492450610u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4351414967u,        0.5492450610u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4999976603u,        0.4663750716u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999976603u,        0.4663750716u) withcolor MidnightBlue;
draw (        0.5000000001u,        0.3611431511u)--(        0.4819755662u,        0.2574663607u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3767441149u,        0.2577806084u)--(        0.4819755662u,        0.2574663607u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3767441149u,        0.2577806084u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3767441149u,        0.2577806084u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4819755662u,        0.2574663607u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4819755662u,        0.2574663607u) withcolor Green;
draw (        0.5000000001u,        0.3611431511u)--(        0.3947685488u,        0.3614573988u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3767441149u,        0.2577806084u)--(        0.3947685488u,        0.3614573988u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3767441149u,        0.2577806084u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3767441149u,        0.2577806084u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3947685488u,        0.3614573988u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3947685488u,        0.3614573988u) withcolor MidnightBlue;
draw (        0.5000000001u,        0.3611431511u)--(        0.3947685513u,        0.3614582482u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3299092410u,        0.4443257748u)--(        0.3947685513u,        0.3614582482u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3299092410u,        0.4443257748u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3299092410u,        0.4443257748u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3947685513u,        0.3614582482u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3947685513u,        0.3614582482u) withcolor MidnightBlue;
draw (        0.5000000001u,        0.3611431511u)--(        0.4351406897u,        0.4440106777u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3299092410u,        0.4443257748u)--(        0.4351406897u,        0.4440106777u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.3611431511u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3299092410u,        0.4443257748u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3299092410u,        0.4443257748u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4351406897u,        0.4440106777u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4351406897u,        0.4440106777u) withcolor MidnightBlue;
draw (        0.4351414967u,        0.5492450610u)--(        0.3692743185u,        0.6313137635u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3571129903u,        0.7358405983u)--(        0.3692743185u,        0.6313137635u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4351414967u,        0.5492450610u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4351414967u,        0.5492450610u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3571129903u,        0.7358405983u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3571129903u,        0.7358405983u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3692743185u,        0.6313137635u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3692743185u,        0.6313137635u) withcolor MidnightBlue;
draw (        0.4351414967u,        0.5492450610u)--(        0.4229801685u,        0.6537718958u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3571129903u,        0.7358405983u)--(        0.4229801685u,        0.6537718958u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4351414967u,        0.5492450610u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4351414967u,        0.5492450610u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3571129903u,        0.7358405983u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3571129903u,        0.7358405983u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4229801685u,        0.6537718958u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4229801685u,        0.6537718958u) withcolor MidnightBlue;
draw (        0.4351414967u,        0.5492450610u)--(        0.4351406971u,        0.4440131405u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3299092410u,        0.4443257748u)--(        0.4351406971u,        0.4440131405u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4351414967u,        0.5492450610u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4351414967u,        0.5492450610u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3299092410u,        0.4443257748u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3299092410u,        0.4443257748u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4351406971u,        0.4440131405u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4351406971u,        0.4440131405u) withcolor MidnightBlue;
draw (        0.4351414967u,        0.5492450610u)--(        0.3692725791u,        0.6313123674u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2643030202u,        0.6238861431u)--(        0.3692725791u,        0.6313123674u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4351414967u,        0.5492450610u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4351414967u,        0.5492450610u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2643030202u,        0.6238861431u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2643030202u,        0.6238861431u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3692725791u,        0.6313123674u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3692725791u,        0.6313123674u) withcolor MidnightBlue;
draw (        0.3571129903u,        0.7358405983u)--(        0.2822343156u,        0.8097793647u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1773246241u,        0.8015505302u)--(        0.2822343156u,        0.8097793647u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3571129903u,        0.7358405983u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3571129903u,        0.7358405983u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1773246241u,        0.8015505302u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1773246241u,        0.8015505302u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2822343156u,        0.8097793647u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2822343156u,        0.8097793647u) withcolor MidnightBlue;
draw (        0.3571129903u,        0.7358405983u)--(        0.3692724953u,        0.6313135514u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2643030202u,        0.6238861431u)--(        0.3692724953u,        0.6313135514u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3571129903u,        0.7358405983u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3571129903u,        0.7358405983u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2643030202u,        0.6238861431u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2643030202u,        0.6238861431u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3692724953u,        0.6313135514u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3692724953u,        0.6313135514u) withcolor MidnightBlue;
draw (        0.1773246241u,        0.8015505302u)--(        0.1885394555u,        0.6969179127u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2643030202u,        0.6238861431u)--(        0.1885394555u,        0.6969179127u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1773246241u,        0.8015505302u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1773246241u,        0.8015505302u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2643030202u,        0.6238861431u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2643030202u,        0.6238861431u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1885394555u,        0.6969179127u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1885394555u,        0.6969179127u) withcolor MidnightBlue;
draw (        0.1773246241u,        0.8015505302u)--(        0.1885385848u,        0.6969178194u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0850045740u,        0.7157450790u)--(        0.1885385848u,        0.6969178194u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1773246241u,        0.8015505302u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1773246241u,        0.8015505302u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0850045740u,        0.7157450790u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0850045740u,        0.7157450790u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1885385848u,        0.6969178194u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1885385848u,        0.6969178194u) withcolor MidnightBlue;
draw (        0.1773246241u,        0.8015505302u)--(        0.0737906133u,        0.8203777898u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0850045740u,        0.7157450790u)--(        0.0737906133u,        0.8203777898u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1773246241u,        0.8015505302u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1773246241u,        0.8015505302u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0850045740u,        0.7157450790u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0850045740u,        0.7157450790u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0737906133u,        0.8203777898u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0737906133u,        0.8203777898u) withcolor Green;
draw (        0.7157632804u,        0.0850205270u)--(        0.6111111111u,        0.0962250449u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7157632804u,        0.0850205270u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7157632804u,        0.0850205270u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6111111111u,        0.0962250449u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6111111111u,        0.0962250449u) withcolor Yellow;
draw (        0.7157632804u,        0.0850205270u)--(        0.7777777778u,        0.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7157632804u,        0.0850205270u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7157632804u,        0.0850205270u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7777777778u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7777777778u,        0.0000000000u) withcolor Yellow;
draw (        0.3767441149u,        0.2577806084u)--(        0.3947676628u,        0.3614575528u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3299092410u,        0.4443257748u)--(        0.3947676628u,        0.3614575528u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3767441149u,        0.2577806084u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3767441149u,        0.2577806084u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3299092410u,        0.4443257748u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3299092410u,        0.4443257748u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3947676628u,        0.3614575528u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3947676628u,        0.3614575528u) withcolor MidnightBlue;
draw (        0.3767441149u,        0.2577806084u)--(        0.2966569232u,        0.1895178892u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1914275902u,        0.1902558381u)--(        0.2966569232u,        0.1895178892u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3767441149u,        0.2577806084u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3767441149u,        0.2577806084u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1914275902u,        0.1902558381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1914275902u,        0.1902558381u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2966569232u,        0.1895178892u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2966569232u,        0.1895178892u) withcolor MidnightBlue;
draw (        0.3767441149u,        0.2577806084u)--(        0.2715147818u,        0.2585185572u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1914275902u,        0.1902558381u)--(        0.2715147818u,        0.2585185572u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3767441149u,        0.2577806084u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3767441149u,        0.2577806084u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1914275902u,        0.1902558381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1914275902u,        0.1902558381u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2715147818u,        0.2585185572u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2715147818u,        0.2585185572u) withcolor MidnightBlue;
draw (        0.3767441149u,        0.2577806084u)--(        0.2715148495u,        0.2585281463u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2490332076u,        0.3613305435u)--(        0.2715148495u,        0.2585281463u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3767441149u,        0.2577806084u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3767441149u,        0.2577806084u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2490332076u,        0.3613305435u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2490332076u,        0.3613305435u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2715148495u,        0.2585281463u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2715148495u,        0.2585281463u) withcolor MidnightBlue;
draw (        0.3767441149u,        0.2577806084u)--(        0.2966578448u,        0.1895168080u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2842367195u,        0.0850205270u)--(        0.2966578448u,        0.1895168080u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3767441149u,        0.2577806084u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3767441149u,        0.2577806084u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2842367195u,        0.0850205270u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2842367195u,        0.0850205270u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2966578448u,        0.1895168080u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2966578448u,        0.1895168080u) withcolor MidnightBlue;
draw (        0.3299092410u,        0.4443257748u)--(        0.2557651546u,        0.5190011492u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2643030202u,        0.6238861431u)--(        0.2557651546u,        0.5190011492u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3299092410u,        0.4443257748u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3299092410u,        0.4443257748u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2643030202u,        0.6238861431u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2643030202u,        0.6238861431u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2557651546u,        0.5190011492u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2557651546u,        0.5190011492u) withcolor Green;
draw (        0.3299092410u,        0.4443257748u)--(        0.2265585124u,        0.4641344596u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2490332076u,        0.3613305435u)--(        0.2265585124u,        0.4641344596u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3299092410u,        0.4443257748u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3299092410u,        0.4443257748u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2490332076u,        0.3613305435u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2490332076u,        0.3613305435u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2265585124u,        0.4641344596u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2265585124u,        0.4641344596u) withcolor Green;
draw (        0.2643030202u,        0.6238861431u)--(        0.1885384033u,        0.6969168211u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0850045740u,        0.7157450790u)--(        0.1885384033u,        0.6969168211u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2643030202u,        0.6238861431u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2643030202u,        0.6238861431u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0850045740u,        0.7157450790u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0850045740u,        0.7157450790u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1885384033u,        0.6969168211u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1885384033u,        0.6969168211u) withcolor MidnightBlue;
draw (        0.0850045740u,        0.7157450790u)--(        0.0000000000u,        0.7777777778u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0850045740u,        0.7157450790u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0850045740u,        0.7157450790u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.7777777778u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.7777777778u) withcolor Yellow;
draw (        0.0850045740u,        0.7157450790u)--(        0.0962250449u,        0.6111111111u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0850045740u,        0.7157450790u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0850045740u,        0.7157450790u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0962250449u,        0.6111111111u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0962250449u,        0.6111111111u) withcolor Yellow;
draw (        0.1914275902u,        0.1902558381u)--(        0.2715079026u,        0.2585266273u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2490332076u,        0.3613305435u)--(        0.2715079026u,        0.2585266273u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1914275902u,        0.1902558381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1914275902u,        0.1902558381u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2490332076u,        0.3613305435u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2490332076u,        0.3613305435u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2715079026u,        0.2585266273u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2715079026u,        0.2585266273u) withcolor MidnightBlue;
draw (        0.1914275902u,        0.1902558381u)--(        0.2966569164u,        0.1895169184u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2842367195u,        0.0850205270u)--(        0.2966569164u,        0.1895169184u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1914275902u,        0.1902558381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1914275902u,        0.1902558381u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2842367195u,        0.0850205270u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2842367195u,        0.0850205270u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2966569164u,        0.1895169184u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2966569164u,        0.1895169184u) withcolor MidnightBlue;
draw (        0.1914275902u,        0.1902558381u)--(        0.1790073933u,        0.0857594467u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2842367195u,        0.0850205270u)--(        0.1790073933u,        0.0857594467u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1914275902u,        0.1902558381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1914275902u,        0.1902558381u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2842367195u,        0.0850205270u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2842367195u,        0.0850205270u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1790073933u,        0.0857594467u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1790073933u,        0.0857594467u) withcolor MidnightBlue;
draw (        0.1914275902u,        0.1902558381u)--(        0.0867972717u,        0.1790195778u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0850210082u,        0.2842365060u)--(        0.0867972717u,        0.1790195778u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1914275902u,        0.1902558381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1914275902u,        0.1902558381u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0850210082u,        0.2842365060u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0850210082u,        0.2842365060u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0867972717u,        0.1790195778u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0867972717u,        0.1790195778u) withcolor MidnightBlue;
draw (        0.1914275902u,        0.1902558381u)--(        0.1790014158u,        0.0857601574u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0743719922u,        0.0745155666u)--(        0.1790014158u,        0.0857601574u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1914275902u,        0.1902558381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1914275902u,        0.1902558381u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0743719922u,        0.0745155666u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0743719922u,        0.0745155666u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1790014158u,        0.0857601574u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1790014158u,        0.0857601574u) withcolor MidnightBlue;
draw (        0.1914275902u,        0.1902558381u)--(        0.0867981666u,        0.1790112474u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0743719922u,        0.0745155666u)--(        0.0867981666u,        0.1790112474u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1914275902u,        0.1902558381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1914275902u,        0.1902558381u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0743719922u,        0.0745155666u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0743719922u,        0.0745155666u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0867981666u,        0.1790112474u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0867981666u,        0.1790112474u) withcolor MidnightBlue;
draw (        0.2842367195u,        0.0850205270u)--(        0.2222222222u,        0.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2842367195u,        0.0850205270u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2842367195u,        0.0850205270u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2222222222u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2222222222u,        0.0000000000u) withcolor Yellow;
draw (        0.2842367195u,        0.0850205270u)--(        0.3888888889u,        0.0962250449u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2842367195u,        0.0850205270u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2842367195u,        0.0850205270u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3888888889u,        0.0962250449u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3888888889u,        0.0962250449u) withcolor Yellow;
draw (        0.2842367195u,        0.0850205270u)--(        0.1796016692u,        0.0738284163u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0743719922u,        0.0745155666u)--(        0.1796016692u,        0.0738284163u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2842367195u,        0.0850205270u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2842367195u,        0.0850205270u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0743719922u,        0.0745155666u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0743719922u,        0.0745155666u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1796016692u,        0.0738284163u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1796016692u,        0.0738284163u) withcolor Green;
draw (        0.2842367195u,        0.0850205270u)--(        0.1790070425u,        0.0857076773u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0743719922u,        0.0745155666u)--(        0.1790070425u,        0.0857076773u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2842367195u,        0.0850205270u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2842367195u,        0.0850205270u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0743719922u,        0.0745155666u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0743719922u,        0.0745155666u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1790070425u,        0.0857076773u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1790070425u,        0.0857076773u) withcolor MidnightBlue;
draw (        0.0850210082u,        0.2842365060u)--(        0.0962250449u,        0.3888888889u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0850210082u,        0.2842365060u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0850210082u,        0.2842365060u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0962250449u,        0.3888888889u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0962250449u,        0.3888888889u) withcolor Yellow;
draw (        0.0850210082u,        0.2842365060u)--(        0.0000000000u,        0.2222222222u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0850210082u,        0.2842365060u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0850210082u,        0.2842365060u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.2222222222u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.2222222222u) withcolor Yellow;
draw (        0.0850210082u,        0.2842365060u)--(        0.0867364455u,        0.1790185685u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0743719922u,        0.0745155666u)--(        0.0867364455u,        0.1790185685u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0850210082u,        0.2842365060u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0850210082u,        0.2842365060u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0743719922u,        0.0745155666u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0743719922u,        0.0745155666u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0867364455u,        0.1790185685u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0867364455u,        0.1790185685u) withcolor MidnightBlue;
draw (        0.0743719922u,        0.0745155666u)--(        0.0000000000u,        0.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0743719922u,        0.0745155666u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0743719922u,        0.0745155666u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor Yellow;
% r =         0.1052319205
endfig;
end;
