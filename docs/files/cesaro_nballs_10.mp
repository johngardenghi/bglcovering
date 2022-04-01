prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(10);
u := 20cm;
path pols[];
pair p[];
R := .2127686460u;
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

x1 = .9036525769u;
y1 = .8214460753u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .6399817856u;
y2 = .5075920159u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .6450449914u;
y3 = .8369312544u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .8821822380u;
y4 = .1772017865u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .6365034200u;
y5 = .1700875822u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .3600182144u;
y6 = .4924079842u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .3634965800u;
y7 = .8299124179u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .1178177621u;
y8 = .8227982134u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .3549550086u;
y9 = .1630687456u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = .0937650637u;
y10 = .1792387800u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.6399817856u,        0.5075920159u)--(        0.5086726292u,        0.3401751645u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6365034200u,        0.1700875822u)--(        0.5086726292u,        0.3401751645u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6399817856u,        0.5075920159u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6399817856u,        0.5075920159u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6365034200u,        0.1700875822u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6365034200u,        0.1700875822u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5086726292u,        0.3401751645u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5086726292u,        0.3401751645u) withcolor MidnightBlue;
draw (        0.6399817856u,        0.5075920159u)--(        0.5086680224u,        0.3401787778u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3600182144u,        0.4924079842u)--(        0.5086680224u,        0.3401787778u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6399817856u,        0.5075920159u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6399817856u,        0.5075920159u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3600182144u,        0.4924079842u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3600182144u,        0.4924079842u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5086680224u,        0.3401787778u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5086680224u,        0.3401787778u) withcolor MidnightBlue;
draw (        0.6399817856u,        0.5075920159u)--(        0.4913319776u,        0.6598212223u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3600182144u,        0.4924079842u)--(        0.4913319776u,        0.6598212223u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6399817856u,        0.5075920159u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6399817856u,        0.5075920159u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3600182144u,        0.4924079842u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3600182144u,        0.4924079842u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4913319776u,        0.6598212223u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4913319776u,        0.6598212223u) withcolor MidnightBlue;
draw (        0.6399817856u,        0.5075920159u)--(        0.4913634901u,        0.6598519875u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3634965800u,        0.8299124179u)--(        0.4913634901u,        0.6598519875u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6399817856u,        0.5075920159u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6399817856u,        0.5075920159u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3634965800u,        0.8299124179u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3634965800u,        0.8299124179u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4913634901u,        0.6598519875u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4913634901u,        0.6598519875u) withcolor MidnightBlue;
draw (        0.8821822380u,        0.1772017865u)--(        1.0000000000u,        0.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8821822380u,        0.1772017865u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8821822380u,        0.1772017865u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.8821822380u,        0.1772017865u)--(        0.9037749551u,        0.3888888889u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8821822380u,        0.1772017865u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8821822380u,        0.1772017865u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9037749551u,        0.3888888889u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9037749551u,        0.3888888889u) withcolor Yellow;
draw (        0.6365034200u,        0.1700875822u)--(        0.5086365100u,        0.3401480127u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3600182144u,        0.4924079842u)--(        0.5086365100u,        0.3401480127u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6365034200u,        0.1700875822u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6365034200u,        0.1700875822u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3600182144u,        0.4924079842u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3600182144u,        0.4924079842u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5086365100u,        0.3401480127u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5086365100u,        0.3401480127u) withcolor MidnightBlue;
draw (        0.3600182144u,        0.4924079842u)--(        0.4913273708u,        0.6598248357u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3634965800u,        0.8299124179u)--(        0.4913273708u,        0.6598248357u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3600182144u,        0.4924079842u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3600182144u,        0.4924079842u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3634965800u,        0.8299124179u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3634965800u,        0.8299124179u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4913273708u,        0.6598248357u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4913273708u,        0.6598248357u) withcolor MidnightBlue;
draw (        0.1178177621u,        0.8227982134u)--(        0.0000000000u,        1.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1178177621u,        0.8227982134u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1178177621u,        0.8227982134u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        1.0000000000u) withcolor Yellow;
draw (        0.1178177621u,        0.8227982134u)--(        0.0962250449u,        0.6111111111u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1178177621u,        0.8227982134u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1178177621u,        0.8227982134u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0962250449u,        0.6111111111u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0962250449u,        0.6111111111u) withcolor Yellow;
% r =         0.2127686460
endfig;
end;
