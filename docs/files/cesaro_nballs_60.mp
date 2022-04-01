prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(60);
u := 20cm;
path pols[];
pair p[];
R := .0759122527u;
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

x1 = .9421617706u;
y1 = .9508102737u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .9647979534u;
y2 = .8343393648u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .8345069379u;
y3 = .8936964995u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .7242519573u;
y4 = .9648765018u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .9626617137u;
y5 = .7009605849u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .8518610433u;
y6 = .7749028581u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .6870218072u;
y7 = .9030213527u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .8641439517u;
y8 = .6538774276u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .7223610789u;
y9 = .8148557321u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = .7373136762u;
y10 = .7017440596u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;

x11 = .6043023540u;
y11 = .7493539941u;
draw fullcircle scaled D shifted z11 withcolor RoyalBlue withpen pencircle scaled 1.0;

x12 = .2762629011u;
y12 = .9648172966u;
draw fullcircle scaled D shifted z12 withcolor RoyalBlue withpen pencircle scaled 1.0;

x13 = .7390748433u;
y13 = .5883295786u;
draw fullcircle scaled D shifted z13 withcolor RoyalBlue withpen pencircle scaled 1.0;

x14 = .6174635186u;
y14 = .6374726157u;
draw fullcircle scaled D shifted z14 withcolor RoyalBlue withpen pencircle scaled 1.0;

x15 = .5076600968u;
y15 = .6897577252u;
draw fullcircle scaled D shifted z15 withcolor RoyalBlue withpen pencircle scaled 1.0;

x16 = .3994628086u;
y16 = .7463964677u;
draw fullcircle scaled D shifted z16 withcolor RoyalBlue withpen pencircle scaled 1.0;

x17 = .2786014891u;
y17 = .8054877169u;
draw fullcircle scaled D shifted z17 withcolor RoyalBlue withpen pencircle scaled 1.0;

x18 = .3139801932u;
y18 = .9018090336u;
draw fullcircle scaled D shifted z18 withcolor RoyalBlue withpen pencircle scaled 1.0;

x19 = .1704086046u;
y19 = .8924526367u;
draw fullcircle scaled D shifted z19 withcolor RoyalBlue withpen pencircle scaled 1.0;

x20 = .0589860212u;
y20 = .9522060002u;
draw fullcircle scaled D shifted z20 withcolor RoyalBlue withpen pencircle scaled 1.0;

x21 = .9538832839u;
y21 = .2730338845u;
draw fullcircle scaled D shifted z21 withcolor RoyalBlue withpen pencircle scaled 1.0;

x22 = .9131685244u;
y22 = .3150187097u;
draw fullcircle scaled D shifted z22 withcolor RoyalBlue withpen pencircle scaled 1.0;

x23 = .7932563251u;
y23 = .3762646768u;
draw fullcircle scaled D shifted z23 withcolor RoyalBlue withpen pencircle scaled 1.0;

x24 = .7183144197u;
y24 = .4554704764u;
draw fullcircle scaled D shifted z24 withcolor RoyalBlue withpen pencircle scaled 1.0;

x25 = .6187378042u;
y25 = .5164618138u;
draw fullcircle scaled D shifted z25 withcolor RoyalBlue withpen pencircle scaled 1.0;

x26 = .4978230044u;
y26 = .5615893438u;
draw fullcircle scaled D shifted z26 withcolor RoyalBlue withpen pencircle scaled 1.0;

x27 = .3903736952u;
y27 = .6150605422u;
draw fullcircle scaled D shifted z27 withcolor RoyalBlue withpen pencircle scaled 1.0;

x28 = .2813097361u;
y28 = .6812722800u;
draw fullcircle scaled D shifted z28 withcolor RoyalBlue withpen pencircle scaled 1.0;

x29 = .1636398933u;
y29 = .7615006227u;
draw fullcircle scaled D shifted z29 withcolor RoyalBlue withpen pencircle scaled 1.0;

x30 = .0573993480u;
y30 = .8547188200u;
draw fullcircle scaled D shifted z30 withcolor RoyalBlue withpen pencircle scaled 1.0;

x31 = .9426006519u;
y31 = .1452811799u;
draw fullcircle scaled D shifted z31 withcolor RoyalBlue withpen pencircle scaled 1.0;

x32 = .8363601066u;
y32 = .2384993773u;
draw fullcircle scaled D shifted z32 withcolor RoyalBlue withpen pencircle scaled 1.0;

x33 = .7186902638u;
y33 = .3187277201u;
draw fullcircle scaled D shifted z33 withcolor RoyalBlue withpen pencircle scaled 1.0;

x34 = .6096263047u;
y34 = .3849394578u;
draw fullcircle scaled D shifted z34 withcolor RoyalBlue withpen pencircle scaled 1.0;

x35 = .5021769957u;
y35 = .4384106561u;
draw fullcircle scaled D shifted z35 withcolor RoyalBlue withpen pencircle scaled 1.0;

x36 = .3812621958u;
y36 = .4835381863u;
draw fullcircle scaled D shifted z36 withcolor RoyalBlue withpen pencircle scaled 1.0;

x37 = .2816855803u;
y37 = .5445295235u;
draw fullcircle scaled D shifted z37 withcolor RoyalBlue withpen pencircle scaled 1.0;

x38 = .2067436749u;
y38 = .6237353231u;
draw fullcircle scaled D shifted z38 withcolor RoyalBlue withpen pencircle scaled 1.0;

x39 = .0886155066u;
y39 = .6729408046u;
draw fullcircle scaled D shifted z39 withcolor RoyalBlue withpen pencircle scaled 1.0;

x40 = .0461167160u;
y40 = .7269661156u;
draw fullcircle scaled D shifted z40 withcolor RoyalBlue withpen pencircle scaled 1.0;

x41 = .9410139788u;
y41 = .0477939998u;
draw fullcircle scaled D shifted z41 withcolor RoyalBlue withpen pencircle scaled 1.0;

x42 = .8295913953u;
y42 = .1075473633u;
draw fullcircle scaled D shifted z42 withcolor RoyalBlue withpen pencircle scaled 1.0;

x43 = .6860198068u;
y43 = .0981909664u;
draw fullcircle scaled D shifted z43 withcolor RoyalBlue withpen pencircle scaled 1.0;

x44 = .7213985109u;
y44 = .1945122832u;
draw fullcircle scaled D shifted z44 withcolor RoyalBlue withpen pencircle scaled 1.0;

x45 = .6005371915u;
y45 = .2536035324u;
draw fullcircle scaled D shifted z45 withcolor RoyalBlue withpen pencircle scaled 1.0;

x46 = .4923399034u;
y46 = .3102422747u;
draw fullcircle scaled D shifted z46 withcolor RoyalBlue withpen pencircle scaled 1.0;

x47 = .3825364813u;
y47 = .3625273842u;
draw fullcircle scaled D shifted z47 withcolor RoyalBlue withpen pencircle scaled 1.0;

x48 = .2609251568u;
y48 = .4116704214u;
draw fullcircle scaled D shifted z48 withcolor RoyalBlue withpen pencircle scaled 1.0;

x49 = .7237370989u;
y49 = .0351827035u;
draw fullcircle scaled D shifted z49 withcolor RoyalBlue withpen pencircle scaled 1.0;

x50 = .3956976459u;
y50 = .2506460059u;
draw fullcircle scaled D shifted z50 withcolor RoyalBlue withpen pencircle scaled 1.0;

x51 = .2626863237u;
y51 = .2982559405u;
draw fullcircle scaled D shifted z51 withcolor RoyalBlue withpen pencircle scaled 1.0;

x52 = .2776389212u;
y52 = .1851442680u;
draw fullcircle scaled D shifted z52 withcolor RoyalBlue withpen pencircle scaled 1.0;

x53 = .1358560484u;
y53 = .3461225723u;
draw fullcircle scaled D shifted z53 withcolor RoyalBlue withpen pencircle scaled 1.0;

x54 = .3129781928u;
y54 = .0969786473u;
draw fullcircle scaled D shifted z54 withcolor RoyalBlue withpen pencircle scaled 1.0;

x55 = .1481389566u;
y55 = .2250971421u;
draw fullcircle scaled D shifted z55 withcolor RoyalBlue withpen pencircle scaled 1.0;

x56 = .0373382862u;
y56 = .2990394151u;
draw fullcircle scaled D shifted z56 withcolor RoyalBlue withpen pencircle scaled 1.0;

x57 = .2757480427u;
y57 = .0351234982u;
draw fullcircle scaled D shifted z57 withcolor RoyalBlue withpen pencircle scaled 1.0;

x58 = .1654930621u;
y58 = .1063035006u;
draw fullcircle scaled D shifted z58 withcolor RoyalBlue withpen pencircle scaled 1.0;

x59 = .0352020467u;
y59 = .1656606351u;
draw fullcircle scaled D shifted z59 withcolor RoyalBlue withpen pencircle scaled 1.0;

x60 = .0578382294u;
y60 = .0491897262u;
draw fullcircle scaled D shifted z60 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.9421617706u,        0.9508102737u)--(        1.0000000000u,        1.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9421617706u,        0.9508102737u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9421617706u,        0.9508102737u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        1.0000000000u) withcolor Yellow;
draw (        0.9647979534u,        0.8343393648u)--(        0.8891814494u,        0.8410337240u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8345069379u,        0.8936964995u)--(        0.8891814494u,        0.8410337240u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9647979534u,        0.8343393648u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9647979534u,        0.8343393648u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8345069379u,        0.8936964995u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8345069379u,        0.8936964995u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8891814494u,        0.8410337240u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8891814494u,        0.8410337240u) withcolor MidnightBlue;
draw (        0.9647979534u,        0.8343393648u)--(        0.9274845924u,        0.7682304909u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9626617137u,        0.7009605849u)--(        0.9274845924u,        0.7682304909u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9647979534u,        0.8343393648u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9647979534u,        0.8343393648u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9626617137u,        0.7009605849u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9626617137u,        0.7009605849u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9274845924u,        0.7682304909u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9274845924u,        0.7682304909u) withcolor MidnightBlue;
draw (        0.9647979534u,        0.8343393648u)--(        0.9274797330u,        0.7682332339u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8518610433u,        0.7749028581u)--(        0.9274797330u,        0.7682332339u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9647979534u,        0.8343393648u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9647979534u,        0.8343393648u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8518610433u,        0.7749028581u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8518610433u,        0.7749028581u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9274797330u,        0.7682332339u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9274797330u,        0.7682332339u) withcolor MidnightBlue;
draw (        0.9647979534u,        0.8343393648u)--(        0.8891792637u,        0.8410089891u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8518610433u,        0.7749028581u)--(        0.8891792637u,        0.8410089891u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9647979534u,        0.8343393648u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9647979534u,        0.8343393648u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8518610433u,        0.7749028581u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8518610433u,        0.7749028581u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8891792637u,        0.8410089891u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8891792637u,        0.8410089891u) withcolor MidnightBlue;
draw (        0.8345069379u,        0.8936964995u)--(        0.8000824704u,        0.9613546561u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7242519573u,        0.9648765018u)--(        0.8000824704u,        0.9613546561u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8345069379u,        0.8936964995u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8345069379u,        0.8936964995u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7242519573u,        0.9648765018u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7242519573u,        0.9648765018u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8000824704u,        0.9613546561u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8000824704u,        0.9613546561u) withcolor Green;
draw (        0.8345069379u,        0.8936964995u)--(        0.8891652419u,        0.8410169027u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8518610433u,        0.7749028581u)--(        0.8891652419u,        0.8410169027u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8345069379u,        0.8936964995u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8345069379u,        0.8936964995u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8518610433u,        0.7749028581u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8518610433u,        0.7749028581u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8891652419u,        0.8410169027u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8891652419u,        0.8410169027u) withcolor MidnightBlue;
draw (        0.8345069379u,        0.8936964995u)--(        0.7972027393u,        0.8275824549u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8518610433u,        0.7749028581u)--(        0.7972027393u,        0.8275824549u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8345069379u,        0.8936964995u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8345069379u,        0.8936964995u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8518610433u,        0.7749028581u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8518610433u,        0.7749028581u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7972027393u,        0.8275824549u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7972027393u,        0.8275824549u) withcolor MidnightBlue;
draw (        0.8345069379u,        0.8936964995u)--(        0.7596660640u,        0.8809876951u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6870218072u,        0.9030213527u)--(        0.7596660640u,        0.8809876951u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8345069379u,        0.8936964995u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8345069379u,        0.8936964995u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6870218072u,        0.9030213527u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6870218072u,        0.9030213527u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7596660640u,        0.8809876951u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7596660640u,        0.8809876951u) withcolor MidnightBlue;
draw (        0.8345069379u,        0.8936964995u)--(        0.7971985011u,        0.8275848465u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7223610789u,        0.8148557321u)--(        0.7971985011u,        0.8275848465u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8345069379u,        0.8936964995u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8345069379u,        0.8936964995u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7223610789u,        0.8148557321u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7223610789u,        0.8148557321u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7971985011u,        0.8275848465u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7971985011u,        0.8275848465u) withcolor MidnightBlue;
draw (        0.8345069379u,        0.8936964995u)--(        0.7596695157u,        0.8809673851u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7223610789u,        0.8148557321u)--(        0.7596695157u,        0.8809673851u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8345069379u,        0.8936964995u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8345069379u,        0.8936964995u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7223610789u,        0.8148557321u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7223610789u,        0.8148557321u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7596695157u,        0.8809673851u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7596695157u,        0.8809673851u) withcolor MidnightBlue;
draw (        0.7242519573u,        0.9648765018u)--(        0.6484209084u,        0.9683867925u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6870218072u,        0.9030213527u)--(        0.6484209084u,        0.9683867925u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7242519573u,        0.9648765018u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7242519573u,        0.9648765018u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6870218072u,        0.9030213527u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6870218072u,        0.9030213527u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6484209084u,        0.9683867925u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6484209084u,        0.9683867925u) withcolor Green;
draw (        0.9626617137u,        0.7009605849u)--(        0.9274792442u,        0.7682276940u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8518610433u,        0.7749028581u)--(        0.9274792442u,        0.7682276940u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9626617137u,        0.7009605849u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9626617137u,        0.7009605849u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8518610433u,        0.7749028581u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8518610433u,        0.7749028581u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9274792442u,        0.7682276940u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9274792442u,        0.7682276940u) withcolor MidnightBlue;
draw (        0.8518610433u,        0.7749028581u)--(        0.8128125713u,        0.7098038033u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8641439517u,        0.6538774276u)--(        0.8128125713u,        0.7098038033u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8518610433u,        0.7749028581u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8518610433u,        0.7749028581u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8641439517u,        0.6538774276u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8641439517u,        0.6538774276u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8128125713u,        0.7098038033u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8128125713u,        0.7098038033u) withcolor MidnightBlue;
draw (        0.8518610433u,        0.7749028581u)--(        0.7971994825u,        0.8275790756u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7223610789u,        0.8148557321u)--(        0.7971994825u,        0.8275790756u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8518610433u,        0.7749028581u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8518610433u,        0.7749028581u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7223610789u,        0.8148557321u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7223610789u,        0.8148557321u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7971994825u,        0.8275790756u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7971994825u,        0.8275790756u) withcolor MidnightBlue;
draw (        0.8518610433u,        0.7749028581u)--(        0.8127957811u,        0.7098138775u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7373136762u,        0.7017440596u)--(        0.8127957811u,        0.7098138775u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8518610433u,        0.7749028581u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8518610433u,        0.7749028581u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7373136762u,        0.7017440596u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7373136762u,        0.7017440596u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8127957811u,        0.7098138775u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8127957811u,        0.7098138775u) withcolor MidnightBlue;
draw (        0.6870218072u,        0.9030213527u)--(        0.6111111111u,        0.9037749551u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6870218072u,        0.9030213527u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6870218072u,        0.9030213527u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6111111111u,        0.9037749551u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6111111111u,        0.9037749551u) withcolor Yellow;
draw (        0.6870218072u,        0.9030213527u)--(        0.7596613069u,        0.8809720169u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7223610789u,        0.8148557321u)--(        0.7596613069u,        0.8809720169u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6870218072u,        0.9030213527u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6870218072u,        0.9030213527u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7223610789u,        0.8148557321u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7223610789u,        0.8148557321u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7596613069u,        0.8809720169u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7596613069u,        0.8809720169u) withcolor MidnightBlue;
draw (        0.6870218072u,        0.9030213527u)--(        0.6497215793u,        0.8369050679u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7223610789u,        0.8148557321u)--(        0.6497215793u,        0.8369050679u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6870218072u,        0.9030213527u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6870218072u,        0.9030213527u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7223610789u,        0.8148557321u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7223610789u,        0.8148557321u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6497215793u,        0.8369050679u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6497215793u,        0.8369050679u) withcolor Green;
draw (        0.8641439517u,        0.6538774276u)--(        0.7886593687u,        0.6458308234u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7373136762u,        0.7017440596u)--(        0.7886593687u,        0.6458308234u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8641439517u,        0.6538774276u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8641439517u,        0.6538774276u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7373136762u,        0.7017440596u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7373136762u,        0.7017440596u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7886593687u,        0.6458308234u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7886593687u,        0.6458308234u) withcolor MidnightBlue;
draw (        0.8641439517u,        0.6538774276u)--(        0.8127982593u,        0.7097906638u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7373136762u,        0.7017440596u)--(        0.8127982593u,        0.7097906638u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8641439517u,        0.6538774276u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8641439517u,        0.6538774276u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7373136762u,        0.7017440596u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7373136762u,        0.7017440596u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8127982593u,        0.7097906638u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8127982593u,        0.7097906638u) withcolor MidnightBlue;
draw (        0.8641439517u,        0.6538774276u)--(        0.7886617699u,        0.6458083297u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7390748433u,        0.5883295786u)--(        0.7886617699u,        0.6458083297u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8641439517u,        0.6538774276u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8641439517u,        0.6538774276u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7390748433u,        0.5883295786u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7390748433u,        0.5883295786u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7886617699u,        0.6458083297u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7886617699u,        0.6458083297u) withcolor MidnightBlue;
draw (        0.7223610789u,        0.8148557321u)--(        0.6801872606u,        0.7517364865u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7373136762u,        0.7017440596u)--(        0.6801872606u,        0.7517364865u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7223610789u,        0.8148557321u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7223610789u,        0.8148557321u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7373136762u,        0.7017440596u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7373136762u,        0.7017440596u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6801872606u,        0.7517364865u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6801872606u,        0.7517364865u) withcolor MidnightBlue;
draw (        0.7223610789u,        0.8148557321u)--(        0.6801769947u,        0.7517433469u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6043023540u,        0.7493539941u)--(        0.6801769947u,        0.7517433469u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7223610789u,        0.8148557321u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7223610789u,        0.8148557321u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6043023540u,        0.7493539941u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6043023540u,        0.7493539941u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6801769947u,        0.7517433469u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6801769947u,        0.7517433469u) withcolor MidnightBlue;
draw (        0.7373136762u,        0.7017440596u)--(        0.6614384723u,        0.6993726566u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6043023540u,        0.7493539941u)--(        0.6614384723u,        0.6993726566u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7373136762u,        0.7017440596u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7373136762u,        0.7017440596u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6043023540u,        0.7493539941u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6043023540u,        0.7493539941u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6614384723u,        0.6993726566u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6614384723u,        0.6993726566u) withcolor MidnightBlue;
draw (        0.7373136762u,        0.7017440596u)--(        0.6801775579u,        0.7517253971u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6043023540u,        0.7493539941u)--(        0.6801775579u,        0.7517253971u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7373136762u,        0.7017440596u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7373136762u,        0.7017440596u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6043023540u,        0.7493539941u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6043023540u,        0.7493539941u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6801775579u,        0.7517253971u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6801775579u,        0.7517253971u) withcolor MidnightBlue;
draw (        0.7373136762u,        0.7017440596u)--(        0.7886478992u,        0.6458202930u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7390748433u,        0.5883295786u)--(        0.7886478992u,        0.6458202930u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7373136762u,        0.7017440596u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7373136762u,        0.7017440596u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7390748433u,        0.5883295786u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7390748433u,        0.5883295786u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7886478992u,        0.6458202930u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7886478992u,        0.6458202930u) withcolor MidnightBlue;
draw (        0.7373136762u,        0.7017440596u)--(        0.6614391837u,        0.6993500025u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6174635186u,        0.6374726157u)--(        0.6614391837u,        0.6993500025u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7373136762u,        0.7017440596u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7373136762u,        0.7017440596u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6174635186u,        0.6374726157u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6174635186u,        0.6374726157u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6614391837u,        0.6993500025u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6614391837u,        0.6993500025u) withcolor MidnightBlue;
draw (        0.6043023540u,        0.7493539941u)--(        0.5555555556u,        0.8075499103u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6043023540u,        0.7493539941u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6043023540u,        0.7493539941u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5555555556u,        0.8075499103u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5555555556u,        0.8075499103u) withcolor Yellow;
draw (        0.6043023540u,        0.7493539941u)--(        0.6614265221u,        0.6993589990u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6174635186u,        0.6374726157u)--(        0.6614265221u,        0.6993589990u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6043023540u,        0.7493539941u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6043023540u,        0.7493539941u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6174635186u,        0.6374726157u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6174635186u,        0.6374726157u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6614265221u,        0.6993589990u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6614265221u,        0.6993589990u) withcolor MidnightBlue;
draw (        0.6043023540u,        0.7493539941u)--(        0.5295285906u,        0.7624518728u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5076600968u,        0.6897577252u)--(        0.5295285906u,        0.7624518728u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6043023540u,        0.7493539941u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6043023540u,        0.7493539941u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5076600968u,        0.6897577252u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5076600968u,        0.6897577252u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5295285906u,        0.7624518728u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5295285906u,        0.7624518728u) withcolor Green;
draw (        0.2762629011u,        0.9648172966u)--(        0.3521298743u,        0.9674388435u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3139801932u,        0.9018090336u)--(        0.3521298743u,        0.9674388435u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2762629011u,        0.9648172966u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2762629011u,        0.9648172966u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3139801932u,        0.9018090336u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3139801932u,        0.9018090336u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3521298743u,        0.9674388435u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3521298743u,        0.9674388435u) withcolor Green;
draw (        0.2762629011u,        0.9648172966u)--(        0.2003960970u,        0.9621908645u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1704086046u,        0.8924526367u)--(        0.2003960970u,        0.9621908645u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2762629011u,        0.9648172966u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2762629011u,        0.9648172966u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1704086046u,        0.8924526367u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1704086046u,        0.8924526367u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2003960970u,        0.9621908645u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2003960970u,        0.9621908645u) withcolor Green;
draw (        0.7390748433u,        0.5883295786u)--(        0.8075499103u,        0.5555555556u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7390748433u,        0.5883295786u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7390748433u,        0.5883295786u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8075499103u,        0.5555555556u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8075499103u,        0.5555555556u) withcolor Yellow;
draw (        0.7390748433u,        0.5883295786u)--(        0.6639455071u,        0.5774551610u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6174635186u,        0.6374726157u)--(        0.6639455071u,        0.5774551610u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7390748433u,        0.5883295786u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7390748433u,        0.5883295786u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6174635186u,        0.6374726157u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6174635186u,        0.6374726157u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6639455071u,        0.5774551610u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6639455071u,        0.5774551610u) withcolor MidnightBlue;
draw (        0.7390748433u,        0.5883295786u)--(        0.6938746332u,        0.5273409638u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7183144197u,        0.4554704764u)--(        0.6938746332u,        0.5273409638u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7390748433u,        0.5883295786u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7390748433u,        0.5883295786u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7183144197u,        0.4554704764u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7183144197u,        0.4554704764u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6938746332u,        0.5273409638u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6938746332u,        0.5273409638u) withcolor MidnightBlue;
draw (        0.7390748433u,        0.5883295786u)--(        0.6938654715u,        0.5273477548u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6187378042u,        0.5164618138u)--(        0.6938654715u,        0.5273477548u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7390748433u,        0.5883295786u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7390748433u,        0.5883295786u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6187378042u,        0.5164618138u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6187378042u,        0.5164618138u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6938654715u,        0.5273477548u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6938654715u,        0.5273477548u) withcolor MidnightBlue;
draw (        0.7390748433u,        0.5883295786u)--(        0.6639471759u,        0.5774436375u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6187378042u,        0.5164618138u)--(        0.6639471759u,        0.5774436375u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7390748433u,        0.5883295786u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7390748433u,        0.5883295786u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6187378042u,        0.5164618138u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6187378042u,        0.5164618138u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6639471759u,        0.5774436375u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6639471759u,        0.5774436375u) withcolor MidnightBlue;
draw (        0.6174635186u,        0.6374726157u)--(        0.5430251450u,        0.6225864257u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5076600968u,        0.6897577252u)--(        0.5430251450u,        0.6225864257u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6174635186u,        0.6374726157u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6174635186u,        0.6374726157u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5076600968u,        0.6897577252u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5076600968u,        0.6897577252u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5430251450u,        0.6225864257u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5430251450u,        0.6225864257u) withcolor MidnightBlue;
draw (        0.6174635186u,        0.6374726157u)--(        0.6639387202u,        0.5774499053u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6187378042u,        0.5164618138u)--(        0.6639387202u,        0.5774499053u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6174635186u,        0.6374726157u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6174635186u,        0.6374726157u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6187378042u,        0.5164618138u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6187378042u,        0.5164618138u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6639387202u,        0.5774499053u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6639387202u,        0.5774499053u) withcolor MidnightBlue;
draw (        0.6174635186u,        0.6374726157u)--(        0.5722626026u,        0.5764845242u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6187378042u,        0.5164618138u)--(        0.5722626026u,        0.5764845242u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6174635186u,        0.6374726157u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6174635186u,        0.6374726157u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6187378042u,        0.5164618138u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6187378042u,        0.5164618138u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5722626026u,        0.5764845242u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5722626026u,        0.5764845242u) withcolor MidnightBlue;
draw (        0.6174635186u,        0.6374726157u)--(        0.5722590530u,        0.5764871551u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4978230044u,        0.5615893438u)--(        0.5722590530u,        0.5764871551u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6174635186u,        0.6374726157u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6174635186u,        0.6374726157u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4978230044u,        0.5615893438u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4978230044u,        0.5615893438u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5722590530u,        0.5764871551u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5722590530u,        0.5764871551u) withcolor MidnightBlue;
draw (        0.6174635186u,        0.6374726157u)--(        0.5430274700u,        0.6225748044u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4978230044u,        0.5615893438u)--(        0.5430274700u,        0.6225748044u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6174635186u,        0.6374726157u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6174635186u,        0.6374726157u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4978230044u,        0.5615893438u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4978230044u,        0.5615893438u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5430274700u,        0.6225748044u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5430274700u,        0.6225748044u) withcolor MidnightBlue;
draw (        0.5076600968u,        0.6897577252u)--(        0.4326451265u,        0.6781205283u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3994628086u,        0.7463964677u)--(        0.4326451265u,        0.6781205283u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5076600968u,        0.6897577252u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5076600968u,        0.6897577252u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3994628086u,        0.7463964677u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3994628086u,        0.7463964677u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4326451265u,        0.6781205283u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4326451265u,        0.6781205283u) withcolor MidnightBlue;
draw (        0.5076600968u,        0.6897577252u)--(        0.5430173359u,        0.6225823148u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4978230044u,        0.5615893438u)--(        0.5430173359u,        0.6225823148u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5076600968u,        0.6897577252u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5076600968u,        0.6897577252u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4978230044u,        0.5615893438u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4978230044u,        0.5615893438u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5430173359u,        0.6225823148u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5430173359u,        0.6225823148u) withcolor MidnightBlue;
draw (        0.5076600968u,        0.6897577252u)--(        0.4326461780u,        0.6781137524u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3903736952u,        0.6150605422u)--(        0.4326461780u,        0.6781137524u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5076600968u,        0.6897577252u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5076600968u,        0.6897577252u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3903736952u,        0.6150605422u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3903736952u,        0.6150605422u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4326461780u,        0.6781137524u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4326461780u,        0.6781137524u) withcolor MidnightBlue;
draw (        0.3994628086u,        0.7463964677u)--(        0.4444444444u,        0.8075499103u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3994628086u,        0.7463964677u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3994628086u,        0.7463964677u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4444444444u,        0.8075499103u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4444444444u,        0.8075499103u) withcolor Yellow;
draw (        0.3994628086u,        0.7463964677u)--(        0.3235785688u,        0.7443343653u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2786014891u,        0.8054877169u)--(        0.3235785688u,        0.7443343653u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3994628086u,        0.7463964677u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3994628086u,        0.7463964677u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2786014891u,        0.8054877169u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2786014891u,        0.8054877169u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3235785688u,        0.7443343653u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3235785688u,        0.7443343653u) withcolor MidnightBlue;
draw (        0.3994628086u,        0.7463964677u)--(        0.3544857289u,        0.8075498194u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2786014891u,        0.8054877169u)--(        0.3544857289u,        0.8075498194u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3994628086u,        0.7463964677u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3994628086u,        0.7463964677u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2786014891u,        0.8054877169u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2786014891u,        0.8054877169u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3544857289u,        0.8075498194u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3544857289u,        0.8075498194u) withcolor Green;
draw (        0.3994628086u,        0.7463964677u)--(        0.4326398761u,        0.6781179769u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3903736952u,        0.6150605422u)--(        0.4326398761u,        0.6781179769u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3994628086u,        0.7463964677u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3994628086u,        0.7463964677u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3903736952u,        0.6150605422u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3903736952u,        0.6150605422u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4326398761u,        0.6781179769u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4326398761u,        0.6781179769u) withcolor MidnightBlue;
draw (        0.3994628086u,        0.7463964677u)--(        0.3571966276u,        0.6833390330u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3903736952u,        0.6150605422u)--(        0.3571966276u,        0.6833390330u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3994628086u,        0.7463964677u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3994628086u,        0.7463964677u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3903736952u,        0.6150605422u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3903736952u,        0.6150605422u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3571966276u,        0.6833390330u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3571966276u,        0.6833390330u) withcolor MidnightBlue;
draw (        0.3994628086u,        0.7463964677u)--(        0.3571937977u,        0.6833409299u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2813097361u,        0.6812722800u)--(        0.3571937977u,        0.6833409299u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3994628086u,        0.7463964677u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3994628086u,        0.7463964677u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2813097361u,        0.6812722800u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2813097361u,        0.6812722800u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3571937977u,        0.6833409299u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3571937977u,        0.6833409299u) withcolor MidnightBlue;
draw (        0.3994628086u,        0.7463964677u)--(        0.3235787470u,        0.7443278178u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2813097361u,        0.6812722800u)--(        0.3235787470u,        0.7443278178u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3994628086u,        0.7463964677u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3994628086u,        0.7463964677u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2813097361u,        0.6812722800u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2813097361u,        0.6812722800u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3235787470u,        0.7443278178u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3235787470u,        0.7443278178u) withcolor MidnightBlue;
draw (        0.2786014891u,        0.8054877169u)--(        0.2437722431u,        0.8729383933u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3139801932u,        0.9018090336u)--(        0.2437722431u,        0.8729383933u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2786014891u,        0.8054877169u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2786014891u,        0.8054877169u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3139801932u,        0.9018090336u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3139801932u,        0.9018090336u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2437722431u,        0.8729383933u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2437722431u,        0.8729383933u) withcolor MidnightBlue;
draw (        0.2786014891u,        0.8054877169u)--(        0.3488094392u,        0.8343583572u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3139801932u,        0.9018090336u)--(        0.3488094392u,        0.8343583572u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2786014891u,        0.8054877169u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2786014891u,        0.8054877169u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3139801932u,        0.9018090336u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3139801932u,        0.9018090336u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3488094392u,        0.8343583572u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3488094392u,        0.8343583572u) withcolor Green;
draw (        0.2786014891u,        0.8054877169u)--(        0.2052406873u,        0.8250034252u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1704086046u,        0.8924526367u)--(        0.2052406873u,        0.8250034252u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2786014891u,        0.8054877169u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2786014891u,        0.8054877169u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1704086046u,        0.8924526367u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1704086046u,        0.8924526367u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2052406873u,        0.8250034252u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2052406873u,        0.8250034252u) withcolor MidnightBlue;
draw (        0.2786014891u,        0.8054877169u)--(        0.2437694065u,        0.8729369285u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1704086046u,        0.8924526367u)--(        0.2437694065u,        0.8729369285u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2786014891u,        0.8054877169u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2786014891u,        0.8054877169u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1704086046u,        0.8924526367u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1704086046u,        0.8924526367u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2437694065u,        0.8729369285u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2437694065u,        0.8729369285u) withcolor MidnightBlue;
draw (        0.2786014891u,        0.8054877169u)--(        0.3235739960u,        0.7443310023u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2813097361u,        0.6812722800u)--(        0.3235739960u,        0.7443310023u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2786014891u,        0.8054877169u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2786014891u,        0.8054877169u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2813097361u,        0.6812722800u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2813097361u,        0.6812722800u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3235739960u,        0.7443310023u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3235739960u,        0.7443310023u) withcolor MidnightBlue;
draw (        0.2786014891u,        0.8054877169u)--(        0.2052396851u,        0.8249996576u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1636398933u,        0.7615006227u)--(        0.2052396851u,        0.8249996576u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2786014891u,        0.8054877169u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2786014891u,        0.8054877169u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1636398933u,        0.7615006227u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1636398933u,        0.7615006227u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2052396851u,        0.8249996576u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2052396851u,        0.8249996576u) withcolor MidnightBlue;
draw (        0.3139801932u,        0.9018090336u)--(        0.2437707573u,        0.8729420069u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1704086046u,        0.8924526367u)--(        0.2437707573u,        0.8729420069u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3139801932u,        0.9018090336u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3139801932u,        0.9018090336u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1704086046u,        0.8924526367u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1704086046u,        0.8924526367u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2437707573u,        0.8729420069u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2437707573u,        0.8729420069u) withcolor MidnightBlue;
draw (        0.1704086046u,        0.8924526367u)--(        0.1345597623u,        0.9593669896u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0589860212u,        0.9522060002u)--(        0.1345597623u,        0.9593669896u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1704086046u,        0.8924526367u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1704086046u,        0.8924526367u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0589860212u,        0.9522060002u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0589860212u,        0.9522060002u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1345597623u,        0.9593669896u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1345597623u,        0.9593669896u) withcolor Green;
draw (        0.1704086046u,        0.8924526367u)--(        0.2052369108u,        0.8250014751u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1636398933u,        0.7615006227u)--(        0.2052369108u,        0.8250014751u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1704086046u,        0.8924526367u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1704086046u,        0.8924526367u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1636398933u,        0.7615006227u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1636398933u,        0.7615006227u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2052369108u,        0.8250014751u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2052369108u,        0.8250014751u) withcolor MidnightBlue;
draw (        0.1704086046u,        0.8924526367u)--(        0.1288115872u,        0.8289517843u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1636398933u,        0.7615006227u)--(        0.1288115872u,        0.8289517843u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1704086046u,        0.8924526367u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1704086046u,        0.8924526367u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1636398933u,        0.7615006227u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1636398933u,        0.7615006227u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1288115872u,        0.8289517843u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1288115872u,        0.8289517843u) withcolor MidnightBlue;
draw (        0.1704086046u,        0.8924526367u)--(        0.1288060490u,        0.8289554125u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0573993480u,        0.8547188200u)--(        0.1288060490u,        0.8289554125u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1704086046u,        0.8924526367u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1704086046u,        0.8924526367u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0573993480u,        0.8547188200u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0573993480u,        0.8547188200u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1288060490u,        0.8289554125u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1288060490u,        0.8289554125u) withcolor MidnightBlue;
draw (        0.0589860212u,        0.9522060002u)--(        0.0000000000u,        1.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0589860212u,        0.9522060002u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0589860212u,        0.9522060002u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        1.0000000000u) withcolor Yellow;
draw (        0.0589860212u,        0.9522060002u)--(        0.0000101578u,        0.9044093722u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0573993480u,        0.8547188200u)--(        0.0000101578u,        0.9044093722u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0589860212u,        0.9522060002u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0589860212u,        0.9522060002u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0573993480u,        0.8547188200u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0573993480u,        0.8547188200u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000101578u,        0.9044093722u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000101578u,        0.9044093722u) withcolor Green;
draw (        0.9538832839u,        0.2730338845u)--(        1.0000000000u,        0.3333333333u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9538832839u,        0.2730338845u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9538832839u,        0.2730338845u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        0.3333333333u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        0.3333333333u) withcolor Yellow;
draw (        0.9538832839u,        0.2730338845u)--(        0.9077712104u,        0.2127317553u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9426006519u,        0.1452811799u)--(        0.9077712104u,        0.2127317553u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9538832839u,        0.2730338845u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9538832839u,        0.2730338845u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9426006519u,        0.1452811799u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9426006519u,        0.1452811799u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9077712104u,        0.2127317553u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9077712104u,        0.2127317553u) withcolor MidnightBlue;
draw (        0.9538832839u,        0.2730338845u)--(        0.9077665680u,        0.2127353055u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8363601066u,        0.2384993773u)--(        0.9077665680u,        0.2127353055u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9538832839u,        0.2730338845u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9538832839u,        0.2730338845u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8363601066u,        0.2384993773u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8363601066u,        0.2384993773u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9077665680u,        0.2127353055u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9077665680u,        0.2127353055u) withcolor MidnightBlue;
draw (        0.9131685244u,        0.3150187097u)--(        0.8372587419u,        0.3144063108u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7932563251u,        0.3762646768u)--(        0.8372587419u,        0.3144063108u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9131685244u,        0.3150187097u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9131685244u,        0.3150187097u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7932563251u,        0.3762646768u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7932563251u,        0.3762646768u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8372587419u,        0.3144063108u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8372587419u,        0.3144063108u) withcolor MidnightBlue;
draw (        0.9131685244u,        0.3150187097u)--(        0.8372587419u,        0.3144063109u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8363601066u,        0.2384993773u)--(        0.8372587419u,        0.3144063109u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9131685244u,        0.3150187097u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9131685244u,        0.3150187097u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8363601066u,        0.2384993773u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8363601066u,        0.2384993773u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8372587419u,        0.3144063109u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8372587419u,        0.3144063109u) withcolor MidnightBlue;
draw (        0.7932563251u,        0.3762646768u)--(        0.7941549603u,        0.4521716104u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7183144197u,        0.4554704764u)--(        0.7941549603u,        0.4521716104u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7932563251u,        0.3762646768u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7932563251u,        0.3762646768u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7183144197u,        0.4554704764u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7183144197u,        0.4554704764u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7941549603u,        0.4521716104u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7941549603u,        0.4521716104u) withcolor Green;
draw (        0.7932563251u,        0.3762646768u)--(        0.8372587420u,        0.3144063109u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8363601066u,        0.2384993773u)--(        0.8372587420u,        0.3144063109u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7932563251u,        0.3762646768u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7932563251u,        0.3762646768u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8363601066u,        0.2384993773u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8363601066u,        0.2384993773u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8372587420u,        0.3144063109u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8372587420u,        0.3144063109u) withcolor MidnightBlue;
draw (        0.7932563251u,        0.3762646768u)--(        0.7923576898u,        0.3003577431u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8363601066u,        0.2384993773u)--(        0.7923576898u,        0.3003577431u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7932563251u,        0.3762646768u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7932563251u,        0.3762646768u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8363601066u,        0.2384993773u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8363601066u,        0.2384993773u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7923576898u,        0.3003577431u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7923576898u,        0.3003577431u) withcolor MidnightBlue;
draw (        0.7932563251u,        0.3762646768u)--(        0.7923463500u,        0.3003578782u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7186902638u,        0.3187277201u)--(        0.7923463500u,        0.3003578782u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7932563251u,        0.3762646768u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7932563251u,        0.3762646768u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7186902638u,        0.3187277201u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7186902638u,        0.3187277201u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7923463500u,        0.3003578782u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7923463500u,        0.3003578782u) withcolor MidnightBlue;
draw (        0.7183144197u,        0.4554704764u)--(        0.6938668391u,        0.5273383129u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6187378042u,        0.5164618138u)--(        0.6938668391u,        0.5273383129u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7183144197u,        0.4554704764u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7183144197u,        0.4554704764u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6187378042u,        0.5164618138u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6187378042u,        0.5164618138u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6938668391u,        0.5273383129u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6938668391u,        0.5273383129u) withcolor MidnightBlue;
draw (        0.7183144197u,        0.4554704764u)--(        0.6855177813u,        0.3870084386u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7186902638u,        0.3187277201u)--(        0.6855177813u,        0.3870084386u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7183144197u,        0.4554704764u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7183144197u,        0.4554704764u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7186902638u,        0.3187277201u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7186902638u,        0.3187277201u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6855177813u,        0.3870084386u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6855177813u,        0.3870084386u) withcolor MidnightBlue;
draw (        0.7183144197u,        0.4554704764u)--(        0.6855102589u,        0.3870120427u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6096263047u,        0.3849394578u)--(        0.6855102589u,        0.3870120427u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7183144197u,        0.4554704764u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7183144197u,        0.4554704764u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6096263047u,        0.3849394578u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6096263047u,        0.3849394578u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6855102589u,        0.3870120427u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6855102589u,        0.3870120427u) withcolor MidnightBlue;
draw (        0.6187378042u,        0.5164618138u)--(        0.5443008235u,        0.5015686603u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4978230044u,        0.5615893438u)--(        0.5443008235u,        0.5015686603u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6187378042u,        0.5164618138u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6187378042u,        0.5164618138u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4978230044u,        0.5615893438u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4978230044u,        0.5615893438u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5443008235u,        0.5015686603u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5443008235u,        0.5015686603u) withcolor MidnightBlue;
draw (        0.6187378042u,        0.5164618138u)--(        0.5722599851u,        0.5764824973u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4978230044u,        0.5615893438u)--(        0.5722599851u,        0.5764824973u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6187378042u,        0.5164618138u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6187378042u,        0.5164618138u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4978230044u,        0.5615893438u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4978230044u,        0.5615893438u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5722599851u,        0.5764824973u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5722599851u,        0.5764824973u) withcolor MidnightBlue;
draw (        0.6187378042u,        0.5164618138u)--(        0.5766239297u,        0.4533025572u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6096263047u,        0.3849394578u)--(        0.5766239297u,        0.4533025572u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6187378042u,        0.5164618138u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6187378042u,        0.5164618138u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6096263047u,        0.3849394578u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6096263047u,        0.3849394578u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5766239297u,        0.4533025572u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5766239297u,        0.4533025572u) withcolor MidnightBlue;
draw (        0.6187378042u,        0.5164618138u)--(        0.5766127326u,        0.4533100247u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5021769957u,        0.4384106561u)--(        0.5766127326u,        0.4533100247u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6187378042u,        0.5164618138u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6187378042u,        0.5164618138u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5021769957u,        0.4384106561u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5021769957u,        0.4384106561u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5766127326u,        0.4533100247u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5766127326u,        0.4533100247u) withcolor MidnightBlue;
draw (        0.6187378042u,        0.5164618138u)--(        0.5443020672u,        0.5015624451u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5021769957u,        0.4384106561u)--(        0.5443020672u,        0.5015624451u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6187378042u,        0.5164618138u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6187378042u,        0.5164618138u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5021769957u,        0.4384106561u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5021769957u,        0.4384106561u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5443020672u,        0.5015624451u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5443020672u,        0.5015624451u) withcolor MidnightBlue;
draw (        0.4978230044u,        0.5615893438u)--(        0.4233849185u,        0.5467017151u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3903736952u,        0.6150605422u)--(        0.4233849185u,        0.5467017151u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4978230044u,        0.5615893438u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4978230044u,        0.5615893438u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3903736952u,        0.6150605422u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3903736952u,        0.6150605422u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4233849185u,        0.5467017151u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4233849185u,        0.5467017151u) withcolor MidnightBlue;
draw (        0.4978230044u,        0.5615893438u)--(        0.5442970869u,        0.5015657670u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5021769957u,        0.4384106561u)--(        0.5442970869u,        0.5015657670u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4978230044u,        0.5615893438u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4978230044u,        0.5615893438u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5021769957u,        0.4384106561u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5021769957u,        0.4384106561u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5442970869u,        0.5015657670u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5442970869u,        0.5015657670u) withcolor MidnightBlue;
draw (        0.4978230044u,        0.5615893438u)--(        0.4557029132u,        0.4984342329u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5021769957u,        0.4384106561u)--(        0.4557029132u,        0.4984342329u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4978230044u,        0.5615893438u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4978230044u,        0.5615893438u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5021769957u,        0.4384106561u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5021769957u,        0.4384106561u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4557029132u,        0.4984342329u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4557029132u,        0.4984342329u) withcolor MidnightBlue;
draw (        0.4978230044u,        0.5615893438u)--(        0.4556979328u,        0.4984375548u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3812621958u,        0.4835381863u)--(        0.4556979328u,        0.4984375548u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4978230044u,        0.5615893438u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4978230044u,        0.5615893438u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3812621958u,        0.4835381863u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3812621958u,        0.4835381863u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4556979328u,        0.4984375548u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4556979328u,        0.4984375548u) withcolor MidnightBlue;
draw (        0.4978230044u,        0.5615893438u)--(        0.4233872674u,        0.5466899754u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3812621958u,        0.4835381863u)--(        0.4233872674u,        0.5466899754u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4978230044u,        0.5615893438u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4978230044u,        0.5615893438u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3812621958u,        0.4835381863u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3812621958u,        0.4835381863u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4233872674u,        0.5466899754u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4233872674u,        0.5466899754u) withcolor MidnightBlue;
draw (        0.3903736952u,        0.6150605422u)--(        0.3144895456u,        0.6129951216u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2813097361u,        0.6812722800u)--(        0.3144895456u,        0.6129951216u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3903736952u,        0.6150605422u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3903736952u,        0.6150605422u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2813097361u,        0.6812722800u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2813097361u,        0.6812722800u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3144895456u,        0.6129951216u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3144895456u,        0.6129951216u) withcolor MidnightBlue;
draw (        0.3903736952u,        0.6150605422u)--(        0.3571938857u,        0.6833377006u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2813097361u,        0.6812722800u)--(        0.3571938857u,        0.6833377006u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3903736952u,        0.6150605422u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3903736952u,        0.6150605422u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2813097361u,        0.6812722800u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2813097361u,        0.6812722800u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3571938857u,        0.6833377006u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3571938857u,        0.6833377006u) withcolor MidnightBlue;
draw (        0.3903736952u,        0.6150605422u)--(        0.4233760703u,        0.5466974429u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3812621958u,        0.4835381863u)--(        0.4233760703u,        0.5466974429u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3903736952u,        0.6150605422u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3903736952u,        0.6150605422u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3812621958u,        0.4835381863u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3812621958u,        0.4835381863u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4233760703u,        0.5466974429u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4233760703u,        0.5466974429u) withcolor MidnightBlue;
draw (        0.3903736952u,        0.6150605422u)--(        0.3144897409u,        0.6129879573u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2816855803u,        0.5445295235u)--(        0.3144897409u,        0.6129879573u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3903736952u,        0.6150605422u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3903736952u,        0.6150605422u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2816855803u,        0.5445295235u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2816855803u,        0.5445295235u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3144897409u,        0.6129879573u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3144897409u,        0.6129879573u) withcolor MidnightBlue;
draw (        0.2813097361u,        0.6812722800u)--(        0.2076561381u,        0.6996520955u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1636398933u,        0.7615006227u)--(        0.2076561381u,        0.6996520955u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2813097361u,        0.6812722800u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2813097361u,        0.6812722800u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1636398933u,        0.7615006227u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1636398933u,        0.7615006227u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2076561381u,        0.6996520955u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2076561381u,        0.6996520955u) withcolor MidnightBlue;
draw (        0.2813097361u,        0.6812722800u)--(        0.3144822185u,        0.6129915614u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2816855803u,        0.5445295235u)--(        0.3144822185u,        0.6129915614u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2813097361u,        0.6812722800u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2813097361u,        0.6812722800u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2816855803u,        0.5445295235u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2816855803u,        0.5445295235u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3144822185u,        0.6129915614u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3144822185u,        0.6129915614u) withcolor MidnightBlue;
draw (        0.2813097361u,        0.6812722800u)--(        0.2076536499u,        0.6996421217u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2067436749u,        0.6237353231u)--(        0.2076536499u,        0.6996421217u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2813097361u,        0.6812722800u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2813097361u,        0.6812722800u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2067436749u,        0.6237353231u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2067436749u,        0.6237353231u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2076536499u,        0.6996421217u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2076536499u,        0.6996421217u) withcolor MidnightBlue;
draw (        0.1636398933u,        0.7615006227u)--(        0.0922360700u,        0.7872720046u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0573993480u,        0.8547188200u)--(        0.0922360700u,        0.7872720046u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1636398933u,        0.7615006227u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1636398933u,        0.7615006227u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0573993480u,        0.8547188200u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0573993480u,        0.8547188200u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0922360700u,        0.7872720046u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0922360700u,        0.7872720046u) withcolor MidnightBlue;
draw (        0.1636398933u,        0.7615006227u)--(        0.1288031713u,        0.8289474381u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0573993480u,        0.8547188200u)--(        0.1288031713u,        0.8289474381u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1636398933u,        0.7615006227u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1636398933u,        0.7615006227u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0573993480u,        0.8547188200u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0573993480u,        0.8547188200u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1288031713u,        0.8289474381u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1288031713u,        0.8289474381u) withcolor MidnightBlue;
draw (        0.1636398933u,        0.7615006227u)--(        0.2076423101u,        0.6996422567u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2067436749u,        0.6237353231u)--(        0.2076423101u,        0.6996422567u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1636398933u,        0.7615006227u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1636398933u,        0.7615006227u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2067436749u,        0.6237353231u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2067436749u,        0.6237353231u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2076423101u,        0.6996422567u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2076423101u,        0.6996422567u) withcolor MidnightBlue;
draw (        0.1636398933u,        0.7615006227u)--(        0.0922334320u,        0.7872646945u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0461167160u,        0.7269661156u)--(        0.0922334320u,        0.7872646945u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1636398933u,        0.7615006227u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1636398933u,        0.7615006227u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0461167160u,        0.7269661156u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0461167160u,        0.7269661156u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0922334320u,        0.7872646945u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0922334320u,        0.7872646945u) withcolor MidnightBlue;
draw (        0.0573993480u,        0.8547188200u)--(        0.0922287897u,        0.7872682447u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0461167160u,        0.7269661156u)--(        0.0922287897u,        0.7872682447u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0573993480u,        0.8547188200u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0573993480u,        0.8547188200u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0461167160u,        0.7269661156u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0461167160u,        0.7269661156u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0922287897u,        0.7872682447u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0922287897u,        0.7872682447u) withcolor MidnightBlue;
draw (        0.9426006519u,        0.1452811799u)--(        0.8711968286u,        0.1710525618u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8363601066u,        0.2384993773u)--(        0.8711968286u,        0.1710525618u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9426006519u,        0.1452811799u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9426006519u,        0.1452811799u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8363601066u,        0.2384993773u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8363601066u,        0.2384993773u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8711968286u,        0.1710525618u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8711968286u,        0.1710525618u) withcolor MidnightBlue;
draw (        0.9426006519u,        0.1452811799u)--(        0.9077639300u,        0.2127279954u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8363601066u,        0.2384993773u)--(        0.9077639300u,        0.2127279954u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9426006519u,        0.1452811799u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9426006519u,        0.1452811799u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8363601066u,        0.2384993773u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8363601066u,        0.2384993773u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9077639300u,        0.2127279954u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9077639300u,        0.2127279954u) withcolor MidnightBlue;
draw (        0.9426006519u,        0.1452811799u)--(        0.9999898422u,        0.0955906279u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.9410139788u,        0.0477939998u)--(        0.9999898422u,        0.0955906279u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9426006519u,        0.1452811799u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9426006519u,        0.1452811799u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9410139788u,        0.0477939998u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9410139788u,        0.0477939998u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9999898422u,        0.0955906279u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9999898422u,        0.0955906279u) withcolor Green;
draw (        0.9426006519u,        0.1452811799u)--(        0.8711939509u,        0.1710445875u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8295913953u,        0.1075473633u)--(        0.8711939509u,        0.1710445875u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9426006519u,        0.1452811799u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9426006519u,        0.1452811799u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8295913953u,        0.1075473633u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8295913953u,        0.1075473633u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8711939509u,        0.1710445875u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8711939509u,        0.1710445875u) withcolor MidnightBlue;
draw (        0.8363601066u,        0.2384993773u)--(        0.7923438617u,        0.3003479044u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7186902638u,        0.3187277201u)--(        0.7923438617u,        0.3003479044u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8363601066u,        0.2384993773u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8363601066u,        0.2384993773u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7186902638u,        0.3187277201u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7186902638u,        0.3187277201u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7923438617u,        0.3003479044u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7923438617u,        0.3003479044u) withcolor MidnightBlue;
draw (        0.8363601066u,        0.2384993773u)--(        0.8711884128u,        0.1710482156u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8295913953u,        0.1075473633u)--(        0.8711884128u,        0.1710482156u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8363601066u,        0.2384993773u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8363601066u,        0.2384993773u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8295913953u,        0.1075473633u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8295913953u,        0.1075473633u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8711884128u,        0.1710482156u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8711884128u,        0.1710482156u) withcolor MidnightBlue;
draw (        0.8363601066u,        0.2384993773u)--(        0.7947630891u,        0.1749985249u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8295913953u,        0.1075473633u)--(        0.7947630891u,        0.1749985249u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8363601066u,        0.2384993773u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8363601066u,        0.2384993773u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8295913953u,        0.1075473633u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8295913953u,        0.1075473633u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7947630891u,        0.1749985249u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7947630891u,        0.1749985249u) withcolor MidnightBlue;
draw (        0.8363601066u,        0.2384993773u)--(        0.7947603149u,        0.1750003423u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7213985109u,        0.1945122832u)--(        0.7947603149u,        0.1750003423u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8363601066u,        0.2384993773u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8363601066u,        0.2384993773u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7213985109u,        0.1945122832u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7213985109u,        0.1945122832u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7947603149u,        0.1750003423u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7947603149u,        0.1750003423u) withcolor MidnightBlue;
draw (        0.7186902638u,        0.3187277201u)--(        0.6428061142u,        0.3166622995u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6096263047u,        0.3849394578u)--(        0.6428061142u,        0.3166622995u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7186902638u,        0.3187277201u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7186902638u,        0.3187277201u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6096263047u,        0.3849394578u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6096263047u,        0.3849394578u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6428061142u,        0.3166622995u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6428061142u,        0.3166622995u) withcolor MidnightBlue;
draw (        0.7186902638u,        0.3187277201u)--(        0.6855104542u,        0.3870048784u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6096263047u,        0.3849394578u)--(        0.6855104542u,        0.3870048784u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7186902638u,        0.3187277201u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7186902638u,        0.3187277201u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6096263047u,        0.3849394578u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6096263047u,        0.3849394578u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6855104542u,        0.3870048784u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6855104542u,        0.3870048784u) withcolor MidnightBlue;
draw (        0.7186902638u,        0.3187277201u)--(        0.6764260040u,        0.2556689977u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7213985109u,        0.1945122832u)--(        0.6764260040u,        0.2556689977u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7186902638u,        0.3187277201u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7186902638u,        0.3187277201u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7213985109u,        0.1945122832u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7213985109u,        0.1945122832u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6764260040u,        0.2556689977u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6764260040u,        0.2556689977u) withcolor MidnightBlue;
draw (        0.7186902638u,        0.3187277201u)--(        0.6764212531u,        0.2556721822u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6005371915u,        0.2536035324u)--(        0.6764212531u,        0.2556721822u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7186902638u,        0.3187277201u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7186902638u,        0.3187277201u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6005371915u,        0.2536035324u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6005371915u,        0.2536035324u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6764212531u,        0.2556721822u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6764212531u,        0.2556721822u) withcolor MidnightBlue;
draw (        0.7186902638u,        0.3187277201u)--(        0.6428062022u,        0.3166590703u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6005371915u,        0.2536035324u)--(        0.6428062022u,        0.3166590703u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7186902638u,        0.3187277201u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7186902638u,        0.3187277201u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6005371915u,        0.2536035324u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6005371915u,        0.2536035324u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6428062022u,        0.3166590703u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6428062022u,        0.3166590703u) withcolor MidnightBlue;
draw (        0.6096263047u,        0.3849394578u)--(        0.5766150815u,        0.4532982851u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.5021769957u,        0.4384106561u)--(        0.5766150815u,        0.4532982851u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6096263047u,        0.3849394578u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6096263047u,        0.3849394578u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5021769957u,        0.4384106561u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5021769957u,        0.4384106561u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5766150815u,        0.4532982851u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5766150815u,        0.4532982851u) withcolor MidnightBlue;
draw (        0.6096263047u,        0.3849394578u)--(        0.6428033724u,        0.3166609671u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6005371915u,        0.2536035324u)--(        0.6428033724u,        0.3166609671u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6096263047u,        0.3849394578u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6096263047u,        0.3849394578u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6005371915u,        0.2536035324u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6005371915u,        0.2536035324u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6428033724u,        0.3166609671u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6428033724u,        0.3166609671u) withcolor MidnightBlue;
draw (        0.6096263047u,        0.3849394578u)--(        0.5673601238u,        0.3218820231u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6005371915u,        0.2536035324u)--(        0.5673601238u,        0.3218820231u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6096263047u,        0.3849394578u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6096263047u,        0.3849394578u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6005371915u,        0.2536035324u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6005371915u,        0.2536035324u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5673601238u,        0.3218820231u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5673601238u,        0.3218820231u) withcolor MidnightBlue;
draw (        0.6096263047u,        0.3849394578u)--(        0.5673538222u,        0.3218862474u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4923399034u,        0.3102422747u)--(        0.5673538222u,        0.3218862474u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6096263047u,        0.3849394578u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6096263047u,        0.3849394578u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4923399034u,        0.3102422747u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4923399034u,        0.3102422747u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5673538222u,        0.3218862474u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5673538222u,        0.3218862474u) withcolor MidnightBlue;
draw (        0.5021769957u,        0.4384106561u)--(        0.4277400149u,        0.4235175028u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3812621958u,        0.4835381863u)--(        0.4277400149u,        0.4235175028u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5021769957u,        0.4384106561u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5021769957u,        0.4384106561u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3812621958u,        0.4835381863u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3812621958u,        0.4835381863u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4277400149u,        0.4235175028u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4277400149u,        0.4235175028u) withcolor MidnightBlue;
draw (        0.5021769957u,        0.4384106561u)--(        0.4556991766u,        0.4984313396u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3812621958u,        0.4835381863u)--(        0.4556991766u,        0.4984313396u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5021769957u,        0.4384106561u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5021769957u,        0.4384106561u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3812621958u,        0.4835381863u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3812621958u,        0.4835381863u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4556991766u,        0.4984313396u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4556991766u,        0.4984313396u) withcolor MidnightBlue;
draw (        0.5021769957u,        0.4384106561u)--(        0.4569826642u,        0.3774176850u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4923399034u,        0.3102422747u)--(        0.4569826642u,        0.3774176850u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5021769957u,        0.4384106561u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5021769957u,        0.4384106561u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4923399034u,        0.3102422747u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4923399034u,        0.3102422747u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4569826642u,        0.3774176850u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4569826642u,        0.3774176850u) withcolor MidnightBlue;
draw (        0.5021769957u,        0.4384106561u)--(        0.4569725299u,        0.3774251956u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3825364813u,        0.3625273842u)--(        0.4569725299u,        0.3774251956u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5021769957u,        0.4384106561u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5021769957u,        0.4384106561u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3825364813u,        0.3625273842u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3825364813u,        0.3625273842u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4569725299u,        0.3774251956u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4569725299u,        0.3774251956u) withcolor MidnightBlue;
draw (        0.5021769957u,        0.4384106561u)--(        0.4277409471u,        0.4235128447u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3825364813u,        0.3625273842u)--(        0.4277409471u,        0.4235128447u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5021769957u,        0.4384106561u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5021769957u,        0.4384106561u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3825364813u,        0.3625273842u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3825364813u,        0.3625273842u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4277409471u,        0.4235128447u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4277409471u,        0.4235128447u) withcolor MidnightBlue;
draw (        0.3812621958u,        0.4835381863u)--(        0.3061331610u,        0.4726616870u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2816855803u,        0.5445295235u)--(        0.3061331610u,        0.4726616870u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3812621958u,        0.4835381863u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3812621958u,        0.4835381863u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2816855803u,        0.5445295235u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2816855803u,        0.5445295235u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3061331610u,        0.4726616870u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3061331610u,        0.4726616870u) withcolor MidnightBlue;
draw (        0.3812621958u,        0.4835381863u)--(        0.4277373972u,        0.4235154758u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3825364813u,        0.3625273842u)--(        0.4277373972u,        0.4235154758u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3812621958u,        0.4835381863u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3812621958u,        0.4835381863u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3825364813u,        0.3625273842u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3825364813u,        0.3625273842u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4277373972u,        0.4235154758u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4277373972u,        0.4235154758u) withcolor MidnightBlue;
draw (        0.3812621958u,        0.4835381863u)--(        0.3360612800u,        0.4225500947u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3825364813u,        0.3625273842u)--(        0.3360612800u,        0.4225500947u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3812621958u,        0.4835381863u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3812621958u,        0.4835381863u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3825364813u,        0.3625273842u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3825364813u,        0.3625273842u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3360612800u,        0.4225500947u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3360612800u,        0.4225500947u) withcolor MidnightBlue;
draw (        0.3812621958u,        0.4835381863u)--(        0.3360528241u,        0.4225563626u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2609251568u,        0.4116704214u)--(        0.3360528241u,        0.4225563626u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3812621958u,        0.4835381863u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3812621958u,        0.4835381863u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2609251568u,        0.4116704214u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2609251568u,        0.4116704214u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3360528241u,        0.4225563626u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3360528241u,        0.4225563626u) withcolor MidnightBlue;
draw (        0.3812621958u,        0.4835381863u)--(        0.3061345285u,        0.4726522452u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2609251568u,        0.4116704214u)--(        0.3061345285u,        0.4726522452u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3812621958u,        0.4835381863u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3812621958u,        0.4835381863u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2609251568u,        0.4116704214u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2609251568u,        0.4116704214u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3061345285u,        0.4726522452u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3061345285u,        0.4726522452u) withcolor MidnightBlue;
draw (        0.2816855803u,        0.5445295235u)--(        0.2058450397u,        0.5478283895u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2067436749u,        0.6237353231u)--(        0.2058450397u,        0.5478283895u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2816855803u,        0.5445295235u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2816855803u,        0.5445295235u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2067436749u,        0.6237353231u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2067436749u,        0.6237353231u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2058450397u,        0.5478283895u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2058450397u,        0.5478283895u) withcolor Green;
draw (        0.2816855803u,        0.5445295235u)--(        0.3061253668u,        0.4726590362u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2609251568u,        0.4116704214u)--(        0.3061253668u,        0.4726590362u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2816855803u,        0.5445295235u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2816855803u,        0.5445295235u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2609251568u,        0.4116704214u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2609251568u,        0.4116704214u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3061253668u,        0.4726590362u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3061253668u,        0.4726590362u) withcolor MidnightBlue;
draw (        0.0461167160u,        0.7269661156u)--(        0.0000000000u,        0.6666666667u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0461167160u,        0.7269661156u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0461167160u,        0.7269661156u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.6666666667u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.6666666667u) withcolor Yellow;
draw (        0.9410139788u,        0.0477939998u)--(        1.0000000000u,        0.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9410139788u,        0.0477939998u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9410139788u,        0.0477939998u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.9410139788u,        0.0477939998u)--(        0.8654402377u,        0.0406330104u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.8295913953u,        0.1075473633u)--(        0.8654402377u,        0.0406330104u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.9410139788u,        0.0477939998u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9410139788u,        0.0477939998u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8295913953u,        0.1075473633u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8295913953u,        0.1075473633u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8654402377u,        0.0406330104u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8654402377u,        0.0406330104u) withcolor Green;
draw (        0.8295913953u,        0.1075473633u)--(        0.7562292427u,        0.1270579931u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6860198068u,        0.0981909664u)--(        0.7562292427u,        0.1270579931u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8295913953u,        0.1075473633u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8295913953u,        0.1075473633u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6860198068u,        0.0981909664u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6860198068u,        0.0981909664u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7562292427u,        0.1270579931u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7562292427u,        0.1270579931u) withcolor MidnightBlue;
draw (        0.8295913953u,        0.1075473633u)--(        0.7562305935u,        0.1270630716u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7213985109u,        0.1945122832u)--(        0.7562305935u,        0.1270630716u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8295913953u,        0.1075473633u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8295913953u,        0.1075473633u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7213985109u,        0.1945122832u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7213985109u,        0.1945122832u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7562305935u,        0.1270630716u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7562305935u,        0.1270630716u) withcolor MidnightBlue;
draw (        0.8295913953u,        0.1075473633u)--(        0.7947593128u,        0.1749965749u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7213985109u,        0.1945122832u)--(        0.7947593128u,        0.1749965749u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8295913953u,        0.1075473633u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8295913953u,        0.1075473633u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7213985109u,        0.1945122832u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7213985109u,        0.1945122832u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7947593128u,        0.1749965749u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7947593128u,        0.1749965749u) withcolor MidnightBlue;
draw (        0.8295913953u,        0.1075473633u)--(        0.7996039030u,        0.0378091354u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7237370989u,        0.0351827035u)--(        0.7996039030u,        0.0378091354u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8295913953u,        0.1075473633u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8295913953u,        0.1075473633u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7237370989u,        0.0351827035u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7237370989u,        0.0351827035u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7996039030u,        0.0378091354u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7996039030u,        0.0378091354u) withcolor Green;
draw (        0.6860198068u,        0.0981909664u)--(        0.6511905609u,        0.1656416429u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7213985109u,        0.1945122832u)--(        0.6511905609u,        0.1656416429u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6860198068u,        0.0981909664u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6860198068u,        0.0981909664u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7213985109u,        0.1945122832u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7213985109u,        0.1945122832u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6511905609u,        0.1656416429u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6511905609u,        0.1656416429u) withcolor Green;
draw (        0.6860198068u,        0.0981909664u)--(        0.7562277568u,        0.1270616068u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7213985109u,        0.1945122832u)--(        0.7562277568u,        0.1270616068u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6860198068u,        0.0981909664u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6860198068u,        0.0981909664u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7213985109u,        0.1945122832u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7213985109u,        0.1945122832u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7562277568u,        0.1270616068u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7562277568u,        0.1270616068u) withcolor MidnightBlue;
draw (        0.6860198068u,        0.0981909664u)--(        0.6478701258u,        0.0325611565u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.7237370989u,        0.0351827035u)--(        0.6478701258u,        0.0325611565u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6860198068u,        0.0981909664u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6860198068u,        0.0981909664u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7237370989u,        0.0351827035u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7237370989u,        0.0351827035u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6478701258u,        0.0325611565u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6478701258u,        0.0325611565u) withcolor Green;
draw (        0.7213985109u,        0.1945122832u)--(        0.6455142711u,        0.1924501807u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6005371915u,        0.2536035324u)--(        0.6455142711u,        0.1924501807u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7213985109u,        0.1945122832u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7213985109u,        0.1945122832u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6005371915u,        0.2536035324u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6005371915u,        0.2536035324u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6455142711u,        0.1924501807u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6455142711u,        0.1924501807u) withcolor Green;
draw (        0.7213985109u,        0.1945122832u)--(        0.6764214313u,        0.2556656349u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.6005371915u,        0.2536035324u)--(        0.6764214313u,        0.2556656349u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7213985109u,        0.1945122832u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7213985109u,        0.1945122832u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6005371915u,        0.2536035324u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6005371915u,        0.2536035324u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6764214313u,        0.2556656349u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6764214313u,        0.2556656349u) withcolor MidnightBlue;
draw (        0.6005371915u,        0.2536035324u)--(        0.5555555556u,        0.1924500897u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6005371915u,        0.2536035324u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6005371915u,        0.2536035324u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5555555556u,        0.1924500897u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5555555556u,        0.1924500897u) withcolor Yellow;
draw (        0.6005371915u,        0.2536035324u)--(        0.5673548736u,        0.3218794717u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.4923399034u,        0.3102422747u)--(        0.5673548736u,        0.3218794717u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6005371915u,        0.2536035324u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6005371915u,        0.2536035324u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4923399034u,        0.3102422747u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4923399034u,        0.3102422747u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5673548736u,        0.3218794717u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5673548736u,        0.3218794717u) withcolor MidnightBlue;
draw (        0.4923399034u,        0.3102422747u)--(        0.4569748549u,        0.3774135741u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3825364813u,        0.3625273842u)--(        0.4569748549u,        0.3774135741u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4923399034u,        0.3102422747u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4923399034u,        0.3102422747u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3825364813u,        0.3625273842u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3825364813u,        0.3625273842u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4569748549u,        0.3774135741u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4569748549u,        0.3774135741u) withcolor MidnightBlue;
draw (        0.4923399034u,        0.3102422747u)--(        0.4704714092u,        0.2375481272u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3956976459u,        0.2506460059u)--(        0.4704714092u,        0.2375481272u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4923399034u,        0.3102422747u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4923399034u,        0.3102422747u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3956976459u,        0.2506460059u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3956976459u,        0.2506460059u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4704714092u,        0.2375481272u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4704714092u,        0.2375481272u) withcolor Green;
draw (        0.3825364813u,        0.3625273842u)--(        0.3360544930u,        0.4225448390u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2609251568u,        0.4116704214u)--(        0.3360544930u,        0.4225448390u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3825364813u,        0.3625273842u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3825364813u,        0.3625273842u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2609251568u,        0.4116704214u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2609251568u,        0.4116704214u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3360544930u,        0.4225448390u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3360544930u,        0.4225448390u) withcolor MidnightBlue;
draw (        0.3825364813u,        0.3625273842u)--(        0.3385734777u,        0.3006410009u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3956976459u,        0.2506460059u)--(        0.3385734777u,        0.3006410009u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3825364813u,        0.3625273842u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3825364813u,        0.3625273842u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3956976459u,        0.2506460059u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3956976459u,        0.2506460059u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3385734777u,        0.3006410009u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3385734777u,        0.3006410009u) withcolor MidnightBlue;
draw (        0.3825364813u,        0.3625273842u)--(        0.3385608162u,        0.3006499974u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2626863237u,        0.2982559405u)--(        0.3385608162u,        0.3006499974u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3825364813u,        0.3625273842u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3825364813u,        0.3625273842u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2626863237u,        0.2982559405u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2626863237u,        0.2982559405u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3385608162u,        0.3006499974u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3385608162u,        0.3006499974u) withcolor MidnightBlue;
draw (        0.2609251568u,        0.4116704214u)--(        0.1924500897u,        0.4444444444u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2609251568u,        0.4116704214u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2609251568u,        0.4116704214u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1924500897u,        0.4444444444u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1924500897u,        0.4444444444u) withcolor Yellow;
draw (        0.2609251568u,        0.4116704214u)--(        0.2113521008u,        0.3541797071u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2626863237u,        0.2982559405u)--(        0.2113521008u,        0.3541797071u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2609251568u,        0.4116704214u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2609251568u,        0.4116704214u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2626863237u,        0.2982559405u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2626863237u,        0.2982559405u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2113521008u,        0.3541797071u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2113521008u,        0.3541797071u) withcolor MidnightBlue;
draw (        0.2609251568u,        0.4116704214u)--(        0.2113382303u,        0.3541916702u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1358560484u,        0.3461225723u)--(        0.2113382303u,        0.3541916702u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2609251568u,        0.4116704214u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2609251568u,        0.4116704214u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1358560484u,        0.3461225723u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1358560484u,        0.3461225723u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2113382303u,        0.3541916702u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2113382303u,        0.3541916702u) withcolor MidnightBlue;
draw (        0.3956976459u,        0.2506460059u)--(        0.4444444444u,        0.1924500897u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3956976459u,        0.2506460059u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3956976459u,        0.2506460059u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4444444444u,        0.1924500897u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4444444444u,        0.1924500897u) withcolor Yellow;
draw (        0.3956976459u,        0.2506460059u)--(        0.3198224420u,        0.2482746029u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2626863237u,        0.2982559405u)--(        0.3198224420u,        0.2482746029u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3956976459u,        0.2506460059u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3956976459u,        0.2506460059u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2626863237u,        0.2982559405u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2626863237u,        0.2982559405u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3198224420u,        0.2482746029u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3198224420u,        0.2482746029u) withcolor MidnightBlue;
draw (        0.3956976459u,        0.2506460059u)--(        0.3385615276u,        0.3006273435u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2626863237u,        0.2982559405u)--(        0.3385615276u,        0.3006273435u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3956976459u,        0.2506460059u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3956976459u,        0.2506460059u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2626863237u,        0.2982559405u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2626863237u,        0.2982559405u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3385615276u,        0.3006273435u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3385615276u,        0.3006273435u) withcolor MidnightBlue;
draw (        0.3956976459u,        0.2506460059u)--(        0.3198230051u,        0.2482566533u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2776389212u,        0.1851442680u)--(        0.3198230051u,        0.2482566533u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3956976459u,        0.2506460059u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3956976459u,        0.2506460059u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2776389212u,        0.1851442680u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2776389212u,        0.1851442680u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3198230051u,        0.2482566533u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3198230051u,        0.2482566533u) withcolor MidnightBlue;
draw (        0.2626863237u,        0.2982559405u)--(        0.3198127395u,        0.2482635136u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2776389212u,        0.1851442680u)--(        0.3198127395u,        0.2482635136u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2626863237u,        0.2982559405u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2626863237u,        0.2982559405u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2776389212u,        0.1851442680u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2776389212u,        0.1851442680u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3198127395u,        0.2482635136u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3198127395u,        0.2482635136u) withcolor MidnightBlue;
draw (        0.2626863237u,        0.2982559405u)--(        0.1872017407u,        0.2902093359u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1358560484u,        0.3461225723u)--(        0.1872017407u,        0.2902093359u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2626863237u,        0.2982559405u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2626863237u,        0.2982559405u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1358560484u,        0.3461225723u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1358560484u,        0.3461225723u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1872017407u,        0.2902093359u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1872017407u,        0.2902093359u) withcolor MidnightBlue;
draw (        0.2626863237u,        0.2982559405u)--(        0.2113406315u,        0.3541691769u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1358560484u,        0.3461225723u)--(        0.2113406315u,        0.3541691769u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2626863237u,        0.2982559405u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2626863237u,        0.2982559405u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1358560484u,        0.3461225723u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1358560484u,        0.3461225723u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2113406315u,        0.3541691769u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2113406315u,        0.3541691769u) withcolor MidnightBlue;
draw (        0.2626863237u,        0.2982559405u)--(        0.1872042189u,        0.2901861227u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1481389566u,        0.2250971421u)--(        0.1872042189u,        0.2901861227u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2626863237u,        0.2982559405u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2626863237u,        0.2982559405u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1481389566u,        0.2250971421u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1481389566u,        0.2250971421u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1872042189u,        0.2901861227u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1872042189u,        0.2901861227u) withcolor MidnightBlue;
draw (        0.2776389212u,        0.1851442680u)--(        0.3502784208u,        0.1630949321u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3129781928u,        0.0969786473u)--(        0.3502784208u,        0.1630949321u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2776389212u,        0.1851442680u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2776389212u,        0.1851442680u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3129781928u,        0.0969786473u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3129781928u,        0.0969786473u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3502784208u,        0.1630949321u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3502784208u,        0.1630949321u) withcolor Green;
draw (        0.2776389212u,        0.1851442680u)--(        0.2403386932u,        0.1190279833u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.3129781928u,        0.0969786473u)--(        0.2403386932u,        0.1190279833u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2776389212u,        0.1851442680u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2776389212u,        0.1851442680u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3129781928u,        0.0969786473u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3129781928u,        0.0969786473u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2403386932u,        0.1190279833u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2403386932u,        0.1190279833u) withcolor MidnightBlue;
draw (        0.2776389212u,        0.1851442680u)--(        0.2028005176u,        0.1724209248u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1481389566u,        0.2250971421u)--(        0.2028005176u,        0.1724209248u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2776389212u,        0.1851442680u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2776389212u,        0.1851442680u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1481389566u,        0.2250971421u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1481389566u,        0.2250971421u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2028005176u,        0.1724209248u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2028005176u,        0.1724209248u) withcolor MidnightBlue;
draw (        0.2776389212u,        0.1851442680u)--(        0.2403304843u,        0.1190326151u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1654930621u,        0.1063035006u)--(        0.2403304843u,        0.1190326151u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2776389212u,        0.1851442680u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2776389212u,        0.1851442680u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1654930621u,        0.1063035006u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1654930621u,        0.1063035006u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2403304843u,        0.1190326151u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2403304843u,        0.1190326151u) withcolor MidnightBlue;
draw (        0.2776389212u,        0.1851442680u)--(        0.2028014990u,        0.1724151535u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1654930621u,        0.1063035006u)--(        0.2028014990u,        0.1724151535u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2776389212u,        0.1851442680u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2776389212u,        0.1851442680u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1654930621u,        0.1063035006u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1654930621u,        0.1063035006u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2028014990u,        0.1724151535u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2028014990u,        0.1724151535u) withcolor MidnightBlue;
draw (        0.1358560484u,        0.3461225723u)--(        0.1871874290u,        0.2901961967u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1481389566u,        0.2250971421u)--(        0.1871874290u,        0.2901961967u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1358560484u,        0.3461225723u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1358560484u,        0.3461225723u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1481389566u,        0.2250971421u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1481389566u,        0.2250971421u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1871874290u,        0.2901961967u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1871874290u,        0.2901961967u) withcolor MidnightBlue;
draw (        0.3129781928u,        0.0969786473u)--(        0.3888888889u,        0.0962250449u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3129781928u,        0.0969786473u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3129781928u,        0.0969786473u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3888888889u,        0.0962250449u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3888888889u,        0.0962250449u) withcolor Yellow;
draw (        0.3129781928u,        0.0969786473u)--(        0.3515790915u,        0.0316132074u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.2757480427u,        0.0351234982u)--(        0.3515790915u,        0.0316132074u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3129781928u,        0.0969786473u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3129781928u,        0.0969786473u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2757480427u,        0.0351234982u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2757480427u,        0.0351234982u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3515790915u,        0.0316132074u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3515790915u,        0.0316132074u) withcolor Green;
draw (        0.3129781928u,        0.0969786473u)--(        0.2403339360u,        0.1190123050u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1654930621u,        0.1063035006u)--(        0.2403339360u,        0.1190123050u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3129781928u,        0.0969786473u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3129781928u,        0.0969786473u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1654930621u,        0.1063035006u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1654930621u,        0.1063035006u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2403339360u,        0.1190123050u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2403339360u,        0.1190123050u) withcolor MidnightBlue;
draw (        0.1481389566u,        0.2250971421u)--(        0.0725207558u,        0.2317723061u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0373382862u,        0.2990394151u)--(        0.0725207558u,        0.2317723061u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1481389566u,        0.2250971421u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1481389566u,        0.2250971421u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0373382862u,        0.2990394151u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0373382862u,        0.2990394151u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0725207558u,        0.2317723061u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0725207558u,        0.2317723061u) withcolor MidnightBlue;
draw (        0.1481389566u,        0.2250971421u)--(        0.2027972605u,        0.1724175452u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1654930621u,        0.1063035006u)--(        0.2027972605u,        0.1724175452u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1481389566u,        0.2250971421u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1481389566u,        0.2250971421u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1654930621u,        0.1063035006u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1654930621u,        0.1063035006u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2027972605u,        0.1724175452u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2027972605u,        0.1724175452u) withcolor MidnightBlue;
draw (        0.1481389566u,        0.2250971421u)--(        0.1108347582u,        0.1589830975u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1654930621u,        0.1063035006u)--(        0.1108347582u,        0.1589830975u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1481389566u,        0.2250971421u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1481389566u,        0.2250971421u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1654930621u,        0.1063035006u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1654930621u,        0.1063035006u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1108347582u,        0.1589830975u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1108347582u,        0.1589830975u) withcolor MidnightBlue;
draw (        0.1481389566u,        0.2250971421u)--(        0.1108207364u,        0.1589910110u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0352020467u,        0.1656606351u)--(        0.1108207364u,        0.1589910110u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1481389566u,        0.2250971421u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1481389566u,        0.2250971421u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0352020467u,        0.1656606351u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0352020467u,        0.1656606351u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1108207364u,        0.1589910110u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1108207364u,        0.1589910110u) withcolor MidnightBlue;
draw (        0.1481389566u,        0.2250971421u)--(        0.0725202669u,        0.2317667661u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0352020467u,        0.1656606351u)--(        0.0725202669u,        0.2317667661u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1481389566u,        0.2250971421u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1481389566u,        0.2250971421u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0352020467u,        0.1656606351u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0352020467u,        0.1656606351u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0725202669u,        0.2317667661u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0725202669u,        0.2317667661u) withcolor MidnightBlue;
draw (        0.0373382862u,        0.2990394151u)--(        0.0725154075u,        0.2317695091u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0352020467u,        0.1656606351u)--(        0.0725154075u,        0.2317695091u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0373382862u,        0.2990394151u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0373382862u,        0.2990394151u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0352020467u,        0.1656606351u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0352020467u,        0.1656606351u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0725154075u,        0.2317695091u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0725154075u,        0.2317695091u) withcolor MidnightBlue;
draw (        0.2757480427u,        0.0351234982u)--(        0.1999175295u,        0.0386453439u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.1654930621u,        0.1063035006u)--(        0.1999175295u,        0.0386453439u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2757480427u,        0.0351234982u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2757480427u,        0.0351234982u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1654930621u,        0.1063035006u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1654930621u,        0.1063035006u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1999175295u,        0.0386453439u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1999175295u,        0.0386453439u) withcolor Green;
draw (        0.1654930621u,        0.1063035006u)--(        0.1108185507u,        0.1589662760u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw (        0.0352020467u,        0.1656606351u)--(        0.1108185507u,        0.1589662760u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1654930621u,        0.1063035006u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1654930621u,        0.1063035006u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0352020467u,        0.1656606351u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0352020467u,        0.1656606351u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1108185507u,        0.1589662760u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1108185507u,        0.1589662760u) withcolor MidnightBlue;
draw (        0.0578382294u,        0.0491897262u)--(        0.0000000000u,        0.0000000000u) withpen pencircle scaled 1.5 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.0578382294u,        0.0491897262u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0578382294u,        0.0491897262u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor Yellow;
% r =         0.0759122527
endfig;
end;
