prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(9);
u := 20cm;
path pols[];
pair p[];
R := .2306250950u;
D := 2*R;

pols[1] = 
(        0.0000000000u,        0.0000000000u)--
(        1.0000000000u,        0.0000000000u)--
(        1.0000000000u,        1.0000000000u)--
(        0.0000000000u,        1.0000000000u)--
(        0.0000000000u,        0.0000000000u)--
cycle;
fill pols[1] withcolor 0.9Dandelion;
draw (        0.0000000000u,        0.0000000000u)--(        1.0000000000u,        0.0000000000u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        0.0000000000u)--(        1.0000000000u,        1.0000000000u) withpen pencircle scaled 1.0;
draw (        1.0000000000u,        1.0000000000u)--(        0.0000000000u,        1.0000000000u) withpen pencircle scaled 1.0;
draw (        0.0000000000u,        1.0000000000u)--(        0.0000000000u,        0.0000000000u) withpen pencircle scaled 1.0;

x1 = .8158038377u;
y1 = .8612206464u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .8565318223u;
y2 = .5418384185u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .4999999999u;
y3 = .8105800123u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .1841961624u;
y4 = .8612206465u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .8565767900u;
y5 = .1806168144u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .5000000000u;
y6 = .4492984073u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .1434681778u;
y7 = .5418384185u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .5000000001u;
y8 = .0880899664u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .1434232100u;
y9 = .1806168144u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.8158038377u,        0.8612206464u)--(        1.0000000000u,        1.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8158038377u,        0.8612206464u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8158038377u,        0.8612206464u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        1.0000000000u) withcolor Yellow;
draw (        0.8158038377u,        0.8612206464u)--(        0.9999835788u,        0.7224195029u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.8565318223u,        0.5418384185u)--(        0.9999835788u,        0.7224195029u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8158038377u,        0.8612206464u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8158038377u,        0.8612206464u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8565318223u,        0.5418384185u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8565318223u,        0.5418384185u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9999835788u,        0.7224195029u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9999835788u,        0.7224195029u) withcolor Green;
draw (        0.8158038377u,        0.8612206464u)--(        0.6315912581u,        0.9999782045u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.4999999999u,        0.8105800123u)--(        0.6315912581u,        0.9999782045u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8158038377u,        0.8612206464u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8158038377u,        0.8612206464u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4999999999u,        0.8105800123u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999999u,        0.8105800123u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6315912581u,        0.9999782045u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6315912581u,        0.9999782045u) withcolor Green;
draw (        0.8565318223u,        0.5418384185u)--(        0.6434079249u,        0.6299641173u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.4999999999u,        0.8105800123u)--(        0.6434079249u,        0.6299641173u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8565318223u,        0.5418384185u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8565318223u,        0.5418384185u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4999999999u,        0.8105800123u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999999u,        0.8105800123u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6434079249u,        0.6299641173u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6434079249u,        0.6299641173u) withcolor MidnightBlue;
draw (        0.8565318223u,        0.5418384185u)--(        0.9999686423u,        0.3612454698u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.8565767900u,        0.1806168144u)--(        0.9999686423u,        0.3612454698u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8565318223u,        0.5418384185u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8565318223u,        0.5418384185u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8565767900u,        0.1806168144u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8565767900u,        0.1806168144u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.9999686423u,        0.3612454698u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.9999686423u,        0.3612454698u) withcolor Green;
draw (        0.8565318223u,        0.5418384185u)--(        0.7131399700u,        0.3612097631u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.8565767900u,        0.1806168144u)--(        0.7131399700u,        0.3612097631u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8565318223u,        0.5418384185u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8565318223u,        0.5418384185u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8565767900u,        0.1806168144u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8565767900u,        0.1806168144u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7131399700u,        0.3612097631u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7131399700u,        0.3612097631u) withcolor MidnightBlue;
draw (        0.8565318223u,        0.5418384185u)--(        0.7131394020u,        0.3612102141u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5000000000u,        0.4492984073u)--(        0.7131394020u,        0.3612102141u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8565318223u,        0.5418384185u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8565318223u,        0.5418384185u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7131394020u,        0.3612102141u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7131394020u,        0.3612102141u) withcolor MidnightBlue;
draw (        0.8565318223u,        0.5418384185u)--(        0.6433924204u,        0.6299266117u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5000000000u,        0.4492984073u)--(        0.6433924204u,        0.6299266117u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8565318223u,        0.5418384185u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8565318223u,        0.5418384185u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6433924204u,        0.6299266117u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6433924204u,        0.6299266117u) withcolor MidnightBlue;
draw (        0.4999999999u,        0.8105800123u)--(        0.3684087420u,        0.9999782047u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.1841961624u,        0.8612206465u)--(        0.3684087420u,        0.9999782047u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4999999999u,        0.8105800123u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999999u,        0.8105800123u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1841961624u,        0.8612206465u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1841961624u,        0.8612206465u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3684087420u,        0.9999782047u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3684087420u,        0.9999782047u) withcolor Green;
draw (        0.4999999999u,        0.8105800123u)--(        0.6433765494u,        0.6299392098u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5000000000u,        0.4492984073u)--(        0.6433765494u,        0.6299392098u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4999999999u,        0.8105800123u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999999u,        0.8105800123u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6433765494u,        0.6299392098u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6433765494u,        0.6299392098u) withcolor MidnightBlue;
draw (        0.4999999999u,        0.8105800123u)--(        0.3566234505u,        0.6299392098u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5000000000u,        0.4492984073u)--(        0.3566234505u,        0.6299392098u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4999999999u,        0.8105800123u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999999u,        0.8105800123u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3566234505u,        0.6299392098u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3566234505u,        0.6299392098u) withcolor MidnightBlue;
draw (        0.4999999999u,        0.8105800123u)--(        0.3565920754u,        0.6299641170u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.1434681778u,        0.5418384185u)--(        0.3565920754u,        0.6299641170u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4999999999u,        0.8105800123u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999999u,        0.8105800123u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1434681778u,        0.5418384185u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1434681778u,        0.5418384185u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3565920754u,        0.6299641170u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3565920754u,        0.6299641170u) withcolor MidnightBlue;
draw (        0.1841961624u,        0.8612206465u)--(        0.0000000000u,        1.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1841961624u,        0.8612206465u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1841961624u,        0.8612206465u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        1.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        1.0000000000u) withcolor Yellow;
draw (        0.1841961624u,        0.8612206465u)--(        0.0000164214u,        0.7224195029u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.1434681778u,        0.5418384185u)--(        0.0000164214u,        0.7224195029u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1841961624u,        0.8612206465u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1841961624u,        0.8612206465u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1434681778u,        0.5418384185u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1434681778u,        0.5418384185u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000164214u,        0.7224195029u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000164214u,        0.7224195029u) withcolor Green;
draw (        0.8565767900u,        0.1806168144u)--(        1.0000000000u,        0.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8565767900u,        0.1806168144u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8565767900u,        0.1806168144u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.8565767900u,        0.1806168144u)--(        0.6434379431u,        0.2687063506u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5000000000u,        0.4492984073u)--(        0.6434379431u,        0.2687063506u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8565767900u,        0.1806168144u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8565767900u,        0.1806168144u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6434379431u,        0.2687063506u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6434379431u,        0.2687063506u) withcolor MidnightBlue;
draw (        0.8565767900u,        0.1806168144u)--(        0.7131388469u,        0.3612088710u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5000000000u,        0.4492984073u)--(        0.7131388469u,        0.3612088710u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8565767900u,        0.1806168144u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8565767900u,        0.1806168144u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7131388469u,        0.3612088710u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7131388469u,        0.3612088710u) withcolor MidnightBlue;
draw (        0.8565767900u,        0.1806168144u)--(        0.7131464160u,        0.0000187461u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5000000001u,        0.0880899664u)--(        0.7131464160u,        0.0000187461u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8565767900u,        0.1806168144u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8565767900u,        0.1806168144u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.0880899664u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.0880899664u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7131464160u,        0.0000187461u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7131464160u,        0.0000187461u) withcolor Green;
draw (        0.8565767900u,        0.1806168144u)--(        0.6434303741u,        0.2686880347u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5000000001u,        0.0880899664u)--(        0.6434303741u,        0.2686880347u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8565767900u,        0.1806168144u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8565767900u,        0.1806168144u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.0880899664u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.0880899664u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6434303741u,        0.2686880347u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6434303741u,        0.2686880347u) withcolor MidnightBlue;
draw (        0.5000000000u,        0.4492984073u)--(        0.2868605981u,        0.3612102140u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.1434681778u,        0.5418384185u)--(        0.2868605981u,        0.3612102140u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1434681778u,        0.5418384185u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1434681778u,        0.5418384185u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2868605981u,        0.3612102140u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2868605981u,        0.3612102140u) withcolor MidnightBlue;
draw (        0.5000000000u,        0.4492984073u)--(        0.3566075797u,        0.6299266118u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.1434681778u,        0.5418384185u)--(        0.3566075797u,        0.6299266118u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1434681778u,        0.5418384185u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1434681778u,        0.5418384185u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3566075797u,        0.6299266118u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3566075797u,        0.6299266118u) withcolor MidnightBlue;
draw (        0.5000000000u,        0.4492984073u)--(        0.6434226273u,        0.2686941869u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5000000001u,        0.0880899664u)--(        0.6434226273u,        0.2686941869u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.0880899664u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.0880899664u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6434226273u,        0.2686941869u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6434226273u,        0.2686941869u) withcolor MidnightBlue;
draw (        0.5000000000u,        0.4492984073u)--(        0.3565773728u,        0.2686941868u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5000000001u,        0.0880899664u)--(        0.3565773728u,        0.2686941868u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.0880899664u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.0880899664u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3565773728u,        0.2686941868u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3565773728u,        0.2686941868u) withcolor MidnightBlue;
draw (        0.5000000000u,        0.4492984073u)--(        0.3565620569u,        0.2687063506u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.1434232100u,        0.1806168144u)--(        0.3565620569u,        0.2687063506u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1434232100u,        0.1806168144u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1434232100u,        0.1806168144u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3565620569u,        0.2687063506u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3565620569u,        0.2687063506u) withcolor MidnightBlue;
draw (        0.5000000000u,        0.4492984073u)--(        0.2868611531u,        0.3612088710u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.1434232100u,        0.1806168144u)--(        0.2868611531u,        0.3612088710u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.4492984073u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1434232100u,        0.1806168144u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1434232100u,        0.1806168144u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2868611531u,        0.3612088710u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2868611531u,        0.3612088710u) withcolor MidnightBlue;
draw (        0.1434681778u,        0.5418384185u)--(        0.2868600300u,        0.3612097631u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.1434232100u,        0.1806168144u)--(        0.2868600300u,        0.3612097631u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1434681778u,        0.5418384185u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1434681778u,        0.5418384185u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1434232100u,        0.1806168144u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1434232100u,        0.1806168144u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2868600300u,        0.3612097631u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2868600300u,        0.3612097631u) withcolor MidnightBlue;
draw (        0.1434681778u,        0.5418384185u)--(        0.0000313577u,        0.3612454698u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.1434232100u,        0.1806168144u)--(        0.0000313577u,        0.3612454698u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1434681778u,        0.5418384185u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1434681778u,        0.5418384185u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1434232100u,        0.1806168144u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1434232100u,        0.1806168144u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000313577u,        0.3612454698u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000313577u,        0.3612454698u) withcolor Green;
draw (        0.5000000001u,        0.0880899664u)--(        0.2868535841u,        0.0000187462u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.1434232100u,        0.1806168144u)--(        0.2868535841u,        0.0000187462u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.0880899664u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.0880899664u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1434232100u,        0.1806168144u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1434232100u,        0.1806168144u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2868535841u,        0.0000187462u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2868535841u,        0.0000187462u) withcolor Green;
draw (        0.5000000001u,        0.0880899664u)--(        0.3565696259u,        0.2686880346u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.1434232100u,        0.1806168144u)--(        0.3565696259u,        0.2686880346u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.0880899664u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.0880899664u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1434232100u,        0.1806168144u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1434232100u,        0.1806168144u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3565696259u,        0.2686880346u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3565696259u,        0.2686880346u) withcolor MidnightBlue;
draw (        0.1434232100u,        0.1806168144u)--(        0.0000000000u,        0.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1434232100u,        0.1806168144u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1434232100u,        0.1806168144u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor Yellow;
% r =         0.2306250950
endfig;
end;
