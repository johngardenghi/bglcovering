prologues := 3;
outputformat   := "svg";
outputtemplate := "%j.%o";
input mpcolornames;
beginfig(10);
u := 20cm;
path pols[];
pair p[];
R := .1443146290u;
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

x1 = .8749602999u;
y1 = .0721917670u;
draw fullcircle scaled D shifted z1 withcolor RoyalBlue withpen pencircle scaled 1.0;

x2 = .7499734972u;
y2 = .2886751973u;
draw fullcircle scaled D shifted z2 withcolor RoyalBlue withpen pencircle scaled 1.0;

x3 = .6249867124u;
y3 = .5051585962u;
draw fullcircle scaled D shifted z3 withcolor RoyalBlue withpen pencircle scaled 1.0;

x4 = .4999999456u;
y4 = .7216419635u;
draw fullcircle scaled D shifted z4 withcolor RoyalBlue withpen pencircle scaled 1.0;

x5 = .6249868207u;
y5 = .0721917359u;
draw fullcircle scaled D shifted z5 withcolor RoyalBlue withpen pencircle scaled 1.0;

x6 = .5000000001u;
y6 = .2886751348u;
draw fullcircle scaled D shifted z6 withcolor RoyalBlue withpen pencircle scaled 1.0;

x7 = .3750131973u;
y7 = .5051585020u;
draw fullcircle scaled D shifted z7 withcolor RoyalBlue withpen pencircle scaled 1.0;

x8 = .3750133057u;
y8 = .0721917045u;
draw fullcircle scaled D shifted z8 withcolor RoyalBlue withpen pencircle scaled 1.0;

x9 = .2500264668u;
y9 = .2886750720u;
draw fullcircle scaled D shifted z9 withcolor RoyalBlue withpen pencircle scaled 1.0;

x10 = .1250397542u;
y10 = .0721916732u;
draw fullcircle scaled D shifted z10 withcolor RoyalBlue withpen pencircle scaled 1.0;
draw (        0.8749602999u,        0.0721917670u)--(        1.0000000000u,        0.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8749602999u,        0.0721917670u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8749602999u,        0.0721917670u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        1.0000000000u,        0.0000000000u) withcolor Yellow;
draw (        0.8749602999u,        0.0721917670u)--(        0.7499867981u,        0.1443605689u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.7499734972u,        0.2886751973u)--(        0.7499867981u,        0.1443605689u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8749602999u,        0.0721917670u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8749602999u,        0.0721917670u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7499734972u,        0.2886751973u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7499734972u,        0.2886751973u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7499867981u,        0.1443605689u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7499867981u,        0.1443605689u) withcolor MidnightBlue;
draw (        0.8749602999u,        0.0721917670u)--(        0.8749469990u,        0.2165063953u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.7499734972u,        0.2886751973u)--(        0.8749469990u,        0.2165063953u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8749602999u,        0.0721917670u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8749602999u,        0.0721917670u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7499734972u,        0.2886751973u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7499734972u,        0.2886751973u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.8749469990u,        0.2165063953u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8749469990u,        0.2165063953u) withcolor Green;
draw (        0.8749602999u,        0.0721917670u)--(        0.7499735693u,        0.0000458780u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.6249868207u,        0.0721917359u)--(        0.7499735693u,        0.0000458780u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8749602999u,        0.0721917670u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8749602999u,        0.0721917670u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6249868207u,        0.0721917359u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6249868207u,        0.0721917359u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7499735693u,        0.0000458780u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7499735693u,        0.0000458780u) withcolor Green;
draw (        0.8749602999u,        0.0721917670u)--(        0.7499735513u,        0.1443376249u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.6249868207u,        0.0721917359u)--(        0.7499735513u,        0.1443376249u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.8749602999u,        0.0721917670u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.8749602999u,        0.0721917670u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6249868207u,        0.0721917359u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6249868207u,        0.0721917359u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7499735513u,        0.1443376249u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7499735513u,        0.1443376249u) withcolor MidnightBlue;
draw (        0.7499734972u,        0.2886751973u)--(        0.6249999772u,        0.3608439678u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.6249867124u,        0.5051585962u)--(        0.6249999772u,        0.3608439678u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7499734972u,        0.2886751973u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7499734972u,        0.2886751973u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6249867124u,        0.5051585962u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6249867124u,        0.5051585962u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6249999772u,        0.3608439678u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6249999772u,        0.3608439678u) withcolor MidnightBlue;
draw (        0.7499734972u,        0.2886751973u)--(        0.7499602323u,        0.4329898256u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.6249867124u,        0.5051585962u)--(        0.7499602323u,        0.4329898256u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7499734972u,        0.2886751973u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7499734972u,        0.2886751973u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6249867124u,        0.5051585962u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6249867124u,        0.5051585962u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7499602323u,        0.4329898256u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7499602323u,        0.4329898256u) withcolor Green;
draw (        0.7499734972u,        0.2886751973u)--(        0.7499603046u,        0.1443605689u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.6249868207u,        0.0721917359u)--(        0.7499603046u,        0.1443605689u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7499734972u,        0.2886751973u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7499734972u,        0.2886751973u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6249868207u,        0.0721917359u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6249868207u,        0.0721917359u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.7499603046u,        0.1443605689u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7499603046u,        0.1443605689u) withcolor MidnightBlue;
draw (        0.7499734972u,        0.2886751973u)--(        0.6250000133u,        0.2165063643u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.6249868207u,        0.0721917359u)--(        0.6250000133u,        0.2165063643u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7499734972u,        0.2886751973u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7499734972u,        0.2886751973u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6249868207u,        0.0721917359u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6249868207u,        0.0721917359u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6250000133u,        0.2165063643u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6250000133u,        0.2165063643u) withcolor MidnightBlue;
draw (        0.7499734972u,        0.2886751973u)--(        0.6249867667u,        0.2165293080u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5000000001u,        0.2886751348u)--(        0.6249867667u,        0.2165293080u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7499734972u,        0.2886751973u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7499734972u,        0.2886751973u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6249867667u,        0.2165293080u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6249867667u,        0.2165293080u) withcolor MidnightBlue;
draw (        0.7499734972u,        0.2886751973u)--(        0.6249867306u,        0.3608210241u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5000000001u,        0.2886751348u)--(        0.6249867306u,        0.3608210241u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.7499734972u,        0.2886751973u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.7499734972u,        0.2886751973u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6249867306u,        0.3608210241u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6249867306u,        0.3608210241u) withcolor MidnightBlue;
draw (        0.6249867124u,        0.5051585962u)--(        0.5000131742u,        0.5773273352u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.4999999456u,        0.7216419635u)--(        0.5000131742u,        0.5773273352u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6249867124u,        0.5051585962u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6249867124u,        0.5051585962u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4999999456u,        0.7216419635u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999456u,        0.7216419635u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000131742u,        0.5773273352u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000131742u,        0.5773273352u) withcolor MidnightBlue;
draw (        0.6249867124u,        0.5051585962u)--(        0.6249734839u,        0.6494732246u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.4999999456u,        0.7216419635u)--(        0.6249734839u,        0.6494732246u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6249867124u,        0.5051585962u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6249867124u,        0.5051585962u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4999999456u,        0.7216419635u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999456u,        0.7216419635u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6249734839u,        0.6494732246u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6249734839u,        0.6494732246u) withcolor Green;
draw (        0.6249867124u,        0.5051585962u)--(        0.6249734840u,        0.3608439678u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5000000001u,        0.2886751348u)--(        0.6249734840u,        0.3608439678u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6249867124u,        0.5051585962u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6249867124u,        0.5051585962u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6249734840u,        0.3608439678u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6249734840u,        0.3608439678u) withcolor MidnightBlue;
draw (        0.6249867124u,        0.5051585962u)--(        0.5000132285u,        0.4329897632u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5000000001u,        0.2886751348u)--(        0.5000132285u,        0.4329897632u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6249867124u,        0.5051585962u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6249867124u,        0.5051585962u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000132285u,        0.4329897632u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000132285u,        0.4329897632u) withcolor MidnightBlue;
draw (        0.6249867124u,        0.5051585962u)--(        0.4999999820u,        0.4330127067u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.3750131973u,        0.5051585020u)--(        0.4999999820u,        0.4330127067u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6249867124u,        0.5051585962u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6249867124u,        0.5051585962u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3750131973u,        0.5051585020u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3750131973u,        0.5051585020u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4999999820u,        0.4330127067u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999820u,        0.4330127067u) withcolor MidnightBlue;
draw (        0.6249867124u,        0.5051585962u)--(        0.4999999276u,        0.5773043915u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.3750131973u,        0.5051585020u)--(        0.4999999276u,        0.5773043915u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6249867124u,        0.5051585962u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6249867124u,        0.5051585962u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3750131973u,        0.5051585020u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3750131973u,        0.5051585020u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4999999276u,        0.5773043915u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999276u,        0.5773043915u) withcolor MidnightBlue;
draw (        0.4999999456u,        0.7216419635u)--(        0.5000000000u,        0.8660254038u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4999999456u,        0.7216419635u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999456u,        0.7216419635u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000000u,        0.8660254038u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000000u,        0.8660254038u) withcolor Yellow;
draw (        0.4999999456u,        0.7216419635u)--(        0.4999866811u,        0.5773273352u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.3750131973u,        0.5051585020u)--(        0.4999866811u,        0.5773273352u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4999999456u,        0.7216419635u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999456u,        0.7216419635u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3750131973u,        0.5051585020u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3750131973u,        0.5051585020u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4999866811u,        0.5773273352u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999866811u,        0.5773273352u) withcolor MidnightBlue;
draw (        0.4999999456u,        0.7216419635u)--(        0.3750264618u,        0.6494731304u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.3750131973u,        0.5051585020u)--(        0.3750264618u,        0.6494731304u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.4999999456u,        0.7216419635u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999999456u,        0.7216419635u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3750131973u,        0.5051585020u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3750131973u,        0.5051585020u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3750264618u,        0.6494731304u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3750264618u,        0.6494731304u) withcolor Green;
draw (        0.6249868207u,        0.0721917359u)--(        0.5000133007u,        0.1443605064u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5000000001u,        0.2886751348u)--(        0.5000133007u,        0.1443605064u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6249868207u,        0.0721917359u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6249868207u,        0.0721917359u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000133007u,        0.1443605064u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000133007u,        0.1443605064u) withcolor MidnightBlue;
draw (        0.6249868207u,        0.0721917359u)--(        0.6249735201u,        0.2165063643u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.5000000001u,        0.2886751348u)--(        0.6249735201u,        0.2165063643u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6249868207u,        0.0721917359u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6249868207u,        0.0721917359u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.6249735201u,        0.2165063643u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6249735201u,        0.2165063643u) withcolor MidnightBlue;
draw (        0.6249868207u,        0.0721917359u)--(        0.5000000723u,        0.0000458777u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.3750133057u,        0.0721917045u)--(        0.5000000723u,        0.0000458777u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6249868207u,        0.0721917359u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6249868207u,        0.0721917359u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3750133057u,        0.0721917045u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3750133057u,        0.0721917045u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000723u,        0.0000458777u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000723u,        0.0000458777u) withcolor Green;
draw (        0.6249868207u,        0.0721917359u)--(        0.5000000542u,        0.1443375627u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.3750133057u,        0.0721917045u)--(        0.5000000542u,        0.1443375627u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.6249868207u,        0.0721917359u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.6249868207u,        0.0721917359u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3750133057u,        0.0721917045u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3750133057u,        0.0721917045u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.5000000542u,        0.1443375627u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000542u,        0.1443375627u) withcolor MidnightBlue;
draw (        0.5000000001u,        0.2886751348u)--(        0.3750264618u,        0.3608438736u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.3750131973u,        0.5051585020u)--(        0.3750264618u,        0.3608438736u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3750131973u,        0.5051585020u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3750131973u,        0.5051585020u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3750264618u,        0.3608438736u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3750264618u,        0.3608438736u) withcolor MidnightBlue;
draw (        0.5000000001u,        0.2886751348u)--(        0.4999867355u,        0.4329897632u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.3750131973u,        0.5051585020u)--(        0.4999867355u,        0.4329897632u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3750131973u,        0.5051585020u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3750131973u,        0.5051585020u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4999867355u,        0.4329897632u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999867355u,        0.4329897632u) withcolor MidnightBlue;
draw (        0.5000000001u,        0.2886751348u)--(        0.4999868076u,        0.1443605064u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.3750133057u,        0.0721917045u)--(        0.4999868076u,        0.1443605064u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3750133057u,        0.0721917045u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3750133057u,        0.0721917045u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.4999868076u,        0.1443605064u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.4999868076u,        0.1443605064u) withcolor MidnightBlue;
draw (        0.5000000001u,        0.2886751348u)--(        0.3750264983u,        0.2165063329u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.3750133057u,        0.0721917045u)--(        0.3750264983u,        0.2165063329u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3750133057u,        0.0721917045u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3750133057u,        0.0721917045u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3750264983u,        0.2165063329u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3750264983u,        0.2165063329u) withcolor MidnightBlue;
draw (        0.5000000001u,        0.2886751348u)--(        0.3750132516u,        0.2165292768u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2500264668u,        0.2886750720u)--(        0.3750132516u,        0.2165292768u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2500264668u,        0.2886750720u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500264668u,        0.2886750720u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3750132516u,        0.2165292768u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3750132516u,        0.2165292768u) withcolor MidnightBlue;
draw (        0.5000000001u,        0.2886751348u)--(        0.3750132153u,        0.3608209299u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2500264668u,        0.2886750720u)--(        0.3750132153u,        0.3608209299u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.5000000001u,        0.2886751348u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2500264668u,        0.2886750720u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500264668u,        0.2886750720u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3750132153u,        0.3608209299u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3750132153u,        0.3608209299u) withcolor MidnightBlue;
draw (        0.3750131973u,        0.5051585020u)--(        0.3749999688u,        0.3608438736u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2500264668u,        0.2886750720u)--(        0.3749999688u,        0.3608438736u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3750131973u,        0.5051585020u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3750131973u,        0.5051585020u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2500264668u,        0.2886750720u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500264668u,        0.2886750720u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3749999688u,        0.3608438736u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3749999688u,        0.3608438736u) withcolor MidnightBlue;
draw (        0.3750131973u,        0.5051585020u)--(        0.2500396952u,        0.4329897004u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2500264668u,        0.2886750720u)--(        0.2500396952u,        0.4329897004u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3750131973u,        0.5051585020u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3750131973u,        0.5051585020u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2500264668u,        0.2886750720u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500264668u,        0.2886750720u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2500396952u,        0.4329897004u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500396952u,        0.4329897004u) withcolor Green;
draw (        0.3750133057u,        0.0721917045u)--(        0.2500397676u,        0.1443604436u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2500264668u,        0.2886750720u)--(        0.2500397676u,        0.1443604436u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3750133057u,        0.0721917045u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3750133057u,        0.0721917045u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2500264668u,        0.2886750720u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500264668u,        0.2886750720u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2500397676u,        0.1443604436u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500397676u,        0.1443604436u) withcolor MidnightBlue;
draw (        0.3750133057u,        0.0721917045u)--(        0.3750000049u,        0.2165063329u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.2500264668u,        0.2886750720u)--(        0.3750000049u,        0.2165063329u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3750133057u,        0.0721917045u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3750133057u,        0.0721917045u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2500264668u,        0.2886750720u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500264668u,        0.2886750720u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.3750000049u,        0.2165063329u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3750000049u,        0.2165063329u) withcolor MidnightBlue;
draw (        0.3750133057u,        0.0721917045u)--(        0.2500265390u,        0.0000458780u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.1250397542u,        0.0721916732u)--(        0.2500265390u,        0.0000458780u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3750133057u,        0.0721917045u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3750133057u,        0.0721917045u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1250397542u,        0.0721916732u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1250397542u,        0.0721916732u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2500265390u,        0.0000458780u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500265390u,        0.0000458780u) withcolor Green;
draw (        0.3750133057u,        0.0721917045u)--(        0.2500265209u,        0.1443374996u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.1250397542u,        0.0721916732u)--(        0.2500265209u,        0.1443374996u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.3750133057u,        0.0721917045u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.3750133057u,        0.0721917045u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1250397542u,        0.0721916732u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1250397542u,        0.0721916732u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2500265209u,        0.1443374996u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500265209u,        0.1443374996u) withcolor MidnightBlue;
draw (        0.2500264668u,        0.2886750720u)--(        0.2500132742u,        0.1443604436u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.1250397542u,        0.0721916732u)--(        0.2500132742u,        0.1443604436u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2500264668u,        0.2886750720u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500264668u,        0.2886750720u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1250397542u,        0.0721916732u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1250397542u,        0.0721916732u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.2500132742u,        0.1443604436u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500132742u,        0.1443604436u) withcolor MidnightBlue;
draw (        0.2500264668u,        0.2886750720u)--(        0.1250529467u,        0.2165063015u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw (        0.1250397542u,        0.0721916732u)--(        0.1250529467u,        0.2165063015u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.2500264668u,        0.2886750720u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.2500264668u,        0.2886750720u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1250397542u,        0.0721916732u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1250397542u,        0.0721916732u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.1250529467u,        0.2165063015u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1250529467u,        0.2165063015u) withcolor Green;
draw (        0.1250397542u,        0.0721916732u)--(        0.0000000000u,        0.0000000000u) withpen pencircle scaled 1.8 dashed evenly withcolor 0.6white;
draw fullcircle scaled 8.0 shifted (        0.1250397542u,        0.0721916732u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.1250397542u,        0.0721916732u) withcolor 0.9red;
draw fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor 0.4white;
fill fullcircle scaled 8.0 shifted (        0.0000000000u,        0.0000000000u) withcolor Yellow;
% r =         0.1443146290
endfig;
end;
