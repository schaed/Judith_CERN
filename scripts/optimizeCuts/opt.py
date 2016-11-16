import ROOT

my = [748,
3430,
7421,
18622,
20874,
21415,
23712,
24202,
25264,
26296,
27368,
30325,
34590,
38125,
38668,
40174,
43723,
48700,
49610,
53600,
58102,
60094,
61291,
61392,
61550,
61707,
62356,
67403,
68238,
70106,
70241,
71884,
72586,
74575,
75266,
75950,
76468,
76479,
79834,
81754,
82076,
88566,
92982,
93421,
93654,
103079,
104906,
106846,
106935,
107308,
108467,
116822,
118197,
123378,
127684,
137673,
143868,
144479,
145509,
150557,
152984,
154536,
154581,
157673,
157988,
158757,
159239,
159688,
159909,
162877,
164708,
166912,
168346,
172303,
172492,
177263,
177515,
180617,
181373,
181578,
181880,
182284,
182285,
182468,
182583,
182934,
183714,
185314,
187363,
189329,
191429,
191443,
193025,
193472,
193677,
193729,
193823,
194984,
195843,
196320,
196988,
197445,
199339,
200330,
200538,
200948,
201181,
201198,
202644,
203433,
203832,
205534,
205934,
207315,
207692,
208186,
209466,
210476,
210867,
212253,
212482,
212681,
212886,
213866,
214684,
215761,
217279,
217512,
217938,
218404,
218488,
220783,
221570,
223461,
225220,
225553,
225680,
225976,
226362,
226809,
228006,
228219,
229895,
230381,
231046,
231113,
231488,
231986,
232484,
233379,
233810,
233910,
234801,
235768,
235947,
237415,
238138,
238577,
238610,
240071,
240290,
241262,
242588,
242969,
243171,
243860,
244531,
245167,
245781,
247538,
247659,
247835,
248991,
249831,
250191,
250323,
251180,
251244,
251742,
251802,
252282,
253412,
253776,
254418,
255749,
256899,
257192,
257236,
257432,
258089,
259378,
259534,
260063,
260496,
260575,
261072,
262324,
262828,
264073,
264108,
265465,
266708,
267407,
267801,
268304,
269067,
271317,
271795,
271896,
271964,
272513,
272545,
272878,
274380,
274753,
275152,
275906,
275927,
276572,
276652,
276856,
277190,
277269,
278549,
278602,
279003,
279073,
280358,
281798,
281909,
282204,
282236,
282542,
283728,
284231,
284252,
284280,
284789,
285069,
289558,
290320,
291018,
291052,
291672,
292803,
293450,
293730,
294264,
295673,
296592,
297502,
297853,
298693,
300071,
301332,]

myb = [129,
204,
351,
358,
539,
636,
680,
1205,
1763,
1824,
1908,
2179,
2188,
2566,
2584,
2866,
3015,
3046,
3188,
3440,
3677,
3890,
4342,
4658,
4879,
5113,
5905,
6028,
6173,
6305,
6337,
6372,
6393,
6699,
6795,
7236,
7271,
7381,
7706,
7759,
7826,
7947,
8658,
8724,
8902,
8922,
9365,
9397,
10322,
10670,
10781,
10909,
10913,
11249,
11626,
11651,
11992,
12291,
12427,
12843,
12989,
13118,
13402,
13776,
13783,
14158,
14269,
14365,
14391,
14924,
14963,
15120,
15804,
15924,
16214,
16345,
16984,
17204,
17865,
17960,
18245,
18765,
19185,
19242,
20715,
20728,
20770,
20836,
21281,
21487,
21718,
21735,
21790,
22407,
22553,
22588,
22597,
22698,
23511,
23554,
23733,
23811,
23904,
23943,
23977,
24082,
24092,
24535,
24637,
24682,
24944,
25112,
25338,
25350,
25867,
26127,
26150,
26239,
26574,
27017,
27119,
27492,
27599,
27624,
27642,
27725,
27737,
27935,
28123,
28157,
28584,
28630,
29043,
29201,
29290,
29389,
29431,
29455,
29545,
30090,
30371,
30474,
30847,
31115,
31177,
31516,
31558,
31606,
31697,
31842,
31932,
32155,
32315,
32880,
33060,
33262,
33350,
33870,
34574,
34647,
34910,
35571,
35920,
36063,
36168,
36576,
36585,
36782,
37022,
37128,
37201,
37284,
37355,
37363,
37398,
37534,
38510,
38639,
38653,
38887,
39574,
39902,
40257,
40505,
40840,
41593,
41650,
41828,
41957,
42129,
42225,
42463,
42647,
43056,
43101,
43179,
43288,
43390,
43674,
43887,
43912,
44092,
44125,
44374,
44442,
44712,
44982,
45004,
45087,
45170,
45363,
45603,
45623,
45758,
45947,
46197,
46211,
46621,
46636,
46647,
46830,
47196,
47228,
47829,
47870,
48022,
48064,
48269,
48352,
48457,
48579,
48609,
48740,
49069,
49078,
49156,
49239,
49330,
49549,
49968,
50079,
50092,
50097,
50445,
50480,
50546,
50672,
50798,
50800,
51276,
51327,
51422,
51725,
51758,
51812,
51825,
51871,
51877,
52178,
52827,
52863,
52919,
53008,
53204,
53397,
53449,
53534,
53616,
53676,
53882,
54140,
54340,
54465,
54614,
54766,
55124,
55210,
55429,
55580,
55593,
55672,
55913,
56074,
56296,
57006,
57157,
57225,
57281,
57324,
57360,
57545,
57561,
57841,
58476,
58697,
58732,
58743,
58950,
59028,
59313,
59459,
59750,
59947,
59963,
60214,
60368,
60372,
60666,
61190,
61282,
61486,
61559,
61646,
61791,
61809,
61970,
62136,
62150,
62201,
62398,
62817,
62878,
63204,
63230,
63319,
63446,
63513,
64038,
64483,
64488,
64766,
65116,
65247,
65572,
65576,
65697,
65742,
65762,
66175,
66215,
66580,
66668,
66788,
66969,
66989,
67294,
67572,
67874,
67895,
68140,
68165,
68232,
68944,
69436,
69754,
69789,
69988,
70067,
70776,
71181,
71343,
71508,
71857,
72007,
72180,
72387,
72528,
72648,
72845,
72886,
72894,
73094,
73407,
73448,
73612,
73675,
74108,
74109,
74498,
74920,
75036,
75058,
75466,
75579,
76039,
76278,
76294,
76942,
77077,
77189,
77262,
77595,
77708,
77778,
77936,
78137,
78315,
78596,
79064,
79140,
79173,
79316,
80450,
80515,
80740,
80854,
81250,
81345,
81558,
81746,
82222,
82881,
83466,
83599,
83996,
84149,
84256,
84327,
84364,
84532,
84705,
84882,
85190,
85208,
85258,
85336,
85343,
85344,
85505,
85573,
85696,
85731,
85940,
85995,
86103,
86123,
86598,
86710,
86976,
87465,
87658,
87747,
87918,
88150,
88204,
88451,
88689,
88770,
88871,
89236,
89237,
89403,
89477,
89588,
89967,
91790,
92038,
92235,
92543,
92653,
92952,
93018,
93243,
93387,
93439,
93451,
93754,
93770,
93809,
93811,
93992,
94097,
94146,
94484,
94582,
94654,
94701,
94757,
94809,
94854,
94896,
95052,
95191,
95215,
95339,
95399,
95739,
95793,
96210,
96381,
96392,
96525,
96655,
96685,
96958,
97054,
97481,
97550,
97585,
97641,
97713,
98169,
98389,
98541,
98594,
98803,
99160,
99801,
100025,
100144,
100154,
100780,
100791,
101024,
101530,
101596,
101611,
101749,
102131,
102291,
102400,
102782,
102820,
103156,
103639,
103708,
104088,
104147,
104413,
104644,
104962,
105058,
105321,
105342,
105419,
105632,
105633,
105651,
105694,
105882,
105981,
106032,
106075,
106242,
106287,
106540,
106803,
106820,
106932,
106959,
106983,
107132,
107381,
107454,
107495,
107731,
107811,
107914,
108557,
108682,
108736,
109059,
109168,
109177,
109672,
109821,
109832,
109838,
109942,
110320,
110339,
110479,
110808,
110967,
111005,
111182,
111322,
111533,
111800,
112237,
112279,
112342,
112469,
112934,
113197,
113505,
113612,
113658,
113750,
113857,
114081,
114120,
114265,
114430,
114502,
114774,
114946,
115005,
115494,
115757,
115892,
116306,
116413,
116491,
116525,
116721,
116776,
116970,
117107,
117147,
117263,
117606,
117680,
117790,
117796,
118078,
118444,
118458,
118668,
118933,
119101,
119330,
119356,
119383,
119495,
119530,
119538,
119666,
120531,
120884,
121159,
121429,
122459,
122567,
122624,
122879,
122904,
122994,
123101,
123440,
123494,
123756,
123775,
123834,
123906,
124028,
124040,
124140,
124225,
124355,
124660,
124880,
124954,
125104,
125284,
125326,
125376,
125444,
125486,
125977,
126041,
126042,
126204,
126273,
126444,
126628,
126920,
127163,
127227,
127375,
127680,
127707,
127724,
128042,
128122,
128143,
128431,
128491,
128495,
128557,
129145,
129202,
129423,
129884,
130044,
130321,
130377,
130404,
130550,
130620,
130662,
130989,
131151,
131507,
131865,
131901,
131966,
132150,
132348,
132394,
132570,
132734,
132758,
132975,
133057,
133611,
133640,
133695,
133886,
133912,
134199,
134406,
134713,
134767,
134790,
134793,
134944,
135180,
135216,
135262,
135380,
135394,
136328,
136432,
136864,
136935,
136936,
137280,
137753,
137759,
137791,
137878,
137903,
138201,
138588,
138596,
138643,
139465,
139714,
139767,
139794,
139862,
139995,
139998,
140328,
140533,
140789,
140831,
141091,
141181,
141726,
141805,
141942,
142071,
142125,
142324,
142594,
142614,
142942,
142954,
142985,
143088,
143273,
143292,
143520,
143652,
143718,
143800,
143824,
143972,
144001,
144233,
144359,
144419,
144911,
145115,
145701,
146253,
146436,
146441,
146559,
146935,
147048,
147090,
147361,
147593,
147712,
147718,
147901,
147970,
148168,
148285,
148557,
148643,
149056,
149176,
149231,
149339,
150548,
150584,
150898,
150926,
150963,
151308,
151355,
151409,
151484,
151586,
151674,
151679,
151710,
151737,
152137,
152354,
152369,
152696,
152980,
153259,
153807,
154057,
154547,
154934,
155076,
155822,
155951,
155955,
156040,
156044,
156046,
156267,
156585,
156966,
157065,
157143,
157187,
157201,
157420,
157525,
157725,
157898,
158081,
158094,
158304,
158379,
158646,
158680,
158694,
158716,
158962,
159003,
159009,
159124,
159455,
159500,
159553,
159674,
160007,
160235,
160392,
160439,
160492,
160719,
160739,
160753,
160856,
160992,
160996,
161068,
161116,
161222,
161529,
161838,
161839,
162002,
162151,
162241,
162272,
162504,
162683,
162692,
162913,
162997,
163163,
163498,
163519,
164022,
164502,
164616,
164863,
164987,
165178,
165219,
165467,
165612,
165655,
165737,
165935,
166070,
166336,
166351,
166386,
166415,
166474,
166495,
166741,
166985,
167145,
167325,
167662,
167993,
168599,
168713,
169168,
169345,
170076,
170469,
170483,
170511,
170654,
171232,
171332,
171459,
171778,
171929,
172360,
173394,
173506,
173565,
173721,
173741,
173742,
173745,
173996,
174045,
174500,
174546,
175006,
175134,
175155,
175256,
175342,
175358,
175438,
175570,
175713,
175784,
175794,
175849,
176021,
176342,
176665,
176684,
176745,
176836,
177032,
177092,
177351,
177485,
177519,
177564,
177576,
177667,
177976,
178056,
178315,
178458,
178503,
178549,
178726,
178728,
178776,
178829,
179003,
179114,
179338,
179416,
179442,
179541,
179780,
179836,
180111,
180181,
180220,
180248,
180370,
180481,
180731,
180795,
181453,
181641,
182148,
182187,
182567,
182667,
182748,
182751,
182987,
182988,
183405,
183545,
183860,
183871,
183954,
183994,
184003,
184410,
184817,
184832,
184877,
184899,
184985,
185086,
185266,
185298,
185542,
185765,
185876,
186513,
186777,
186911,
187044,
187352,
187417,
187509,
187647,
187749,
187863,
188156,
188722,
189003,
189077,
189150,
189500,
189588,
190403,
190504,
190811,
190886,
191130,
191830,
191921,
192044,
192117,
192190,
192295,
192347,
192422,
192666,
192782,
192923,
193250,
193468,
193643,
193685,
194088,
194371,
194491,
194538,
194657,
194802,
194972,
195000,
195046,
195100,
195117,
195174,
195355,
195384,
195543,
195586,
195665,
195873,
195951,
195959,
196121,
196455,
196460,
196595,
196682,
196890,
197150,
197220,
197261,
197674,
198036,
198494,
198547,
198594,
198628,
198800,
199088,
199115,
199627,
199649,
199892,
199906,
200004,
200353,
200437,
200696,
200752,
201059,
201211,
201343,
201541,
201602,
201655,
201829,
201900,
202089,
202462,
202649,
203020,
203127,
203208,
203214,
203370,
203431,
203545,
203663,
203778,
203787,
204065,
204243,
204407,
204724,
204817,
204857,
205023,
205314,
205671,
205869,
205870,
205998,
206088,
206134,
206217,
206312,
206824,
206933,
207073,
207090,
207125,
207384,
207448,
207570,
208122,
208385,
208516,
208740,
208835,
208952,
209027,
209079,
209155,
209230,
209337,
209348,
209361,
209369,
209397,
209405,
209475,
209588,
209680,
209732,
209770,
209936,
209941,
210115,
210254,
210323,
210355,
210435,
210514,
210568,
210653,
210756,
210858,
210954,
211077,
211116,
211245,
211539,
211574,
212007,
212042,
212214,
212710,
212778,
213081,
213525,
213548,
213662,
213956,
214278,
214645,
214733,
214761,
214812,
214815,
214842,
214974,
215185,
215260,
215269,
215284,
215342,
215684,
215992,
216001,
216110,
216227,
216292,
216318,
216508,
216543,
216607,
216632,
216694,
216818,
216914,
217061,
217066,
217080,
217081,
217131,
217137,
217215,
217260,
217328,
217633,
217812,
217891,
217978,
217990,
218097,
218464,
218608,
219273,
219305,
219448,
219629,
220138,
220163,
220332,
220429,
220755,
221082,
221305,
221365,
221502,
221563,
221564,
221712,
221878,
222032,
222253,
222271,
222497,
222679,
223248,
223362,
223527,
223606,
223733,
223876,
223972,
224242,
224417,
224476,
224601,
224690,
224723,
224994,
225288,
225647,
225766,
225768,
225859,
225873,
225989,
226132,
226192,
226234,
226423,
226474,
226799,
226943,
227013,
227215,
227347,
227440,
227742,
227778,
227789,
227842,
228017,
228501,
228565,
228580,
228800,
228930,
229030,
229045,
229252,
229500,
229510,
229729,
229740,
229951,
230049,
230435,
230800,
231115,
231197,
231686,
231887,
232049,
232376,
232655,
232693,
232983,
233063,
233262,
233418,
233506,
233902,
234024,
234135,
234151,
234861,
234945,
235083,
235282,
235377,
235622,
235827,
235830,
235899,
235956,
236148,
236277,
236769,
236874,
236960,
237256,
237346,
237507,
237517,
237668,
237890,
237913,
238267,
238689,
238853,
239422,
239597,
239791,
239926,
239937,
240102,
240321,
240576,
240954,
241656,
241765,
241852,
241905,
241996,
241999,
242069,
242673,
242704,
242836,
242857,
242912,
243001,
243020,
243102,
243181,
243395,
243718,
243802,
244368,
244405,
245231,
245330,
245371,
245587,
245610,
245625,
246270,
246346,
246552,
246754,
247039,
247258,
247545,
247569,
247833,
247879,
248041,
248296,
248328,
248434,
248472,
248516,
248752,
248775,
248934,
249210,
249332,
249534,
249579,
249606,
249774,
249903,
250079,
250104,
250160,
250163,
250514,
250728,
251028,
251189,
251400,
251450,
252185,
252191,
252386,
252448,
252747,
253008,
253083,
253259,
253361,
253514,
253585,
253684,
254226,
254256,
254282,
254292,
254386,
254439,
254688,
254934,
255215,
255237,
255394,
255841,
255919,
255929,
255978,
255988,
256310,
256947,
257425,
257434,
258027,
258144,
258467,
258489,
258502,
258583,
258652,
258764,
258785,
259000,
259037,
259090,
259263,
259365,
259443,
259449,
259506,
259787,
259808,
259981,
260045,
260179,
260192,
260206,
260328,
260734,
261030,
261392,
261429,
261691,
261972,
262084,
262328,
262415,
262577,
262578,
262605,
262617,
262633,
262655,
262725,
262973,
263052,
263384,
263479,
263624,
263679,
263919,
263951,
264089,
264156,
264172,
264289,
264315,
264488,
264501,
264607,
264677,
264772,
264965,
265183,
265240,
265321,
265413,
265757,
265776,
265958,
266152,
266511,
266526,
266710,
267272,
267467,
267558,
267636,
267681,
267828,
268181,
268628,
268885,
268958,
269065,
269537,
269655,
269851,
269914,
269918,
269961,
269988,
270011,
270096,
270108,
270151,
270372,
270754,
271013,
271790,
272347,
272444,
272556,
272577,
272767,
272772,
272814,
272826,
272941,
272948,
273179,
273207,
273494,
273546,
273911,
273930,
273980,
274119,
274179,
274438,
274544,
274607,
274808,
274816,
274977,
275093,
275118,
275470,
275681,
275926,
275986,
276081,
276094,
276227,
276291,
276442,
276476,
276812,
277130,
277208,
277535,
277859,
277887,
278180,
278213,
278486,
278507,
278564,
278606,
278772,
278796,
278996,
279033,
279242,
279479,
279813,
280064,
280264,
280397,
280429,
280697,
280746,
280779,
280821,
281005,
281089,
281100,
281114,
281246,
281268,
281287,
281371,
281380,
281516,
281994,
282300,
282422,
282431,
282629,
282817,
282845,
282853,
283181,
283541,
283555,
283593,
283878,
284027,
284128,
284154,
284253,
284301,
284723,
284763,
284983,
285589,
285833,
286328,
286799,
286929,
287019,
287091,
287177,
287483,
287862,
288058,
288191,
288554,
288724,
288894,
289192,
289234,
289290,
289624,
289781,
289799,
289928,
289972,
290001,
290060,
290566,
290726,
290744,
290759,
290764,
290936,
291080,
291749,
291775,
291948,
292008,
292084,
292408,
292447,
292618,
292816,
293045,
293100,
293207,
293208,
293656,
293720,
294008,
294238,
294440,
294533,
294565,
294714,
294947,
295142,
295280,
295379,
295543,
295879,
295927,
295933,
296546,
296820,
296920,
297125,
297932,
298411,
298471,
298746,
298858,
299455,
299529,
299531,
299861,
300537,
300660,
300705,
300972,
301291,
301368,
301481,
301484,]
class Pixel:

    def __init__(self, num=0, cut=''):

        self.num=num
        self.hCharge = ROOT.TH1F('charge%s_%s' %(cut,num),'charge%s_%s' %(cut,num),2000,0.0,0.2)
        self.hT0 = ROOT.TH1F('t0%s_%s' %(cut,num),'t0%s_%s'%(cut,num),1000,0.0,1000.0)
        #self.hTime = ROOT.TH1F('timing%s_%s' %(cut,num),'timing%s_%s' %(cut,num),20000,-10000.0,10000.0)
        self.hTime = ROOT.TH1F('timing%s_%s' %(cut,num),'timing%s_%s' %(cut,num),500,0.0,500.0)
        self.hMag = ROOT.TH1F('fftmag%s_%s' %(cut,num),'fftmag%s_%s' %(cut,num),500,0.0,15.0)
        self.hPhs = ROOT.TH1F('fftphs%s_%s' %(cut,num),'fftphs%s_%s' %(cut,num),64,-3.2,3.2) 

        self.hCharge.GetXaxis().SetTitle('Charge [V]')
        self.hCharge.GetYaxis().SetTitle('Events')
        self.hT0.GetXaxis().SetTitle('Time after Trigger [ns]')
        self.hT0.GetYaxis().SetTitle('Events')
        self.hTime.GetXaxis().SetTitle('Charge Collection Time [ns]')
        self.hTime.GetYaxis().SetTitle('Events')
        self.hMag.GetXaxis().SetTitle('FFT Magnitude of low freq. [mV]')
        self.hMag.GetYaxis().SetTitle('Events')
        self.hPhs.GetXaxis().SetTitle('FFT Phase of low freq. [rad]')
        self.hPhs.GetYaxis().SetTitle('Events')
    def Fill(self, ch, t0, timing, mag, phs):
        self.hCharge.Fill(ch)
        self.hT0.Fill(t0)
        self.hTime.Fill(timing)
        self.hMag.Fill(mag)
        self.hPhs.Fill(phs)
    def Write(self,fout):
        self.hCharge.SetDirectory(fout)
        self.hT0.SetDirectory(fout)
        self.hTime.SetDirectory(fout)
        self.hMag.SetDirectory(fout)
        self.hPhs.SetDirectory(fout)
        #fout.Write()
        #fout.Close()
        
#f = ROOT.TFile.Open('../../TowerJazz/run_904_tj_W3R13_50um_6V_1DRS.root')
#f = ROOT.TFile.Open('../../TowerJazz/run_804_tj_W3R15_50um_6V.root')
#f = ROOT.TFile.Open('../../TowerJazz/run_804_tj_W3R15_50um_6V_new.root')
#f = ROOT.TFile.Open('../../TowerJazz/run_868_tj_W3R13_50um_6V_1DRS_new.root')
#f = ROOT.TFile.Open('../../TowerJazz/test_lowpass.root')
#f = ROOT.TFile.Open('../../TowerJazz/test_new.root')
#f = ROOT.TFile.Open('../../TowerJazz/test_new.root')
#f = ROOT.TFile.Open('../../TowerJazz/run804_new5_sync-analysis-result.root')
#f = ROOT.TFile.Open('../../TowerJazz/run_32877_tj_W3R13_50um_6V_3DRS.root')
#fout = ROOT.TFile.Open('run_32877_tj_W3R13_50um_6V_3DRS_fft-analysis-result.root','RECREATE')
#f = ROOT.TFile.Open('../../TowerJazz/run_33043_tj_W3R13_50um_6V_3DRS_fft.root')
#f = ROOT.TFile.Open('../../TowerJazz/run_33192_tj_W3R15_M129_6V_3DRS_fft10.root')
#f = ROOT.TFile.Open('../../TowerJazz/run_33213_tj_W3R15_M108_6V_3DRS_fft12.root')
f = ROOT.TFile.Open('../../TowerJazz/run_33192_tj_W3R15_M129_6V_3DRS_fft10.root')
#f = ROOT.TFile.Open('../../TowerJazz/test_3147.root')
fout = ROOT.TFile.Open('run_33192_tj_sig.root','RECREATE')

hits = f.Get('Plane0/Hits')
events = f.Get('Event')
fnum=[]
for i in events:
    fnum+=[i.FrameNumber]
ii=0
pixels=[]
pixels_15V=[]
pixels_15V_time=[]
pixels_15V_t0=[]
pixels_FFT2=[]
for h in hits:

    if (ii%10000)==0:
        print 'Event: ',ii

    #if fnum[ii]==7421 or fnum[ii]==3430 or fnum[ii]==748 or fnum[ii] in my:
    #if not (fnum[ii] in myb):
    #    ii+=1
    #    continue
    #if fnum[ii]==636 or fnum[ii]==1763 or fnum[ii]==748:    
    #    print ii,fnum[ii]
    #    for p in range(0,h.NHits):
    #        print '     ',p, h.Value[p], h.T0[p], h.Timing[p], h.LowFreqFFT[p], h.LowFreqFFTPhase[p], h.IsHit[p], h.ValidFit[p]
    ii+=1
    #if h.NHits==4:
    #    if h.Value[3]<0.0:
    #        continue
    #print 'Event: ',ii
    if len(pixels)==0 and h.NHits>0:
        for p in range(0,h.NHits):
            pixels+=[Pixel(p)]
            pixels_15V+=[Pixel(p,'_Val15V')]
            pixels_15V_time+=[Pixel(p,'_Val15VTime')]
            pixels_15V_t0+=[Pixel(p,'_Val15VT0')]
            pixels_FFT2+=[Pixel(p,'_FFT2')]            
    for p in range(0,h.NHits):
        if h.Value[p]<0.0005:
        #if h.Value[p]<0.08:
            continue
        #print fnum[ii]
        #print 'pixel: ',p,' Charge: ',h.Value[p], h.T0[p], h.Timing[p]
        #print 'Charge: ',h.Value[p]
        # Cuts:
        #if h.Value[p]<0.05:
        #    continue
        #if not (h.Timing[p]<380.0 and h.Timing[p]>370.0):
        #if not (h.Timing[p]<220.0 and h.Timing[p]>180.0):
        #    continue
        #if not (h.Timing[p]<230.0 and h.Timing[p]>218.0):
        #    continue
        #print ' PASS Charge pixel: ',p,' Charge: ',h.Value[p], h.T0[p], h.Timing[p]
        pixels[p].Fill(h.Value[p], h.T0[p], h.Timing[p], h.LowFreqFFT[p], h.LowFreqFFTPhase[p])
        if h.Value[p]>0.001: #0.003
            pixels_15V[p].Fill(h.Value[p], h.T0[p], h.Timing[p], h.LowFreqFFT[p], h.LowFreqFFTPhase[p])
            if h.Timing[p]>15.0 and  h.Timing[p]<60.0:
                pixels_15V_time[p].Fill(h.Value[p], h.T0[p], h.Timing[p], h.LowFreqFFT[p], h.LowFreqFFTPhase[p])
            if (h.T0[p]<280.0 and h.T0[p]>200.0) and (h.Timing[p]>15.0 and  h.Timing[p]<60.0): 
                pixels_15V_t0[p].Fill(h.Value[p], h.T0[p], h.Timing[p], h.LowFreqFFT[p], h.LowFreqFFTPhase[p])
            if (h.LowFreqFFT[p]>2.0 and h.T0[p]>5.0) and (h.Timing[p]>10.0 and  h.Timing[p]<60.0): 
                pixels_FFT2[p].Fill(h.Value[p], h.T0[p], h.Timing[p], h.LowFreqFFT[p], h.LowFreqFFTPhase[p])

for p in pixels:
    p.Write(fout)
for p in pixels_15V:
    p.Write(fout)
for p in pixels_15V_time:
    p.Write(fout)
for p in pixels_15V_t0:
    p.Write(fout)
for p in pixels_FFT2:
    p.Write(fout)        
fout.Write()
fout.Close()
print 'Done'

