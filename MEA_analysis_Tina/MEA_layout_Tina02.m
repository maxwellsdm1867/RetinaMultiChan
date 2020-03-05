function channel = MEA_layout(index)
% This function is just to assign the order to the corresponding channel of
% the "MEA-200", you can find the layout number under 'C:\Program Files (x86)\Multi Channel Systems\MC_Rack\MeaLayouts'
% INPUT
% index: the order of the channel
% OUTPUT
% channel: the layout positon of the corresponding channel
layout = [21
19
16
15
12
10
24
22
20
17
14
11
9
7
26
25
23
18
13
8
6
5
29
30
28
27
4
3
1
2
32
31
33
34
57
58
60
59
35
36
38
43
48
53
55
56
37
39
41
44
47
50
52
54
40
42
45
46
49
51
];
 channel = layout(index);
