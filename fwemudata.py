## event data
jaime_dump = """
-1 -1 -1
-1 -1 -1
-1 -1 -1
-1 -1 -1
24 78904 0
24 79089 1
24 78842 1
25 78817 0
PosEmul 105.275 PosFW 105.275
TimeEmul 78717 TimeFW 78717
TanPsiEmul -0.654053 TanPsiFW -0.653564
chi2Emul 0.00375 chi2FW 384
"""

jaime_dump = """
Wh:1 Se:4 St:1
1 79660 0
1 79453 1
1 79735 0
1 79393 1
-1 -1 -1
-1 -1 -1
-1 -1 -1
-1 -1 -1
PosEmul 4.5 PosFW 4.5
TimeEmul 79383 TimeFW 79383
TanPsiEmul 0.140381 TanPsiFW 0.140381
chi2Emul 0.00159375 chi2FW 163
"""

jaime_dump = """
Wh:1 Se:6 St:1
17 56499 1
-1 -1 1
17 56487 0
17 56400 1
-1 -1 1
14 56460 0
13 56422 1
13 56642 1
PosEmul 64.325 PosFW 64.325
TimeEmul 56329 TimeFW 56329
TanPsiEmul 0.74707 TanPsiFW 0.747314
chi2Emul 0.053987 chi2FW 4508
"""

jaime_dump = """
Wh:-2 Se:6 St:1
18 79428 0
18 79122 1
18 79404 0
18 79149 1
17 79372 1
18 79165 0
17 79398 1
18 79135 0
PosEmul 76.575 PosFW 76.575
TimeEmul 79077 TimeFW 79077
TanPsiEmul -0.0598145 TanPsiFW -0.0595703
chi2Emul 0.00429255 chi2FW 396
"""
#sl_phi3_x_offset = 21

jaime_dump = """
Wh:0 Se:4 St:2
51 79105 0
51 78834 1
51 78911 0
51 79018 1
54 78932 0
54 79014 1
54 78755 0
55 79051 0
PosEmul 219.675 PosFW 219.675
TimeEmul 78727 TimeFW 78727
TanPsiEmul -0.395996 TanPsiFW -0.395752
chi2Emul 0.013926 chi2FW 1230
"""
#sl_phi3_x_offset = -42

jaime_dump = """
Wh:1 Se:14 St:4
56 56969 0
56 56760 1
56 56988 0
56 56753 1
56 57021 1
57 56715 0
56 57020 1
57 56728 0
PosEmul 235.275 PosFW 235.275
TimeEmul 56676 TimeFW 56676
TanPsiEmul 0.0285645 TanPsiFW 0.0288086
chi2Emul 0.0116131 chi2FW 0.0111914
"""
#sl_phi3_x_offset = -42

jaime_dump = """
Wh:-2 Se:6 St:1
18 79428 0
18 79122 1
18 79404 0
18 79149 1
-1 -1 -1
-1 -1 -1
-1 -1 -1
-1 -1 -1
PosEmul 76.5875 PosFW 76.5875
TimeEmul 79077 TimeFW 79077
TanPsiEmul -0.0595703 TanPsiFW -0.0595703
chi2Emul 0.000458984 chi2FW 0.000107422
"""

jaime_dump = """
Wh:2 Se:12 St:1
-1 -1 -1
-1 -1 -1
-1 -1 -1
-1 -1 -1
36 78691 0
37 78978 0
36 78903 1
37 78763 0
PosEmul 152.044 PosFW 152.044
TimeEmul 78690 TimeFW 78690
TanPsiEmul -0.440674 TanPsiFW -0.440674
chi2Emul 0.00103516 chi2FW 0.00103516
"""

jaime_dump = """
Wh:-1 Se:12 St:1
23 79365 0
23 79278 1
23 79336 0
23 79304 1
-1 -1 -1
-1 -1 -1
-1 -1 -1
-1 -1 -1
PosEmul 97.4875 PosFW 97.4875
TimeEmul 79120 TimeFW 79120
TanPsiEmul -0.0561523 TanPsiFW -0.0561523
chi2Emul 9.76563e-05 chi2FW 9.76563e-05
"""

jaime_dump = """
Wh:2 Se:5 St:1
31 10617 1
32 10868 0
31 10552 0
31 10771 1
-1 -1 -1
-1 -1 -1
-1 -1 -1
-1 -1 -1
PosEmul 30.95625 PosFW 0.
TimeEmul 10509 TimeFW 0
TanPsiEmul 0.32128906 TanPsiFW 0.
chi2Emul 0.00101563 chi2FW 0.
"""

jaime_dump = """
Wh:-2 Se:5 St:1
33 10595 1
34 10810 0
33 10832 1
34 10574 0
-1 -1 -1
-1 -1 -1
-1 -1 -1
-1 -1 -1
PosEmul 37.30625 PosFW 0.
TimeEmul 10569 TimeFW 0
TanPsiEmul -0.49389648 TanPsiFW 0.
chi2Emul 0.00035156 chi2FW 0.
"""


def parse_jaime_dump(jaime_dump):
  jaime_dump = jaime_dump.strip('\n') + '\n'
  
  data = {}

  [wh, se, st] =  [ int(jaime_dump.split('\n')[0].replace(':',' ').split(' ')[i]) for i in [1,3,5]]
  data['wh'] = wh
  data['se'] = se
  data['st'] = st

  hitlines = jaime_dump.split('\n')[1:9]
  for (sl, slstart) in [('phi1', 0), ('phi2',4)]:
    data[sl]={}
    data[sl]['valid'] = [ '-' not in hitlines[i]          for i in range(slstart,slstart+4) ]
    data[sl]['wires'] = [ int(hitlines[i].split(' ')[0])  for i in range(slstart,slstart+4) ]
    data[sl]['t0s']   = [ int(hitlines[i].split(' ')[1])  for i in range(slstart,slstart+4) ]
    data[sl]['lat']   = [ int(hitlines[i].split(' ')[2])  for i in range(slstart,slstart+4) ]

  for (info,t) in [('Pos',float),('Time',int),('TanPsi',float),('chi2',float)]:
    data[info]={
      'FW'  : t(jaime_dump.split(info+'FW '  )[1].split('\n')[0]),
      'Emul': t(jaime_dump.split(info+'Emul ')[1].split(' ' )[0]),
      }

  return data

data = parse_jaime_dump(jaime_dump)