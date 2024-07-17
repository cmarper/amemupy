import json,math
from pprint import pprint
from sl3offsets import sl3_ofssets
import amemu

## Hitograms a given data collection
def histogram(data, bin_width, sides = False):
  histo = {}
  for x in data:
    bin_num = int(math.floor(x/bin_width+0.5))
    r = range(bin_num-1,bin_num+2) if sides else [bin_num]
    for i in r: 
      if i*bin_width not in histo: histo[i*bin_width] = 0
    histo[bin_num*bin_width] += 1
    
  return histo



with open('debug.json') as f:
  data = json.load(f)


diffs_phi = []
diffs_phib = []


for ev_i in range(len(data)):
  ev = data[ev_i]
  input_data = {}

  for (sl, slstart) in [('phi1', 0), ('phi2',4)]:
    input_data[sl]={}
    input_data[sl]['valid'] = [ ev['firmware']['tdc%i'%i] >= 0  for i in range(slstart,slstart+4) ]
    input_data[sl]['wires'] = [ ev['firmware']['wi%i'%i]        for i in range(slstart,slstart+4) ]
    input_data[sl]['t0s']   = [ ev['firmware']['tdc%i'%i]       for i in range(slstart,slstart+4) ]

  (wh,se,st) = (ev['firmware']['wheel'],ev['firmware']['sector'],ev['firmware']['station'])
  sl_phi3_x_offset = sl3_ofssets[(wh,se,st)]

  tps = amemu.chamber_datapath(input_data, sl_phi3_x_offset)
  tps = [ amemu.global_coordinates( tpg, wh, se, st ) for tpg in tps ]

  fw_phi = ev['firmware']['phi']
  fw_phib = ev['firmware']['phib']

  if len(tps)!=1: continue

  amemu_phi = tps[0]['phi_float']
  amemu_phib = tps[0]['phi_bending_float']

  diffs_phi += [amemu_phi - fw_phi]
  diffs_phib += [amemu_phib - fw_phib]
  
  
  #if abs(amemu_phi - fw_phi) > 1e-3:
  #  print ev_i
  #  break

print len(data), len(diffs_phi)

histo = histogram(diffs_phi, 1./2**17)
print histo[0.0]

histo = histogram(diffs_phib, 1./2**11)
print histo[0.0]
