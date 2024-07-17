import os, math
from pprint import pprint


###################
## CONGIFURATION ##
###################
params = {}

params['verbose_latcombs']  = 0
params['verbose_path']      = 0

# Superlayer level, Maximun chi2 value for a valid superlayer segment.
# 16 bits, mm^2 / 64
params['chisqr_threshold'] = 64

# Superlayer level, maximum threshold for tan(phi) (absolute value). 
# 14 bits (max 16383)
params['tanphi_x4096_threshold'] = 4096

# In the comparisons to select a segment by its chi2 value, use the perfect value instead of the rounded value
# (obviously, this will stop perfectly emulating the fw)
params['compare_by_perfect_chi2'] = False

def setParam(name, value):
  global params
  params[name] = value

def getParam(name):
  return params[name]

########################
## CONSTANTS AND LUTS ##
########################

## Important: simulate the 2-bits or 4-bits precision trigger module
PRECISION = 4 # 4 bits or 2 bits

TDCTIME_REDUCED_SIZE = 10
MAX_DRIFT_TIME_F = 386.75
MAX_DRIFT_TIME_P = 387
CELL_LENGTH_P = 42 * (2**PRECISION)
CELL_LENGTH   = 42
CELL_HEIGHT = 13
DRIFT_SPEED_P = int(55.5 * 4) if PRECISION == 2 else 889
DRIFT_SPEED_F   = CELL_LENGTH / MAX_DRIFT_TIME_F / 2 # 54.3
TAN_DIVISION_DENOMINATOR_BITS = 16
Z_FACTOR_SEGM = [ CELL_HEIGHT * mult for mult in [-6,-2,2,6] ]
LAYERS = range(4)
VERT_PHI1_PHI3_INV_Q15 = int(139.5 * 4)
CH_CENTER_TO_MID_SL_P = int(117.5 * 4) # defined as 470 in fw
Z_FACTOR_CORR = [ CELL_HEIGHT * [-6,-2,2,6][i%4] + CH_CENTER_TO_MID_SL_P * (-1 if i<4 else +1) for i in range(8) ]

NULL_SEGMENT = {
  'bx_id_prec'                : 0,
  'bx_id'                     : 0,
  'chi_square'                : 0,
  'delta_hit_x'               : [0]*4,
  'horizontal_position_float' : 0.,
  'phi_tangent_float'         : 0.,
  'main_params'               : {
    'bx_time'                 : 0,
    'quality'                 : 0,
    'laterality_combination'  : '0000',
    'hits'                    : [  {
        'super_layer_id'      : 0,
        'layer_id'            : 0,
        'channel_id'          : 0,
        'dt_time'             : { 'valid' : 0, 'tdc_time' : 0 },
      } for i in LAYERS ]
    }
  }


## for direct fit of SL segment
from generate_encoder_vhdl import encoder_dict_slvs as LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER
from generate_encoder_vhdl import encoder_dict_real as LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_FLOAT
from generate_encoder_vhdl import INCREASED_RES, NUMER_COEFF_SIZE, NUMER_CONST_SIZE, DIV_SHR_BITS

to_two_comp   = lambda val,size: val%(2**size)
from_two_comp = lambda val,size: val - ((2*val)&(1<<size))

NULL_ENCODER_LIST_SLVS = [ {mag: {'coeff': [ 0 ]*4, 'add': 0, 'mult': 0 } for mag in INCREASED_RES } ] * 6
NULL_ENCODER_LIST_REAL = [ {mag: (0., [ 0. ]*4, 0.) for mag in INCREASED_RES } ] * 6

# Global coordinates conversion
X_SIZE       = 18
TANPSI_SIZE  = 15

PHI_LUT_ADDR_WIDTH = 12
PHI_B_SHL_BITS = 7
PHI_MULT_SHR_BITS = 10
PHI_LUT_A_BITS = 12
PHI_LUT_B_BITS = 20

PHIB_LUT_ADDR_WIDTH = 9
PHIB_B_SHL_BITS = 7
PHIB_MULT_SHR_BITS = 10
PHIB_LUT_A_BITS = 10
PHIB_LUT_B_BITS = 16

PHI_PHIB_RES_DIFF_BITS = 6


###########################
## ENTITIES TO FUNCTIONS ##
###########################

def segment_fitter(valid, t0s, cell_horiz_layout, latcomb_consts, latcomb_consts_real, xwire_mm, coarse_pos, coarse_offset, verbose):
  """returns dict with result of laterality combination analysis"""
  (latcomb, constants ) = latcomb_consts
  (dummy, constants_real ) = latcomb_consts_real
  
  lat_array = [ 1 if (latcomb>>i)&1 else -1 for i in LAYERS ]
  
  result = {}
  result_perfect = {}
  for mag in constants:
    result[mag] = (
        ( 
          sum([ from_two_comp(int(constants[mag]['coeff'][i],2),NUMER_COEFF_SIZE[mag])*t0s[i] for i in LAYERS]) * (2**INCREASED_RES[mag]) +
          from_two_comp(int(constants[mag]['add'],2),NUMER_CONST_SIZE[mag])
        ) * int(constants[mag]['mult'],2)
      ) >> DIV_SHR_BITS[mag]
    result_perfect[mag] = (sum([ constants_real[mag][1][i]*t0s[i] for i in LAYERS]) +  constants_real[mag][0] ) / constants_real[mag][2] if constants_real[mag][2] else 0
  bx_time = result['t0'] + (coarse_offset<<(TDCTIME_REDUCED_SIZE-1))
  bx_time_perfect = result_perfect['t0'] + (coarse_offset<<(TDCTIME_REDUCED_SIZE-1))
  bx_id = ( ( ( bx_time + 13 ) * 41943 ) >> 20 ) % 3564
  
  pos = result['pos'] + coarse_pos
  pos_perfect = result_perfect['pos'] + coarse_pos/16.
  
  segm_valid = 1
  
  delta_hit_x = [0]*4
  
  if abs(result['slope']) > params['tanphi_x4096_threshold'] or latcomb == 0:
    segm_valid = 0
  
  chi2_mm2_p = 0
  chi2_mm2_perfect = 0
  for i in LAYERS:
    drift_time = t0s[i]-result['t0']
    if valid[i] and ( drift_time < 0 or drift_time > MAX_DRIFT_TIME_P ):
      segm_valid = 0
    
    drift_dist = (( (drift_time * 2**4 + 9) * 445 )>>13)
    
    xdist = xwire_mm[i]*2**4 - result['pos'] + lat_array[i] * drift_dist

    delta_hit_x[i] = xdist
    
    xdist -= (3-2*(3-i))*result['slope_xhh']
    res = xdist*(1 if valid[i] else 0)
    res2 = res*res
    chi2_mm2_p += res2 * 4

    ## Perfect chi2 calculated for perfect path
    xdist_perfect = xwire_mm[i] - result_perfect['pos'] + lat_array[i] * (t0s[i]-result_perfect['t0']) * DRIFT_SPEED_F - (3-2*(3-i))*result_perfect['slope_xhh']
    ## Perfect chi2 calculated for rounded path
    #xdist_perfect = xwire_mm[i] - result['pos']/16. + lat_array[i] * (t0s[i]-result['t0']) * DRIFT_SPEED_F - (3-2*(3-i))*result['slope_xhh']/16.
    res_perfect = xdist_perfect*(1 if valid[i] else 0)
    res2_perfect = res_perfect*res_perfect
    chi2_mm2_perfect += res2_perfect
    
  if chi2_mm2_p > params['chisqr_threshold'] * 2**4:
    segm_valid = 0
  
  segm = {}
  
  segm['segm_valid'] = segm_valid
  
  segm['lat'      ] = lat_array
  segm['latcomb'  ] = latcomb
  segm['q'        ] = 0 if not segm_valid else 4 if all(valid) else 2
  segm['t'        ] = bx_time
  segm['t_perfect'] = bx_time_perfect
  segm['valids'   ] = [valid[i] for i in LAYERS ]
  
  
  segm['xcoor_mm_p']                = pos
  segm['xcoor']                     = pos / 1. / 2**PRECISION
  segm['xcoor_perfect']             = pos_perfect
  segm['horizontal_position_float'] = segm['xcoor'] # ab7manager format compatibility
  
  segm['tanphi_x4096']              = result['slope']
  segm['tanphi']                    = - result['slope'] / 4096.
  segm['tanphi_perfect']            = - result_perfect['slope']
  segm['phi_tangent_float']         = - segm['tanphi'] # ab7manager format compatibility
  
  segm['chi2_mm2_p']    = chi2_mm2_p
  segm['chi2_cm2']      = chi2_mm2_p /102400.
  segm['chi2_perfect']  = chi2_mm2_perfect /100.
  segm['delta_hit_x_p']       = delta_hit_x
  segm['delta_hit_x_perfect'] = delta_hit_x

  segm['chi_square']    = chi2_mm2_p            # ab7manager format compatibility
  segm['delta_hit_x']   = delta_hit_x           # ab7manager format compatibility
  
  segm['bx_id']         = bx_id
  segm['bx_id_p']       = 0
  segm['bx_id_perfect'] = bx_id # should be calculated without approximations
  segm['bx_id_prec']    = 0 # ab7manager format compatibility
  
  #pprint(segm)
  
  ## print results to do... some day
  #if verbose:
  #  print '    - laterality combination:', lat_array, 'T0 = %i'%(bx_time)
  #  print 'Verbose values for path calculation', latcomb
  #  print '  xcoor_mm_p', xcoor_mm_p, segm['xcoor']
  #  print '  tanphi_x4096', tanphi_x4096, segm['tanphi']
  #  print '  chi2_mm2_p', chi2_mm2_p, '(perfect =', chi2_perfect, ')'
  #  
  #  print '  drift_dist_um_p', [x for x in drift_dist_um_p]
  #  print '    converted to mm -->', [x/4./1.024 for x in drift_dist_um_p]
  #  print '  wirepos_mm_p', [x for x in wirepos_mm_p]
  #  print '    converted to mm -->', [x/4. for x in wirepos_mm_p]
  #  print '  pos_mm_p', [x for x in pos_mm_p]
  #  print '    converted to mm -->', [x/4. for x in pos_mm_p]

  return segm

  

# in fact analyzer + latcomb_trunk 
def analyzer(valid, wires, t0s, cell_horiz_layout, verbose): # also includes the analyzer_bx_calculation vhd module
  """returns a list of segment candidates"""
  if verbose:
    print '  Calculating laterality combinations'
    print '    cell_horiz_layout = ', cell_horiz_layout

  # calculate the coarse offset position
  tmp = 1 if valid[1] else 3
  coarse_pos = (wires[tmp]*2 - cell_horiz_layout[tmp]) * 21 * 2**4
  
  # calculate the relative position of wires in mm wrt layer 0's cell wire
  xwire_mm = [ 21*cell_horiz_layout[i] for i in LAYERS ]
  
  # divide the timestamps in coarse + reduced part
  valid_coarse_times = [ ( t0s[i]>>(TDCTIME_REDUCED_SIZE-1) ) for i in LAYERS if valid[i] ]
  max_coarse_time = max( valid_coarse_times )
  min_coarse_time = min( valid_coarse_times )
  
  if max_coarse_time - min_coarse_time >= 2: # make the segment candidate invalid
    return []
  
  coarse_offset = max_coarse_time - 1
  
  reduced_times = [ # omg this is so much simpler in vhdl
    ( (1-((max_coarse_time&1)^((t0s[i]>>(TDCTIME_REDUCED_SIZE-1))&1)))<<(TDCTIME_REDUCED_SIZE-1) ) + 
    ( t0s[i] & int('1'*(TDCTIME_REDUCED_SIZE-1),2) )
    for i in LAYERS ]

  ## In the next two assignments, the else condition should never really be met
  latcomb_consts_arr = (
    LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER[(cell_horiz_layout,tuple(valid))]
    if (cell_horiz_layout,tuple(valid)) in LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER else
    NULL_ENCODER_LIST_SLVS
    )
  latcomb_consts_arr_real = (
    LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_FLOAT[(cell_horiz_layout,tuple(valid))]
    if (cell_horiz_layout,tuple(valid)) in LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_FLOAT else
    NULL_ENCODER_LIST_REAL
    )
  segment_candidates = [
    segment_fitter(valid, reduced_times, cell_horiz_layout, latcomb_consts_arr[i], latcomb_consts_arr_real[i], xwire_mm, coarse_pos, coarse_offset, verbose)
    for i in range(len(latcomb_consts_arr))
    ]

  segment_candidates = [ segm for segm in segment_candidates if segm['segm_valid'] ]

  for segm in segment_candidates:
    segm['q'] = segm['q'] - (0 if len(segment_candidates) == 1 else 1)
    segm['main_params'] = { # ab7manager format compatibility
      'bx_time': segm['t'],
      'quality': segm['q'],
      'laterality_combination': bin(16+segm['latcomb'])[3:],
      'hits': [ { 'super_layer_id':0,
                  'layer_id':i,
                  'channel_id':wires[i],
                  'dt_time': {'valid':segm['valids'][i], 'tdc_time':t0s[i]}
                } for i in LAYERS ]
      }
      
  return segment_candidates
        
def seg_duplicated_shifting_filter( segment_candidates, t0s, wires ):
  result = []
  for segm in segment_candidates:
    not_a_dupe = True
    for pre_segm in result:
      for i in LAYERS:
        if  ( ( t0s[i] == t0s[i] and
                wires[i] == wires[i] and
                segm['lat'][i] == pre_segm['lat'][i] and
                segm['valids'][i] and
                pre_segm['valids'][i]
              ) or (
                not segm['valids'][i] and not pre_segm['valids'][i]
              ) ):
          not_a_dupe = False
    if not_a_dupe:
      result += [segm]
  
  return result

def chain_qenhancer_filter( segment_candidates ):
  result = []
  chi2_field = 'chi2_perfect' if params['compare_by_perfect_chi2'] else 'chi2_mm2_p'
  for segm in segment_candidates:
    if not result or (segm['q']==3 and segm[chi2_field] < result[-1][chi2_field]):
      result = [ segm ]
    elif segm['q'] in [1,2]:
      result += [ segm ]
  return result

def superlayer_datapath( wires, t0s, valid ):
  ### PSEUDO-MIXER: operations performed by the vhdl mixer part
  # Pseudo-mixer: only process if 3 hits or more
  if valid.count(True)<3 : return []
  # Pseudo-mixer: complete missing wires for selecting a valid cell_horiz_layout
  if   wires[0]<0: wires[0] = wires[1]
  elif wires[1]<0: wires[1] = wires[0]
  elif wires[2]<0: wires[2] = wires[1]-1
  elif wires[3]<0: wires[3] = wires[2]
  ## Pseudo-mixer: obtain cell_horiz_layout
  cell_horiz_layout = [ (wires[i]-wires[0])*2 - (1 if i in [1,3] else 0) for i in LAYERS ] #LAYERS 2 and 4 are slightly to the left
  cell_horiz_layout = tuple( cell_horiz_layout )

  #### ANALYZER #### (analysis_datapath.vhd)
  
  ## ANALYZER DIRECT HIT
  segment_candidates = analyzer(valid, wires, t0s, cell_horiz_layout, params['verbose_latcombs'])

  #### SEG_DUPLICATED_SHIFTING_FILTER ####
  #segment_candidates = seg_duplicated_shifting_filter( segment_candidates, t0s, wires )
  #pprint([(seg['q'], seg['lat'], seg['t_perfect']) for seg in segment_candidates])
  
  #### CHAIN_QENHANCER_FILTER ####
  segment_candidates = chain_qenhancer_filter( segment_candidates )

  return segment_candidates

def path_calc_bx( bx_time ):
  bx_id_long = ((bx_time<<3) * 41943 ) >> 20 ## lut value 41943 = int(round(2**20/25.,0))
  bx_id = (bx_id_long>>3)+1 if (bx_id_long>>2)&1 else (bx_id_long>>3)
  bx_id_perfect = int( round( bx_time / 25.,0 ) )
  return (bx_id, bx_id_long & 0x7, bx_id_perfect)

def primitive_path_calculator( segments, sl_phi3_x_offset, verbose ):
  if verbose:
    print '#'*40 + '\nCalculating primitive\n'

  [phi1, phi2] = segments

  ## calc_primitive_tan_pos_bx_time.vhd
  phi3_offset_p = sl_phi3_x_offset << PRECISION
  diff_x_pos          = phi2['xcoor_mm_p'   ] + phi3_offset_p    - phi1['xcoor_mm_p'   ]
  diff_x_pos_perfect  = phi2['xcoor_perfect'] + sl_phi3_x_offset - phi1['xcoor_perfect']

  # position
  base_x          = phi1['xcoor_mm_p'   ] - ( phi3_offset_p     if phi3_offset_p    < 0 else 0 )
  base_x_perfect  = phi1['xcoor_perfect'] - ( sl_phi3_x_offset  if sl_phi3_x_offset < 0 else 0 )
  xcoor_mm_p    = base_x          + ( diff_x_pos >> 1 )
  xcoor_perfect = base_x_perfect  + ( diff_x_pos_perfect / 2. )
  
  # tan
  tanphi_x4096    = ( diff_x_pos * VERT_PHI1_PHI3_INV_Q15 ) >> (5+PRECISION)
  tanphi_perfect  = diff_x_pos_perfect / 235.
  
  # bx_time
  t         = (phi1['t'] + phi2['t'] ) >> 1
  t_perfect = (phi1['t_perfect'] + phi2['t_perfect'] ) /2.
  
  # path_calc_bx
  (bx_id, bx_id_p, bx_id_perfect) = path_calc_bx( t )
  bx_id_perfect = int(round(t_perfect/25.,0))
  
  # calc_primitive_chi_square
  phi3_offset_arr = [-sl_phi3_x_offset<<PRECISION, 0] if sl_phi3_x_offset<0 else [0, sl_phi3_x_offset<<PRECISION]
  
  valid =   [ segments[i/4]['valids'][i%4]                                                                                             for i in range(8) ]
  sum_A   = [ segments[i/4]['xcoor_mm_p'] - xcoor_mm_p + phi3_offset_arr[i/4] + segments[i/4]['delta_hit_x_p'][i%4] if valid[i] else 0 for i in range(8) ]
  sum_B   = [ Z_FACTOR_CORR[i]*tanphi_x4096                                                                         if valid[i] else 0 for i in range(8) ]
  factors = [ (sum_A[i]<<(14-PRECISION)) - sum_B[i]                                                                 if valid[i] else 0 for i in range(8) ]
  prods   = [ (factors[i]**2)>>2                                                                                    if valid[i] else 0 for i in range(8) ]
  chi2_mm2_p = ( sum( prods ) >> 16 ) # mult is 1024
  chi2_perfect = sum( [ ((segments[i/4]['xcoor_perfect'] - xcoor_perfect + phi3_offset_arr[i/4]/2**PRECISION + segments[i/4]['delta_hit_x_perfect'][i%4] - Z_FACTOR_CORR[i]*tanphi_perfect/4)**2) if valid[i] else 0 for i in range(8) ] ) / 100.
  
  if verbose:
    print 'Verbose values for path calculation of primitive'
    print '  phi2[xcoor_mm_p]', phi2['xcoor_mm_p']
    print '  phi3_offset_p', phi3_offset_p
    print '  phi1[xcoor_mm_p]', phi1['xcoor_mm_p']
    print '  diff_x_pos', diff_x_pos
    print '  old tan_phi_q_division', diff_x_pos * VERT_PHI1_PHI3_INV_Q15
    print '  tanphi_x4096', tanphi_x4096

  return {
    't'                                 : t,
    't_perfect'                         : t_perfect,
    'bx_id'                             : bx_id,
    'bx_id_perfect'                     : bx_id_perfect,
            
    'bx_id_chamber'                     : bx_id,                          # ab7manager format compatibility
    'bx_time_chamber'                   : t,                              # ab7manager format compatibility
              
    'chi2_mm2_p'                        : chi2_mm2_p,
    'chi2_cm2'                          : chi2_mm2_p/102400.,
    'chi2_perfect'                      : chi2_perfect,
    'chi_square_chamber'                : chi2_mm2_p,                     # ab7manager format compatibility
              
    'tanphi_x4096'                      : tanphi_x4096,
    'tanphi'                            : -tanphi_x4096/4096.,
    'tanphi_perfect'                    : -tanphi_perfect,
    'phi_tangent_chamber_float'         : tanphi_x4096/4096.,             # ab7manager format compatibility
    
    'xcoor_mm_p'                        : xcoor_mm_p,
    'xcoor'                             : xcoor_mm_p / 1. / 2**PRECISION,
    'xcoor_perfect'                     : xcoor_perfect,
    'horizontal_position_chamber_float' : xcoor_mm_p / 1. / 2**PRECISION, # ab7manager format compatibility
    
    'valid_segments'                    : 3,                              # ab7manager format compatibility
    'paired_segments'                   : segments,                       # ab7manager format compatibility
  }

def primitive_from_single_segment( segm, sl01 ) :
  prim = {  'bx_id_chamber'                     : segm['bx_id'],
            'bx_time_chamber'                   : segm['main_params']['bx_time'],
            'chi_square_chamber'                : segm['chi_square'],
            'phi_tangent_chamber_float'         : segm['phi_tangent_float'],
            'horizontal_position_chamber_float' : segm['horizontal_position_float'],
            'valid_segments'                    : 1<<sl01,
            'paired_segments'                   : [ NULL_SEGMENT, NULL_SEGMENT ],
            
            't'                                 : segm['t'],
            't_perfect'                         : segm['t_perfect'],
            'xcoor'                             : segm['xcoor'],
            'xcoor_perfect'                     : segm['xcoor_perfect'],
            'tanphi'                            : segm['tanphi'],
            'tanphi_perfect'                    : segm['tanphi_perfect'],
            'chi2_cm2'                          : segm['chi2_cm2'],
            'chi2_perfect'                      : segm['chi2_perfect'],
            }
  prim['paired_segments'][sl01] = segm
  return prim

def chamber_datapath( data, sl_phi3_x_offset ) :
  """ expects data = [ {'valid':list of bool, 't0s':list of int, 'wires':list of int} for sl in ['phi1','phi2'] ]"""
  segments = {}
  for sl in data:
    valid = data[sl]['valid']
    t0s   = data[sl]['t0s']
    wires = data[sl]['wires']
    if params['verbose_latcombs'] or params['verbose_path']:
      print '#'*40 + '\nCalculating segment for %s\n'%sl
    segments[sl] = superlayer_datapath( wires, t0s, valid )
    
    # if provided a lat combination, filter the outputs to segments that match
    # this laterality combination.
    # this is prepared for lat to be an array of 4 integers:
    # +1 right
    # 0 left
    # -1 not valid hit
    if 'lat' in data[sl]:
      segments[sl] = [
        seg for seg in segments[sl] if
        all( [not valid[i] or data[sl]['lat'][i]*2-1 == seg['lat'][i] for i in LAYERS] )
        ]

  tps = []
  if not all([segments[sl] for sl in segments]):
    tps = []
    for slname in segments:
      sl01 = 0 if slname == 'phi1' else 1
      tps += [ primitive_from_single_segment( segm, sl01 ) for segm in segments[slname] ]
  else:
    # for the moment just takes the first primitive from each SL and goes ahead with them
    segment_pair = [ segments[sl][0] for sl in ['phi1','phi2'] ]
    
    ## PATH CALCULATOR
    prim = primitive_path_calculator( segment_pair, sl_phi3_x_offset, params['verbose_path'] )

    tps = [ prim ]
    
  return tps
  

import arctanluts
def global_coordinates( primitive, wh, se, st ):
  luts = arctanluts.get(wh,se,st)

  # Depending on the type of primitive (SL1, SL3 or correlated), choose the
  # appropriate input data (x, tanspi) from the input primitive data structure
  # and the corresponding phi-lut from the 3 available options
  if primitive['valid_segments'] == 3 :
    x         = primitive['xcoor_mm_p']
    tanpsi    = primitive['tanphi_x4096']
    phi_lut   = luts['phic']
  elif primitive['valid_segments'] == 1:
    x         = primitive['paired_segments'][0]['xcoor_mm_p']
    tanpsi    = primitive['paired_segments'][0]['tanphi_x4096']
    phi_lut   = luts['phi1']
  else:
    x         = primitive['paired_segments'][1]['xcoor_mm_p']
    tanpsi    = primitive['paired_segments'][1]['tanphi_x4096']
    phi_lut   = luts['phi3']

  
  # x and slope are given in two's complement in fw
  x = to_two_comp(x, X_SIZE)
  tanpsi = to_two_comp(tanpsi, TANPSI_SIZE)

  ## Slice x and tanpsi
  # Both x and tanpsi are represented in vhdl as signed, this means their values
  # are stored as two's complement.
  
  # The MSB part is going to be used to index the luts and obtain a and b parameters
  # In vhdl, the indexing is automatic, I just extract the upper part of the signed
  # and pass it to the address port of the memory
  # In python, I'm converting it to an integer (with sign). The indexing of an array
  # in python with negative indexes results in indexing backwards from the end of the array
  # which, fortunately, returns the same value that the fw will pick from this lut
  x_msb = x >> (X_SIZE-PHI_LUT_ADDR_WIDTH)
  x_msb = from_two_comp(x_msb, PHI_LUT_ADDR_WIDTH)
  
  tanpsi_msb = tanpsi >> (TANPSI_SIZE-PHIB_LUT_ADDR_WIDTH)
  tanpsi_msb = from_two_comp(tanpsi_msb, PHIB_LUT_ADDR_WIDTH)
  
  # The LSB part can be sliced right away because it must yield a positive integer
  x_lsb = x & (2**(X_SIZE - PHI_LUT_ADDR_WIDTH)-1)
  tanpsi_lsb = tanpsi & (2**(TANPSI_SIZE - PHIB_LUT_ADDR_WIDTH)-1)

  ## Index the luts wiht the MSB parts already calculated
  phi_lut_q = phi_lut[to_two_comp(x_msb, PHI_LUT_ADDR_WIDTH)]
  phib_lut_q = luts['phib'][to_two_comp(tanpsi_msb, PHIB_LUT_ADDR_WIDTH)]
  
  ## Separate this data into the coefficients a and b
  ## and convert them to signed integers
  phi_lut_a = phi_lut_q['a&b'] >> PHI_LUT_B_BITS
  phi_lut_a = from_two_comp(phi_lut_a, PHI_LUT_A_BITS)
  
  phib_lut_a = phib_lut_q['a&b'] >> PHIB_LUT_B_BITS
  phib_lut_a = from_two_comp(phib_lut_a, PHIB_LUT_A_BITS)
  
  phi_lut_b = phi_lut_q['a&b'] & (2**PHI_LUT_B_BITS-1)
  phi_lut_b = from_two_comp(phi_lut_b, PHI_LUT_B_BITS)
  
  phib_lut_b = phib_lut_q['a&b'] & (2**PHIB_LUT_B_BITS-1)
  phib_lut_b = from_two_comp(phib_lut_b, PHIB_LUT_B_BITS)
  
  ## Do the linear piece-wise operations
  ## At this point all variables that can be negative have already been converted
  ## so will yield negative values when necessary
  phi_uncut = (phi_lut_b << PHI_B_SHL_BITS) + x_lsb * phi_lut_a 
  psi_uncut = (phib_lut_b << PHIB_B_SHL_BITS) + tanpsi_lsb * phib_lut_a 
  
  ## Trim phi to its final size
  ## This operation in python does what it has to do, the same as in fw:
  ## -5 >> 1 = -3
  phi = phi_uncut >> PHI_MULT_SHR_BITS

  ## Calculate phi_bending from the uncut version of phi and psi, and the trim it to size
  ## Same as before, be careful with the shift right operation
  phib_uncut = psi_uncut - (phi_uncut >>(PHI_PHIB_RES_DIFF_BITS + PHI_MULT_SHR_BITS - PHIB_MULT_SHR_BITS) )
  phib = phib_uncut >> PHIB_MULT_SHR_BITS

  primitive['phi']                = phi
  primitive['phi_float']          = primitive['phi'] * 1./2**17
  primitive['phi_bending']        = phib
  primitive['phi_bending_float']  = primitive['phi_bending'] * 4./2**13

  return primitive

#################
## RUN IF MAIN ##
#################

if __name__ == "__main__":
  from sl3offsets import sl3_ofssets
  from fwemudata import data
  (wh,se,st) = (data['wh'],data['se'],data['st'])
  sl_phi3_x_offset = sl3_ofssets[(wh,se,st)]

  input_data = {sl:data[sl] for sl in ['phi1','phi2']}
  tps = chamber_datapath(input_data, sl_phi3_x_offset)
  
  tps = [ global_coordinates( tpg, wh, se, st ) for tpg in tps ]

  for prim in tps:
    print '#'*50 

    print 'Laterality combinations: %s \n'%', '.join([str(prim['paired_segments'][i]['lat']) for i in range(2) if (prim['valid_segments']>>i)&1])

    print 'Primitive Time'
    print '  perfect', prim['t_perfect']
    print '  python ', prim['t']
    print '  fw     ', data['Time']['FW']
    print '  emul   ', data['Time']['Emul']

    print 'Position'
    print '  perfect', prim['xcoor_perfect']
    print '  python ', prim['xcoor']
    print '  fw     ', data['Pos']['FW']*10
    print '  emul   ', data['Pos']['Emul']*10

    print 'Tan'
    print '  perfect', prim['tanphi_perfect']
    print '  python  %.7f'%prim['tanphi']
    print '  fw     ', data['TanPsi']['FW']
    print '  emul   ', data['TanPsi']['Emul']

    print 'Chi2'
    print '  perfect %.9f'%prim['chi2_perfect']
    print '  python  %.9f'%prim['chi2_cm2']
    print '  fw      %.9f'%data['chi2']['FW']
    print '  emul    %.9f'%data['chi2']['Emul']

