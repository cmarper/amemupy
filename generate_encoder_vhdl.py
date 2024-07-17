#####################################
######## General Description ########
#####################################
## Optimization of the fw parameters for the fit of superlayer 3-hit or 4-hit segments
## We have equations to calculate directly the t0, pos and slope. Slope is actually
## calculated in two different units, the natural ones and also in units of mm/half-cell-height,
## which is convenient for the later calculation of the sum of residuals.
##
## A vhdl LUT is generated to store all these parameters.

###########################################
######## FW implementation details ########
###########################################
## For each magnitude, the result of the fit is calculated as follows
##     result = ( numer_const + sum(t_i * numer_coeff_i) ) / denom
## However, in fw, the division is implemented as a combination of a multiplication and a truncation
## (shift right). This is the number that will be actually stored in the lut for the division.
## Because the truncation is equivalent to a "floor" function, not round, a constant must be added to
## the numerator prior to the multiplication to ensure that it acts as a "round" function. This value
## is integrated into the numer_const number.
## The numer_const value has a decimal part, so the precision of the numerator will be increased
## before the final multiplication/truncation to be more precise.
## The multiplications factors in the numerator (numer_coeff_i) are small, so they will typically
## be implemented in logic, not DSP.


##########################################
##########################################
## IMPORTS, CONSTANTS AND CONFIGURATION ##
##########################################
##########################################

from pprint import pprint
import math

HALF_CELL_HEIGHT = 6.5
MAX_DRIFT_DIST = 21
MAX_DRIFT_TIME = 386.75
VDRIFT = MAX_DRIFT_DIST / MAX_DRIFT_TIME
LAYERS = range(4)

# Units are: time [ns], distance [mm], slope [mm/mm] but some of these fields
# are given with increasde bit resolutions to give more information
INCREASED_RES = {
  't0'        : 0,
  'pos'       : 4,
  'slope'     : 12,
  'slope_xhh' : 4,
  }

# When checking the division operation accuracy, this is the maximum accepted error in the LSB
# the error is considered actually symetrically w.r.t. 0
MAX_ERR_IN_LSB = { 
  't0'        : 1.2,
  'pos'       : 1.2,       
  'slope'     : 7, # slope result is given with a lot of resolution but not physical at the SL-level
  'slope_xhh' : 1.2, 
  }

# the list of magnitudes that are calculated with a direct operation
# in my order
mags = ['t0', 'pos', 'slope', 'slope_xhh']


########################
########################
## simplify constants ##
########################
########################
# Auxiliary function
# Takes a set of constants for one of the 4 direct fit calculations and scales it
# so that all the multiplicative coefficients on the input timestamps are the smallest
# possible integers

int_round = lambda x: int(round(x,0))

## used for simplification of coefficients
primes = [2,3,5,7,11,13,17,19]
prodprimes = 1
for p in primes: prodprimes*=p

def simplify_constants(numer_const, numer_coeff, denom):
  ## force all denominators to be positive
  if denom < 0:
    (numer_const, numer_coeff, denom) = (-numer_const, [-n for n in numer_coeff], -denom)

  # make all coefficients integer
  while any([ abs(num - int_round(num))>1e-5  for num in numer_coeff ]):
    numer_const *= prodprimes
    numer_coeff = [ numer_coeff[i]*prodprimes for i in LAYERS ]
    denom *= prodprimes
  numer_coeff = map(int_round, numer_coeff)
  
  ## simplify the coefficients and denominator if there is a common prime factor
  for factor in primes:
    while all([(num%factor) == 0 for num in numer_coeff]):
      (numer_const, numer_coeff, denom) = (numer_const/factor, [n/factor for n in numer_coeff], denom/factor)

  # round possible float values so that unique values are reached through different starting points
  numer_const = round(numer_const,10)
  denom = round(denom,10)

  #output
  return (numer_const, numer_coeff, denom)

##################################
##################################
## constants for T0, pos, slope ##
##################################
##################################
# The fit constants and coefficients

dotp = lambda x,y: sum([x[i]*y[i] for i in LAYERS])

def constants_from_layout_valid_latcomb(layout, valid, latcomb):
  if latcomb in [0,15]: return {key:(0, [0]*4, 0) for key in mags}

  xwire_ns  = [ layout[i]*MAX_DRIFT_TIME      *valid[i] for i in LAYERS ]
  xwire_mm  = [ layout[i]*MAX_DRIFT_DIST      *valid[i] for i in LAYERS ]
  z         = [ 2*(i-1.5)                     *valid[i] for i in LAYERS ]
  c         = [ (1 if (latcomb>>i)&1 else -1) *valid[i] for i in LAYERS ]
  b         = [ 1                             *valid[i] for i in LAYERS ]

  cz = (dotp(c,b)*dotp(z,b)-dotp(c,z)*dotp(b,b))/(dotp(z,z)*dotp(b,b)-dotp(z,b)*dotp(z,b))
  cb = (dotp(c,z)*dotp(b,z)-dotp(c,b)*dotp(z,z))/(dotp(z,z)*dotp(b,b)-dotp(z,b)*dotp(z,b))
  bc = (dotp(b,z)*dotp(c,z)-dotp(b,c)*dotp(z,z))/(dotp(c,c)*dotp(z,z)-dotp(c,z)*dotp(c,z))
  bz = (dotp(b,c)*dotp(z,c)-dotp(b,z)*dotp(c,c))/(dotp(c,c)*dotp(z,z)-dotp(c,z)*dotp(c,z))
  zb = (dotp(z,c)*dotp(b,c)-dotp(z,b)*dotp(c,c))/(dotp(b,b)*dotp(c,c)-dotp(b,c)*dotp(b,c))
  zc = (dotp(z,b)*dotp(c,b)-dotp(z,c)*dotp(b,b))/(dotp(b,b)*dotp(c,c)-dotp(b,c)*dotp(b,c))

  c_tilde = [ c[i] + cz * z[i] + cb * b[i] for i in LAYERS ]
  b_tilde = [ b[i] + bc * c[i] + bz * z[i] for i in LAYERS ]
  z_tilde = [ z[i] + zb * b[i] + zc * c[i] for i in LAYERS ]

  result = {}

  # T0
  numer_const = dotp(xwire_ns,c_tilde)
  numer_coeff = [ c[i]*c_tilde[i] for i in LAYERS ]
  denom = dotp(c_tilde,c_tilde)
  
  result['t0'] = simplify_constants(numer_const, numer_coeff, denom)
  
  # pos
  numer_const = dotp(xwire_ns,b_tilde)
  numer_coeff = [ c[i]*b_tilde[i] for i in LAYERS ]
  denom = dotp(b_tilde,b_tilde)/VDRIFT

  result['pos'] = simplify_constants(numer_const, numer_coeff, denom)
  
  # slope
  # note: only differs from the slope_xhh in the denominator
  numer_const = dotp(xwire_ns,z_tilde)
  numer_coeff = [ c[i]*z_tilde[i] for i in LAYERS ]
  denom = HALF_CELL_HEIGHT*dotp(z_tilde,z_tilde)/VDRIFT

  result['slope'] = simplify_constants(numer_const, numer_coeff, denom)

  # slope_xhh times half-cell-height
  numer_const = dotp(xwire_ns,z_tilde)
  numer_coeff = [ c[i]*z_tilde[i] for i in LAYERS ]
  denom = dotp(z_tilde,z_tilde)/VDRIFT

  result['slope_xhh'] = simplify_constants(numer_const, numer_coeff, denom)

  return result

########################################
########################################
## LATERALITY COMBINATIONS TO EXPLORE ##
########################################
########################################

LATCOMBS_FROM_LAYOUT_4h = {
# each laterality is an integer, with the msb corresponding to layer 4, and lsb corresponding to layer 1
# 0b0000 values are invalid, they are there just to complete to 6 values
  ( 0, -1, -2, -3 ) : [ 0b1000, 0b1100, 0b1110, 0b0001, 0b0011, 0b0111,                         ],
  ( 0, -1, -2, -1 ) : [ 0b0100, 0b0110, 0b0111,                         0b0000, 0b0000, 0b0000, ],
  ( 0, -1,  0, -1 ) : [ 0b0010, 0b1010, 0b0011, 0b1011,                         0b0000, 0b0000, ],
  ( 0, -1,  0,  1 ) : [ 0b0010, 0b0110, 0b1110,                         0b0000, 0b0000, 0b0000, ],
  ( 0,  1,  0, -1 ) : [ 0b0001, 0b1001, 0b1101,                         0b0000, 0b0000, 0b0000, ],
  ( 0,  1,  0,  1 ) : [ 0b0100, 0b1100, 0b0101, 0b1101,                         0b0000, 0b0000, ],
  ( 0,  1,  2,  1 ) : [ 0b1000, 0b1001, 0b1011,                         0b0000, 0b0000, 0b0000, ],
  ( 0,  1,  2,  3 ) : [ 0b1000, 0b1100, 0b1110, 0b0001, 0b0011, 0b0111,                         ], }

LATCOMBS_FROM_EXCLUDELYR_3h = {}
for exclude_layer in LAYERS:
  LATCOMBS_FROM_EXCLUDELYR_3h[exclude_layer] = []
  for lc in range(1,7):
    LATCOMBS_FROM_EXCLUDELYR_3h[exclude_layer] += [ (((lc&(8-2**exclude_layer)))<<1) + (lc&(2**exclude_layer-1)) ]

def latcombs_from_layout_valid(layout, valid):
  if all(valid): return LATCOMBS_FROM_LAYOUT_4h[layout]
  elif valid.count(1) == 3: return LATCOMBS_FROM_EXCLUDELYR_3h[valid.index(0)]
  else: return [0]*6

def latcombs_consts_encoder_from_layout_valid(layout, valid):
  result = []
  for latcomb in latcombs_from_layout_valid(layout, valid):
    constants = constants_from_layout_valid_latcomb(layout, valid, latcomb)
    result += [ (latcomb, constants) ]

  return result


#################################################
#################################################
## CHECK FW IMPLMENTATION DETAILS AND ACCURACY ##
#################################################
#################################################

# The size for the calculation of the quatifier of the corresponding divisions
# Values obtained through trial and error with the test that is done below 
# of the maximum acceptable error
DIV_SHR_BITS = {
  't0'        : 16,
  'pos'       : 21,
  'slope'     : 21,
  'slope_xhh' : 18,
  }

def mult_from_denom_mag(denom, mag):
  return int(round(2**DIV_SHR_BITS[mag]/1./denom,0))   if denom else 0

def calc_size(values):
  neg = 1 if any([v < 0 for v in values]) else 0
  maxconst = max([abs(v+0.5) for v in values]) # 0.5 to make sure that positive 2**N values require more size, but negative ones don't
  size = int(math.ceil(math.log(maxconst,2)))
  return size+neg, neg, maxconst

# Run through all the possible inputs to the LUT and get the parameters that would be
# stored in it, adding them to different data structures in order to get some useful
# information for the developement of the FW

constset = {key:set([]) for key in mags}
coeffset = {key:set([]) for key in mags}
denomset = {key:set([]) for key in mags}
multset  = {key:set([]) for key in mags}
miscset  = {key:set([]) for key in mags}
misc_from_denom = {key:{} for key in mags}
intervals_from_denom = {key:{} for key in mags}

for exclude_layer in [None] + LAYERS:
  valid = [0 if i== exclude_layer else 1 for i in LAYERS]
  for layout in LATCOMBS_FROM_LAYOUT_4h:
    for (latcomb, constants) in latcombs_consts_encoder_from_layout_valid( layout, valid ):
      for key in mags:
        ( numer_const, numer_coeffs, denom ) = constants[key]
        
        #simple sets with coefficients, denominators
        constset[key].add( numer_const )
        coeffset[key] = coeffset[key].union(numer_coeffs)
        denomset[key].add( denom )
        if denom:
          multset[key].add( mult_from_denom_mag(denom, key) )
        miscset[key].add( tuple([numer_const, denom]) )
        
        #storage of the coefficients corresponding for each denominator
        if denom:
          if denom not in misc_from_denom[key]: misc_from_denom[key][denom] = set([])
          misc_from_denom[key][denom].add( tuple(numer_coeffs+[numer_const]) )

        #storage of the maximum legal extremes for the interval of the numerator
        if denom:
          if denom not in intervals_from_denom[key]: intervals_from_denom[key][denom] = [1e5,-1e5]
          if key=='t0':
            current_interval = [
              # there is one hit with timestamp >= 512 ... minimum numerator that gives a valid drift:
              (512- MAX_DRIFT_TIME) * denom,
              # going for the maximum...
              numer_const + 1023 * sum([c for c in numer_coeffs if c>0]) + (1023-MAX_DRIFT_TIME)*sum([c for c in numer_coeffs if c<0])  ]
          else:
            # the extremes are calculated assuming worst case for input timestamps
            # but that fall within the MAX_DRIFT_TIME interval
            current_interval = [
              numer_const + MAX_DRIFT_TIME * sum([c for c in numer_coeffs if c<0])  ,
              numer_const + MAX_DRIFT_TIME * sum([c for c in numer_coeffs if c>0])  ]
          intervals_from_denom[key][denom] = [
            min( intervals_from_denom[key][denom][0], current_interval[0] )   ,
            max( intervals_from_denom[key][denom][1], current_interval[1] )   ]

#####################################################################################################################
#pprint(constset)
# From constset we can verify that the additive constant has to be stored with additional resolution
# in order to get perfect precision and the size of the necessary field:
# t0       -3094 to  1547, decimals in .5  steps --> signed 1 + 12 + 1
# pos     -13923 to 13923, decimals in .25 steps --> signed 1 + 14 + 2
# slope*   -1547 to  1547, decimals in .25 steps --> signed 1 + 11 + 2
#
# However...
# The numer_const, absorbs also the constant needed to make de division "round" instead of "floor"
# This constant is added after the numer_const is multiplied to get the additional precision, which
# for al magnitudes is bigger than the increased precision of the numer_const field, except for t0.
# For simplicity, for t0 we just disregard the .5 decimal parts, which does not cause a big
# rounding error.

## This constant is calculated from the previous print
NUMER_CONST_BASE_SIZE = {
  't0'        : 13,
  'pos'       : 15,
  'slope'     : 12,
  'slope_xhh' : 12,
  }
  
NUMER_CONST_SIZE = {
  mag:NUMER_CONST_BASE_SIZE[mag]+INCREASED_RES[mag]
  for mag in NUMER_CONST_BASE_SIZE }

# This is a printout and check that the previous hardcoded values are ok
if __name__ == '__main__':
  print '\nCalculate the bit size needed to accomodate the numerator additive constants...'
  for mag in mags:
    ( size, neg, maxvalue ) = calc_size(constset[mag])
    if size != NUMER_CONST_BASE_SIZE[mag]:
      raise Exception('Calculated size for %s (%i) does not match the number hardcoded in NUMER_CONST_BASE_SIZE'%(mag,size))
    print '  %s'%(mag+' '*(10-len(mag))) + '%i bits (%i for sign) + %i for additiona precision --> %i'%( size, neg, INCREASED_RES[mag], NUMER_CONST_SIZE[mag])


#####################################################################################################################
#pprint(coeffset)
# From coeffset we get the ranges needed to store each set of factors:
# t0       -2 to  7 --> signed 1 + 3
# pos     -14 to 14 --> signed 1 + 4
# slope*   -5 to  5 --> signed 1 + 3

NUMER_COEFF_SIZE = {
  't0'        : 4,
  'pos'       : 5,
  'slope'     : 4,
  'slope_xhh' : 4,
  }

# This is a printout and check that the previous hardcoded values are ok
if __name__ == '__main__':
  print '\nCalculate the bit size needed to accomodate the timestamp coefficients...'
  for mag in mags:
    ( size, neg, maxvalue ) = calc_size(coeffset[mag])
    if size != NUMER_COEFF_SIZE[mag]:
      raise Exception('Calculated size for %s (%i) does not match the number hardcoded in NUMER_COEFF_SIZE'%(mag,size))
    print '  %s'%(mag+' '*(10-len(mag))) + '%i bits (%i for sign)'%( size, neg )


#####################################################################################################################
#pprint(multset)
# From multset we get the ranges needed to store each set of multipliers
# the number of different possible numerators is very small, so ti may be more
# efficient to store an index and then convert the index to the multiplier in
# each case.
# The small factor for the Xilinx series 7 DSP is 18 bits (signed) so it's good
# if all of these are 17 bits or smaller.
# t0          --> 32768 --> unsigned 16
# pos         --> 56936 --> unsigned 16
# slope       -->  8759 --> unsigned 14
# slope_xhh   -->  7117 --> unsigned 13

MULT_SIZE = {
  't0'        : 16,
  'pos'       : 16,
  'slope'     : 14,
  'slope_xhh' : 13,
  }

# This is a printout and check that the previous hardcoded values are ok
if __name__ == '__main__':
  print '\nCalculate the bit size needed to accomodate the multiplier...'
  for mag in mags:
    ( size, neg, maxvalue ) = calc_size(multset[mag])
    if size != MULT_SIZE[mag]:
      raise Exception('Calculated size for %s (%i) does not match the number hardcoded in MULT_SIZE'%(mag,size))
    print '  %s'%(mag+' '*(10-len(mag))) + '%i --> %i bits (%i for sign)'%( maxvalue, size, neg )

#####################################################################################################################
# This code is to print the size needed to accomodate some intermediate steps

if __name__ == '__main__':
  print '\nCalculate the bit size needed to some steps of the datapath...'
  for mag in ['pos', 'slope', 'slope_xhh']:
    print '  %s'%mag

    # values for the numerator
    numerator_interval = [1e10,-1e10]
    for denom in intervals_from_denom[mag]:
      interval_denom = intervals_from_denom[mag][denom]
      numerator_interval = [ # this is used to estimate the fw's signal size for the numerator
        int(math.floor(min(interval_denom[0], numerator_interval[0]))) , 
        int(math.ceil(max(interval_denom[1], numerator_interval[1]))) ]
    ( size, neg, maxvalue ) = calc_size(numerator_interval)
    print ' '*12 + 'numerator              %i --> %i bits (%i for sign)'%( maxvalue, size, neg )
      
    # values for (numerator + prec) * mult
    mult_result_interval = [1e10,-1e10]
    for denom in intervals_from_denom[mag]:
      interval_denom = intervals_from_denom[mag][denom]
      mult = mult_from_denom_mag(denom, mag)
      mult_interval = [(i*2**INCREASED_RES[mag])*mult for i in interval_denom]
      mult_result_interval = [ # this is used to estimate the fw's signal size before 
        int(math.floor(min(mult_interval[0], mult_result_interval[0]))) , 
        int(math.ceil(max(mult_interval[1], mult_result_interval[1]))) ]
    ( size, neg, maxvalue ) = calc_size(mult_result_interval)
    print ' '*12 + '(numer*2**prec)*mult   %i --> %i bits (%i for sign)'%( maxvalue, size, neg )

    # values for calculated magnitude (after shr)
    print ' '*12 + 'result (after trunc)   %i --> %i bits (%i for sign)'%( maxvalue/2**DIV_SHR_BITS[mag], size-DIV_SHR_BITS[mag], neg )

    # only for slope_xhh, the product of slope_xhh * layer_h_in_hh
    if mag == 'slope_xhh':
      itvl = [ (value>> DIV_SHR_BITS[mag])*lyrhh for lyrhh in [-3,3] for value in mult_result_interval]
      ( size, neg, maxvalue ) = calc_size(itvl)
      print ' '*12 + 'slope_xhh*layer_h      %i --> %i bits (%i for sign)'%( maxvalue, size, neg )


##
## EVALUATION OF DIVISION PRECISION
##

intervals_from_geometry = {
# In fw, the timestamps have 10 bits, and it is known that at least one of the hits >= 512, So valid output values
# range from 125 to 1023
  't0':[125,1023],
# In fw, the maximum range for pos is [-42 mm, 42 mm]
  'pos': [ -42,42 ],
# In fw, the maximum range for slope is 2 cells over 3 layers (-> 2 x h) => 42 mm/cell_height
# for slope_xhh, result is given in mm/half_cell_height
  'slope_xhh': [ -21, 21 ],
# for slope, given in natural units
  'slope': [ -42./13, 42./13 ],
  }


for magnitude in mags:
  errors_hist = {}
  for denom in intervals_from_denom[magnitude]:
    # what interval to analyze
    interval_denom = intervals_from_denom[magnitude][denom]
    interval_geom  = intervals_from_geometry[magnitude]
    interval = [# take the most restrictive of the two
      int(math.floor(max(interval_denom[0], interval_geom[0]*denom))) , 
      int(math.ceil(min(interval_denom[1], interval_geom[1]*denom))) ] 

    mult = mult_from_denom_mag(denom, magnitude)
    for value in range(interval[0], interval[1]+1):
      postadd = 2**(DIV_SHR_BITS[magnitude]-1) # before the shift right, add half a unit, to ensure proper rounding
      preadd = int(round(postadd/1./mult,0)) # alternatively, this number can be added prior to the multiplication
      preadd = int(round(denom / 2.,0)) # this value is equivalent
      
      fwvalue = (((value<<INCREASED_RES[magnitude])+preadd)*mult ) >> (DIV_SHR_BITS[magnitude])
      error = fwvalue - value*(2**INCREASED_RES[magnitude])/1./denom
      error = round(error,2)
      if error not in errors_hist: errors_hist[error] = 0
      errors_hist[error] +=1
  
  if __name__ == '__main__':
    #print '\n%s --> '%magnitude, min(errors_hist.keys()), max(errors_hist.keys()),'\n'
    if max(errors_hist.keys()) > MAX_ERR_IN_LSB[magnitude]/2 or min(errors_hist.keys()) < -1*MAX_ERR_IN_LSB[magnitude]/2 :
      print '\nThe division width used for %s produces unacceptable error: '%magnitude, min(errors_hist.keys()), max(errors_hist.keys()),'\n'

##
## Drift time -> drift distance
##

DTDD_BITS = 13
mult = int(round(2**DTDD_BITS * VDRIFT ,0)) #445
postadd = 2**(DTDD_BITS-1)
preadd = int(round(postadd/1./mult,0)) #9

prodvalues = []
errors_hist = {}
for dt in range(387):
  prod_fw = ( ( (dt<<INCREASED_RES['pos']) + preadd ) * mult )
  dd_fw = prod_fw >> DTDD_BITS
  prodvalues += [prod_fw]
  error = dd_fw - dt*(2**INCREASED_RES['pos'])*VDRIFT
  error = round(error,2)
  if error not in errors_hist: errors_hist[error] = 0
  errors_hist[error] +=1

#pprint(errors_hist)  

if __name__ == '__main__':
  if max(errors_hist.keys()) > MAX_ERR_IN_LSB['pos']/2 or min(errors_hist.keys()) < -1*MAX_ERR_IN_LSB['pos']/2 :
    print '\nThe division width used for drift time to distance produces unacceptable error: ', min(errors_hist.keys()), max(errors_hist.keys()),'\n'

  ( size, neg, maxvalue ) = calc_size(prodvalues)

  print '\nCalculate the drift times to distances operation coefficients'
  print '  - The shift right number of bits is hardcoded to %i'%DTDD_BITS
  print '  - The pre-multiplication additive constant is %i'%preadd
  print '  - The multiplier is %i'%mult
  print '  - The product fits well in %i bits'%size
  print '    However, you\'ll want to keep it signed, so 1 additional bit'

##
## T0 --> bx_id
##

T0BX_BITS = 20
denom = 25.
mult = int(round(2**T0BX_BITS / 25.,0))
preadd = int(round(denom/2,0))

prodvalues = []
errors_hist = {}
for t0 in range(3564*25):
  prod_fw = ( ( t0 + preadd ) * mult )
  bxid_fw = prod_fw >> T0BX_BITS
  prodvalues += [prod_fw]
  error = bxid_fw - int(round(t0/25.,0))
  if error not in errors_hist: errors_hist[error] = 0
  errors_hist[error] +=1

#pprint(errors_hist)  

if __name__ == '__main__':
  if max(errors_hist.keys()) > 0 :
    print '\nThe division width used for t0 to bx_id produces unacceptable error'

  ( size, neg, maxvalue ) = calc_size(prodvalues)

  print '\nCalculate the T0 to bx_id operation coefficients'
  print '  - The shift right number of bits is hardcoded to %i'%T0BX_BITS
  print '  - The pre-multiplication additive constant is %i'%preadd
  print '  - The multiplier is %i'%mult
  print '  - The product fits well in %i bits'%size
  


########################
## Constants for VHDL ##
########################

def to_signed(value, size):
  if value > 2**(size-1)-1 or value < -2**(size-1):
    raise Exception("The value must be an integer in the range [-2**(size-1),2**(size-1)-1]")
  if value < 0:
    value = 2**size + value
  return bin((1<<size)+value)[3:]

def to_unsigned(value, size):
  if value > 2**size-1 or value < 0:
    raise Exception("The value (%i) must be an natural smaller than 2**%i-1"%(value,size))
  return bin((1<<size)+value)[3:]

def data_for_fw_from_constants( constants ):
  result = {}
  for mag in INCREASED_RES:
    (numer_const, numer_coeff, denom) = constants[mag]
    numer_const = to_signed(
      int(round( numer_const * (2**INCREASED_RES[mag]) + denom/2 ,0)),
      NUMER_CONST_SIZE[mag] )
    numer_coeff = [
      to_signed( coeff, NUMER_COEFF_SIZE[mag] )
      for coeff in numer_coeff ]
    mult = to_unsigned( mult_from_denom_mag(denom, mag), MULT_SIZE[mag] )
    result[mag] = {'coeff':numer_coeff, 'add':numer_const, 'mult': mult}
  #pprint(constants)
  #pprint(result)
  return result

#####################
## convert to VHDL ##
#####################

vhdl = ''
encoder_dict_real = {}
encoder_dict_slvs = {}

from_two_comp = lambda val,size:val - ((2*val)&(1<<size))

for exclude_layer in [None] + LAYERS:
  valid_vhdl = ''.join(['0' if i == exclude_layer else '1' for i in range(3,-1,-1)]) # the vhdl signal is (3 downto 0)
  valid = [0 if i== exclude_layer else 1 for i in LAYERS]

  for layout in LATCOMBS_FROM_LAYOUT_4h:
    vhdl += '    %sif\n'%('els' if vhdl else '')
    for i in LAYERS:
      vhdl += '      cell_horiz_layout(%i) = to_signed(% 2i,3) and\n'%(i,layout[i])
    vhdl += '      valid = "%s"\n'%(valid_vhdl)
    vhdl += '    then\n'
    
    encoded_list = latcombs_consts_encoder_from_layout_valid( layout, valid )
    encoder_dict_real[(layout, tuple(valid))] = encoded_list
    encoder_dict_slvs[(layout, tuple(valid))] = []

    for latid in range(6):
      (latcomb, constants) = encoded_list[latid]
      data_for_fw = data_for_fw_from_constants( constants )
      encoder_dict_slvs[(layout, tuple(valid))] += [(latcomb, data_for_fw)]

      constants = (   ''.join(data_for_fw['t0']['coeff'])
                    + data_for_fw['t0']['add']
                    + data_for_fw['t0']['mult']
                    + ''.join(data_for_fw['pos']['coeff'])
                    + data_for_fw['pos']['add']
                    + data_for_fw['pos']['mult']
                    + ''.join(data_for_fw['slope']['coeff'])
                    + data_for_fw['slope']['add']
                    + data_for_fw['slope']['mult']
                    + data_for_fw['slope_xhh']['add']
                    + data_for_fw['slope_xhh']['mult']
                  )
      
      latcomb = ''.join([str((latcomb>>i)&1) for i in range(3,-1,-1)]) # the vhdl signal is (3 downto 0)
      vhdl += '      lat_consts_arr(%i) <= ( "%s", "%s" );\n'%(latid, latcomb, constants )

def pylist_to_cpp(lista):
    mystr = "{"
    for i, elem in enumerate(lista):
        mystr += "{}".format(elem) + "," * (i != len(lista) - 1)
    mystr += "}"
    return mystr

if __name__ == '__main__':
  pass
  #print '#'*20 + '\nEncoder dict for python\n' + '#'*20 +'\n'
  #pprint(encoder_dict_slvs)
 
  #pprint(encoder_dict_real)

  cpp_slvs = ""

  for val in encoder_dict_slvs:
    cpp_slvs += "LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({\n"
    cpp_slvs +=  "  " + pylist_to_cpp((pylist_to_cpp(val[0]), pylist_to_cpp(val[1]))) + ",\n"
    cpp_slvs +=  "  {\n"
    for latcomb_consts in encoder_dict_slvs[val]:
        latcomb, constants = latcomb_consts
        mystr = "    {" + str(latcomb) + ", { "
        for mag in ["pos", "slope", "slope_xhh", "t0"]:
            mystr += "{"
            mystr += "{}, ".format(from_two_comp(int(constants[mag]['add'],2), NUMER_CONST_SIZE[mag]))
            mystr += "{}, ".format(pylist_to_cpp(   [from_two_comp(int(constants[mag]['coeff'][i],2),NUMER_COEFF_SIZE[mag]) for i in LAYERS]    ))
            mystr += "{}".format(int(constants[mag]["mult"], 2))
            mystr += "}, "
            # print from_two_comp(int(constants[mag]['add'],2), NUMER_CONST_SIZE[mag]),
        cpp_slvs += mystr + "}},\n"
                
    cpp_slvs += "  }});\n"
  
  cpp_real = ""
  for val in encoder_dict_real:
    cpp_real += "LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({\n"
    cpp_real +=  "  " + pylist_to_cpp((pylist_to_cpp(val[0]), pylist_to_cpp(val[1]))) + ",\n"
    cpp_real +=  "  {\n"
    for latcomb_consts in encoder_dict_real[val]:
        latcomb, constants = latcomb_consts
        mystr = "    {" + str(latcomb) + ", { "
        for mag in ["pos", "slope", "slope_xhh", "t0"]:
            mystr += "{"
            mystr += "{}, ".format(constants[mag][0])
            mystr += "{}, ".format(pylist_to_cpp(constants[mag][1]))
            mystr += "{}".format(constants[mag][2])
            mystr += "}, "
            # print from_two_comp(int(constants[mag]['add'],2), NUMER_CONST_SIZE[mag]),
        cpp_real += mystr + "}},\n"
                
    cpp_real += "  }});\n"
  
  
  # print '#'*20 + '\nEncoder vhdl\n' + '#'*20 +'\n'
  # print vhdl
  # print '#'*20 + '\nEncoder cpp\n' + '#'*20 +'\n'
  # print cpp_slvs
  # print cpp_real  
