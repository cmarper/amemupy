## Code in this file emulates the lut generation, which in the system is done
## at configuration time. Calculating all the luts takes quite a while, so instead
## I'll be generating them as needed and storing them in a dict inside this module



import os, math

to_two_comp   = lambda val,size: val%(2**size)
from_two_comp = lambda val,size: val - ((2*val)&(1<<size))

## This function is copied from dtupy's ab7manager.py, **without modification**
def calc_atan_lut( msb_num, lsb_num, in_res, abscissa_0, out_res, a_extra_bits, b_extra_bits, a_size, b_size, sgn):
  """
  Calculates the coefficients needed to calculate the arc-tan function in fw
  by doing a piece-wise linear approximation.
  In fw, input (x) and output (y) are integers, these conversions are needed
    t = x*in_res - abscissa_0
    phi = arctan(t)
    y = phi/out_res
    => y = arctan(x*in_res - abcissa_0)/out_res
  The linear function is approximated as
    y = a*x_lsb + b
  Where a, b = func(x_msb) are the coefficients stored in the lut
  
  a is stored as unsigned, b as signed, with their respective sizes a_size, b_size,
  previously shifted left by a_extra_bits and b_extra_bits, respectively
  """
  
  a_min = -2**(a_size-1)
  a_max =  2**(a_size-1) - 1
  b_min = -2**(b_size-1)
  b_max =  2**(b_size-1) - 1
  
  tuck=False

  lut = {}
  
  for x_msb in range(-2**(msb_num-1),2**(msb_num-1)):
    # the values of x, t, phi, y at the extremes of the interval
    x1 = ( (x_msb    ) << lsb_num )
    x2 = ( (x_msb + 1) << lsb_num ) - 1

    t1 = x1*in_res - abscissa_0
    t2 = x2*in_res - abscissa_0

    phi1 = sgn*math.atan(t1)
    phi2 = sgn*math.atan(t2)
    
    y1 = phi1/out_res
    y2 = phi2/out_res
    
    # we want to find a, b so that the error in the extremes is the same as the error in the center
    # so the error in the extremes will be the same, so the "a" is determined by those points
    a = ( y2 - y1 ) / ( x2 - x1 )

    # by derivating the error function and equaling to 0, you get this is the point
    # towards the interval center with the highest error
    # err_f = y - (a*x+b) = sgn*arctan(x*in_res - abcissa_0)/out_res - (a*x+b)
    # d(err_f)/dx = sgn*(1/(1+(x*in_res - abcissa_0)^2))*in_res/out_res - a
    # d(err_f)/dx = 0 => x_max_err = (sqrt(in_res/out_res/a-1) + abscissa_0)/in_res
    # There is sign ambiguity in the sqrt operation. The sqrt has units of t (adimensional).
    # It is resolved by setting the sqrt to have the same sign as (t1+t2)/2
    t_max_err = math.sqrt(sgn*in_res/out_res/a - 1) * (1 if (t1+t2)>=0 else -1)
    x_max_err = (t_max_err + abscissa_0) / in_res
    phi_max_err = sgn*math.atan(t_max_err)
    y_max_err = phi_max_err / out_res

    # once you have the point of max error, the "b" parameter is chosen as the average between
    # those two numbers, which makes the error at the center be equal in absolute value
    # to the error in the extremes
    # units: rad
    b = (    y1    + y_max_err-a*(x_max_err-x1)    )/2
    
    # increase b in 1/2 of y_lsb, so that fw truncate operation on the of the result 
    # is equivalent to a round function instead of a floor function
    b += 0.5

    # shift left and round
    a = int(round(a*(2**a_extra_bits),0))
    b = int(round(b*(2**b_extra_bits),0))
    
    try:
      assert a >= a_min and a <= a_max and b >= b_min and b <= b_max
    except AssertionError as e:
      if not tuck and False:
        print "Tucked constant into bit size of output (un)signed integers. First time:"
        print x_msb
        print a, a_min, a_max
        print b, b_min, b_max
      tuck = True

    # tuck a, b constants into the bit size of the output (un)signed integer
    a = sorted([a_min, a, a_max])[1]
    b = sorted([b_min, b, b_max])[1]
    # convert a, b to two's complement
    a_signed = a % (2**a_size)
    b_signed = b % (2**b_size)
    
    lut[x_msb % 2**msb_num] = { # convert x_msb to two's complement signed
      'a':a,
      'b':b,
      'a&b': (a_signed<<b_size) + b_signed }
    
  return lut


## This class is defined so that the next function works with less modification
## to ease porting changes to amemupy.
class AB7():
  pass

## This code below is an almost exact copy of the same function from ab7manager
## self in the original is the object of a class function, but here it will 
## instead be a dict to which to store the luts
## params dict must contain keys wheel, sector, station
def load_atan_luts(self, params = None):
  if not params: params = self.getParams() # meaningless here, from orig. AB7
    
  if not 'slcoords' in AB7.__dict__:
    with open(os.path.expandvars('sl_phi0_point_perp_pos.txt'), 'r') as f:
      coords_txt = f.read()
    AB7.slcoords = {}
    for line in coords_txt.split('\n'):
      if not line.strip(): continue
      line = line.replace(':',' ').split()
      AB7.slcoords[tuple(map(int,[line[1],line[5],line[3],line[7]]))] = map(float,line[8:10])

  perps  = {sl:AB7.slcoords[(params['wheel'],params['sector'],params['station'],sl)][0] for sl in [1,3]}
  x_phi0 = {sl:AB7.slcoords[(params['wheel'],params['sector'],params['station'],sl)][1] for sl in [1,3]}
  
  luts = {}
  # calc_atan_lut( msb_num, lsb_num, in_res, abscissa_0, out_res, a_extra_bits, b_extra_bits, a_size, b_size, sgn):
  sgn = (-1 if (params['wheel']>0 or (params['wheel']==0 and params['sector'] in [2,3,6,7,10,11,14])) else 1)
  for sl in perps:
    luts['phi%i'%sl] = calc_atan_lut( 12, 6, (1./16)/(perps[sl]*10), x_phi0[sl]/perps[sl], 1./2**17, 10, 3, 12, 20, sgn)
  luts['phic'] =  calc_atan_lut( 12, 6, (1./16)/((perps[1]+perps[3])/.2), max(x_phi0[1],x_phi0[3])/((perps[1]+perps[3])/2), 1./2**17, 10, 3, 12, 20, sgn)
  luts['phib'] = calc_atan_lut( 9, 6, 1./4096, 0., 4./2**13, 10, 3, 10, 16, sgn)
  
  ## Original from ab7manager.py
  #for lut in luts:
  #  mem_data = [ luts[lut][x_msb]['a&b'] for x_msb in sorted(luts[lut].keys()) ]
  #  self.getNode( 'trigger.global_coordinates.lut_%s'%(lut)).writeBlock(mem_data)
  ## Replaced by:
  self[(params['wheel'],params['sector'],params['station'])] = luts


## This is the interface that this module offers to amemupy

luts_per_chamber = {}

def get(wh,se,st):
  global luts_per_chamber
  if (wh,se,st) not in luts_per_chamber:
    load_atan_luts(luts_per_chamber,  {'wheel':wh, 'sector':se, 'station':st} )
  
  return luts_per_chamber[(wh,se,st)]

def gen_all_luts():
  print 'Generating all luts...',
  for wh in range(-2,3):
    for se in range(1,15):
      for st in range(1,5) if se<13 else [4]:
        get(wh,se,st)
  print 'DONE!'