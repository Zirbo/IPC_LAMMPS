#! /usr/bin/python3


from math import cos, acos, pi

def eccentricity_from_delta(delta, gamma):
  dpm = (1.+delta)*.5
  c = cos(gamma*pi/180.)
  rho = (.25+dpm*(dpm-c))/(delta+1.-c)
  e = dpm - rho
  print("eccentricity: {}, patch radius: {}".format(e, rho))

def delta_from_eccentricity(e, rho):
  d = 2.*(e+rho-.5)
  g = 180.*acos(  (.25 + e*e - rho*rho)/e  )/pi
  print("delta: {}, gamma: {}".format(d, g))


def main():
  import argparse
  parser = argparse.ArgumentParser(description='computes eccentricity and patch radius from patch semiangle and interaction range, or viceversa')
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument('-gd', action='store_true', help='get eccentricity and radius from delta and gamma')
  group.add_argument('-er', action='store_true', help='get delta and gamma from eccentricity and radius')
# parser.add_argument('-m', '--mode', type=str, help='mode: gd to get eccentricity and radius from delta and gamma, er for the opposite ')
  parser.add_argument('-d', '--delta', type=float, help="interaction range delta, in units of hard core diameter", default=0)
  parser.add_argument('-g', '--gamma', type=float, help="patch semiamplitude gamma", default=0)
  parser.add_argument('-r', '--radius', type=float, help="radius", default=0)
  parser.add_argument('-e', '--eccentricity', type=float, help="eccentricity", default=0)

  args = parser.parse_args()
  if args.gd:
    eccentricity_from_delta(args.delta, args.gamma)
  elif args.er:
    delta_from_eccentricity(args.eccentricity, args.radius)
  else:
    parser.print_help()


if __name__ == "__main__":
  main()
