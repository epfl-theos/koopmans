import argparse
from koopmans import workflow

'''
Perform KI/KIPZ calculations
'''


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Perform a KI/KIPZ calculation using cp.x')
    parser.add_argument('orbdep', metavar='orbdep', type=str,
                        help="choice of orbital-dependent functional; must be one of 'ki'/'kipz'")
    parser.add_argument('json', metavar='system.json', type=str,
                        help='a JSON file containing the atomic positions, number of atoms etc.')
    parser.add_argument('-a', '--alpha', default=0.6, type=float,
                        help='starting guess for alpha as a single float')
    parser.add_argument('-c', '--cont', action='store_true',
                        help='continue from a previous calculation')
    parser.add_argument('-f', '--alpha_from_file', action='store_true',
                        help='read in starting guess for alpha from file')
    parser.add_argument('-i', '--maxit', default=1, type=int,
                        help='maximum number of self-consistent iterations')

    args = parser.parse_args()

    workflow.run(args.orbdep, args.json, args.alpha,
                 args.alpha_from_file, args.maxit, not args.cont)
