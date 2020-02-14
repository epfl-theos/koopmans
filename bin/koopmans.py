#!/usr/bin/env python3

import argparse
from koopmans_cp import workflow

'''
Perform KI/KIPZ calculations
'''


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Perform a KI/KIPZ calculation using cp.x')
    parser.add_argument('json', metavar='system.json', type=str,
                        help='a single JSON file containing the workflow and cp.x settings')

    args = parser.parse_args()

    workflow.run(args.json)
