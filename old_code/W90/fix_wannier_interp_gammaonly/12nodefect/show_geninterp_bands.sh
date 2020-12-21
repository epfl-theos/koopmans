#!/bin/bash
cat Na_13chain_geninterp.dat  | grep -v '#' | awk '{print $1, $5;}' | ~/svn/precious/trunk/scripts/ordinabande.pl > .geninterp_deleteme && \
cat Na_13chain_geninterp_flat.dat  | grep -v '#' | awk '{print $1, $5;}' | ~/svn/precious/trunk/scripts/ordinabande.pl > .geninterp_flat_deleteme && \
xmgrace .geninterp_deleteme .geninterp_flat_deleteme