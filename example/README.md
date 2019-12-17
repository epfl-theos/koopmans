#### Example calculation for ki.py
To perform a calculation on water, run `$ python3 ../ki.py h2o-template.cpi -i 3`

The calculation should converge within two steps and you should obtain the following alpha values for each orbital: `0.667182  0.699768  0.699794  0.667177`

The final KI calculation is located in `final/ki.cpo`.

As a further test, you might like to try increasing `empty_states_nbnd`, or giving the calculation a different starting guess for alpha using `--alpha VALUE`.
