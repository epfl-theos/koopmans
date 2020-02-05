#### Example calculation for koopmans.py
To perform a calculation on water, run `$ python3 ../../koopmans.py h2o.json`

The calculation should converge within two steps and you should obtain the following alpha values for each orbital: `0.643232  0.709729  0.695975  0.695972`

The final KI calculation is located in `final/ki.cpo`.

As further tests, you might like to try altering some of the parameters in `h2o.json`. e.g.
 - increasing `empty_states_nbnd` (which will allow you to calculate electron affinities)
 - running KIPZ by setting both `calc_type` and `init_manifold` to `kipz`
