====================
H2 molecule with EMT
====================


Google AI
~~~~~~~~~
can ASE compute anharmonicity constant of diatomic molecule ?
give example for H2 with EMT
add comparison with experimental values for all computed data
write the full code

several round of fixing

run
~~~
python code01.py
==================================================================
Spectroscopic Property    | Computed   | Exp.     | Error (%)
==================================================================
Bond Length (r_e)         | 0.780      | 0.741    | +5.3    %
Harmonic Freq. (omega_e)  | 6922.0     | 4401.0   | +57.3   %
Anharmonicity (omega_e_xe) | 516.7      | 121.3    | +326.0  %
==================================================================

To bring these values within 1-2% of the real-world experimental markers, you must substitute the semi-empirical EMT calculator with a true quantum mechanical electronic structure method.






