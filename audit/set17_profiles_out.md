# Set 17, part 2. Film velocity at T2 against Brown's profile

G = drho g / mu_l = 49.9750. Window: 1.5 D below the nose to 0.5 D above the axis tail.

| run | snapshots | columns | film flux / theory, median (10th to 90th percentile) | interface-cell speed / theory, median (10th to 90th) | mean film in these columns |
|---|---|---|---|---|---|
| OF baseline | 5.823776, 6.655744, 7.487712, 8.31968 | 120 | 0.9785 (0.9767 to 0.9789) | 0.9735 (0.9717 to 0.9739) | 0.1340 |
| OF predictor | 5.2, 6, 7, 8.2 | 123 | 0.9748 (0.9712 to 0.9787) | 0.9679 (0.9640 to 0.9724) | 0.1366 |
| OF geometric and predictor | 5.2, 6, 7, 8.2 | 129 | 0.9734 (0.9703 to 0.9748) | 0.9663 (0.9627 to 0.9679) | 0.1409 |
| Fluent | steps 6160, 7040, 7920, 8800, 9680 | 167 | 0.9733 (0.9709 to 0.9745) | 0.9661 (0.9632 to 0.9672) | 0.1420 |

Middle column of the last snapshot of each run (figure):

- OF baseline at t 8.3197: column 2.0896 D behind the nose, film 0.1340, 41 liquid cells, interface speed 0.3948 against theory 0.4054
- OF predictor at t 8.2000: column 2.1182 D behind the nose, film 0.1365, 42 liquid cells, interface speed 0.4070 against theory 0.4201
- OF geometric and predictor at t 8.2000: column 2.1534 D behind the nose, film 0.1408, 43 liquid cells, interface speed 0.4308 against theory 0.4453
- Fluent at t 8.2000: column 2.1870 D behind the nose, film 0.1422, 43 liquid cells, interface speed 0.4380 against theory 0.4531
