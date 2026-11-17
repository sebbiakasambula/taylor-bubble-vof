# Set 17, part 3. Film velocity at T1 against Brown's profile

G = 999.4999, window from 2.0 D below the nose to 0.5 D above the axis tail, snapshots in 5 to 10 s.

| run | snapshots | columns | film flux / theory, median (10th to 90th) | interface-cell speed / theory, median (10th to 90th) | mean film in these columns | columns by distance behind the nose (D) |
|---|---|---|---|---|---|---|
| OF baseline | 8 | 109 | 0.7629 (0.7128 to 0.8088) | 0.7431 (0.6884 to 0.7915) | 0.0655 | 2.01 to 2.76 |
| OF predictor | 21 | 345 | 0.7734 (0.7190 to 0.8166) | 0.7539 (0.6956 to 0.7997) | 0.0653 | 2.00 to 2.86 |
| OF geometric and predictor | 21 | 397 | 0.8036 (0.7459 to 0.8644) | 0.7861 (0.7239 to 0.8533) | 0.0637 | 2.00 to 2.98 |
| Fluent | 8 | 181 | 0.7864 (0.7064 to 0.8408) | 0.7673 (0.6818 to 0.8254) | 0.0656 | 2.00 to 3.39 |
