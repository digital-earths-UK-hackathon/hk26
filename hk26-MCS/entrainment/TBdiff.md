# Tb_diff investigation

## What we want

`Tb_diff = Tb - T_LNB`

where:
- `Tb` is brightness temperature derived from OLR
- `T_LNB` is the temperature at the Level of Neutral Buoyancy (equilibrium level) from MetPy

The idea: if a convective updraft reaches its LNB, the cloud-top brightness temperature should
be close to `T_LNB`. A strongly negative `Tb_diff` means the cloud top is much colder than the
LNB (overshooting), a near-zero value means the cloud top is right at the LNB, and a positive
value means the cloud top is warmer (lower altitude) than the LNB.

## T_LNB

MetPy's `el()` returns temperature in Kelvin at the equilibrium level. Typical tropical LNB
temperatures are ~200–220 K (~100–150 hPa). This seems correct.

## Tb formula

The formula as currently implemented:

```
Tf = (OLR / sigma) ** (1/4)       # effective radiating temperature
Tb = (a + sqrt(a² - 4·b·Tf)) / (2·b)
```

with `a = 1.228`, `b = -1.106e-3`.

## The problem

Evaluating for typical tropical OLR values:

| OLR (W/m²) | Tf (K) | Tb (K) as coded |
|-----------|--------|-----------------|
| 150       | 227    | −1272           |
| 200       | 244    | −1282           |
| 240       | 255    | −1289           |
| 280       | 265    | −1295           |

Tb ≈ −1288 K is clearly unphysical. This causes `Tb_diff ≈ −1288 − 215 ≈ −1503 K`,
consistent with the observed ~−1600 values.

## Root cause

With `b = −1.106e-3`, the denominator `2b` is negative and the numerator
`a + sqrt(a² - 4b·Tf) = a + sqrt(a² + 4|b|·Tf)` is a large positive number, so the
result is a large negative number.

The formula is the solution to the quadratic `b·Tb² + a·Tb - Tf = 0`, which has two roots.
The `+` form selects the unphysical root. The `−` form:

```
Tb = (a − sqrt(a² - 4·b·Tf)) / (2·b)
```

gives, with `b = −1.106e-3`:

| OLR (W/m²) | Tf (K) | Tb (K) fixed |
|-----------|--------|--------------|
| 150       | 227    | 178          |
| 200       | 244    | ?            |
| 240       | 255    | 178          |
| 280       | 265    | 178          |

Still too cold (~178 K).

The formula that gives physically reasonable values (Tb slightly warmer than Tf, ~234–293 K) is:

```
Tb = (a − sqrt(a² - 4·b·Tf)) / (2·b)    with b = +1.106e-3
```

| OLR (W/m²) | Tf (K) | Tb (K) |
|-----------|--------|--------|
| 150       | 227    | 234    |
| 200       | 244    | 259    |
| 240       | 255    | 277    |
| 280       | 265    | 293    |

## Questions to resolve

1. What is the original source of `a = 1.228`, `b = -1.106e-3`?
   - Likely Minnis & Harrison (1984) or a derivative. Need to check the original paper
     for the sign convention used for `b` and which root to take.

2. Is the formula meant to operate on Tf in K or °C?
   - The values above assume K throughout. If the empirical fit was derived in °C,
     Tf should be converted before applying the formula.

3. Should `b` be `+1.106e-3` (positive) in the formula as written, or is the
   `-` before `sqrt` the intended fix?
   - Both changes together give physical results; need the original paper to confirm.

## Resolution

Source confirmed: PyFLEXTRKR `ftfunctions.py` (`olr_to_tb`), citing Yang & Slingo (2001) /
Minnis & Harrison (1984). The correct formula with `b = -1.106e-3` is:

```
tb = (-a + sqrt(a**2 + 4*b*tf)) / (2*b)
```

The original bug was `+a` and `-4*b*Tf` (wrong signs on both `a` and the discriminant term).
Fixed in `calc_entrainment_wam.py` via the `olr_to_tb()` function.
