# SHE Angle Solver for MATLAB
A nonlinear system solver for Selective Harmonic Elimination (SHE) in a 6-pulse Current Source Inverter (CSI).
Computes 6 switching angles that eliminate harmonics 5, 7, 11, 13, and 17 while controlling the fundamental component.

---

## How It Works

The function defines a system of 6 nonlinear equations solved externally using `fsolve`:

1. **Fundamental Control** - Sets the fundamental harmonic amplitude equal to a target value `I_1_target`.
2. **Harmonic Elimination** - Forces the 5th, 7th, 11th, 13th, and 17th harmonic amplitudes to zero.
3. **Fourier Coefficients** - Each harmonic amplitude is computed using the standard CSI half-wave symmetric Fourier series via the helper function `compute_a_n`.

---

## Usage

```matlab
% Define parameters
I_dc = 200;           % DC link current (A)
I_1_target = 1.0;     % Target fundamental amplitude (normalized or absolute)
x0 = [10 20 30 40 50 60] * (pi/180);  % Initial angle guess (rad)

% Solve
options = optimoptions('fsolve', 'TolFun', 1e-12, 'MaxFunEvals', 999999);
alpha_sol = fsolve(@(x) SHE_nonlinear_system_with_jacobian(x, I_dc, I_1_target), x0, options);

% Convert to degrees
alpha_deg = sort(alpha_sol * 180/pi)
```

| Parameter     | Description |
|---------------|-------------|
| `alpha`       | Vector of 6 switching angles (rad) |
| `I_dc`        | DC link current magnitude |
| `I_1_target`  | Desired fundamental harmonic amplitude |
| `F`           | Residual vector - zero at the SHE solution |

---

## Tuning

| Parameter   | Default | Effect |
|-------------|---------|--------|
| `x0`        | User defined | Initial angle guess - affects convergence, try values spread between 0 and pi/2 |
| `TolFun`    | `1e-12` | Solver tolerance - tighten for higher accuracy |
| Harmonics   | 5, 7, 11, 13, 17 | Modify `F` equations to target different harmonic orders |

---

## Limitations

- Assumes **half-wave symmetry** in the CSI switching pattern
- Convergence depends heavily on the initial guess `x0`
- Valid switching angles must lie in `[0, pi/2]` - solutions outside this range are physically invalid
- Does not include a Jacobian explicitly - `fsolve` uses numerical differentiation by default

---

## Author

Aman Jaswal
