#!python3
# Solve maths problem using SymPy

import sympy as sp
from sympy import cos, sin
from sympy.plotting import plot_parametric
from sympy.plotting.plot import MatplotlibBackend

sp.init_printing(use_latex='mathjax')

# %%
# Define the parametric functions
x, y, theta = sp.symbols('x y theta')
t = sp.symbols('t')

x = cos(3*t) / sin(2*t)
y = sin(3*t) / sin(2*t)


# %%
# Interval of integration
t1, t2 = sp.pi/6, sp.pi/3


# %%
# Compute dx/dt
dx_dt = sp.diff(x, t)


# %%
# Integrate y * dx/dt over [t1, t2] (before simplification)
area = sp.integrate(y * dx_dt, (t, t1, t2))

4 * area


# %%
# Simplify dx/dt
dx_dt = dx_dt.simplify()

# Integrate y * dx/dt over [t1, t2] (after simplification - THIS IS SLOW)
# sp.integrate(y * dx_dt, (t, t1, t2))

# %%
# Define the integrand
y_dx_dt = y * dx_dt

# Simplify the integrand
y_dx_dt_simplify = sp.simplify(y_dx_dt)

# Integrate y * dx/dt over [t1, t2] (after simplification - THIS IS SLOW)

sp.integrate(y_dx_dt_simplify, (t, t1, t2))



# %%
# Rewrite the Jacobian of the transformation base on dx_dt

jacobian = (-3*sin(3*t)*sin(2*t) - 2*cos(3*t)*cos(2*t))/sp.sin(2*t)**2


# Expand sin(3*t)*sin(2*t) using product-to-sum identity
# sin(A)sin(B) = (1/2)[cos(A-B) - cos(A+B)]
sin3t_sin2t_expanded = sp.Rational(1, 2) * (cos(3*t - 2*t) - cos(3*t + 2*t))
sin3t_sin2t_simplified = sp.simplify(sin3t_sin2t_expanded)


# cos(A)cos(B) = (1/2)[cos(A-B) + cos(A+B)]
# Expand cos(3*t)*cos(2*t) using product-to-sum identity
cos3t_cos2t_expanded = sp.Rational(1, 2) * (cos(3*t - 2*t) + cos(3*t + 2*t))
cos3t_cos2t_simplified = sp.simplify(cos3t_cos2t_expanded)

# Substitute expanded products into the Jacobian expression
jacobian_expanded = jacobian.subs({
        sin(3*t)*sin(2*t): sin3t_sin2t_expanded,
        cos(3*t)*cos(2*t): cos3t_cos2t_expanded,
    })

# Optional: simplify the result
jacobian_expanded_simplified = sp.simplify(jacobian_expanded)

# Define the integrand based on y * dx_dt

integrand = y * jacobian_expanded_simplified

integrand_simplified = sp.simplify(integrand)

integrand_simplified == y_dx_dt_simplify


# %%
# Separate integrand_simplified into additive terms (each as its own expression)
expanded_integrand = sp.expand(integrand_simplified, trig=False)
terms = [sp.simplify(t) for t in expanded_integrand.as_ordered_terms()]

# Expose individual terms as term_1, term_2, ...
for i, t_i in enumerate(terms, 1):
    globals()[f"term_{i}"] = t_i

terms

# %%

# This took around 5 minutes to compute
# integral_1 = sp.integrate(term_1.simplify(), (t, t1, t2))

from sympy import sqrt

integral_1 = -5 * sqrt(3) / 12

# This took around 17 minutes to compute
# integral_2 = sp.integrate(term_2.simplify(), (t, t1, t2))

integral_2 == -sqrt(3) / 12

integral_1 + integral_2

