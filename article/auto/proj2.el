(TeX-add-style-hook
 "proj2"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "aa"
    "aa10"
    "graphicx"
    "txfonts"
    "amsmath"
    "mathrsfs"
    "booktabs"
    "hyperref")
   (LaTeX-add-labels
    "eq:hydrostatic_equilibrium"
    "eq:mass_conservation"
    "eq:poly_relation"
    "eq:lane_emden"
    "eq:density"
    "eq:mass"
    "eq:pressure"
    "eq:temperature"
    "eq:luminosity"
    "sec:method"
    "fig:poly_ex"
    "eq:initial_cond_1"
    "eq:initial_cond_2"
    "eq:get_n"
    "sec:code_verification"
    "fig:sun_params"
    "fig:sun_emissivity"
    "fig:sun_luminosity"
    "eq:grav_analytical_energy"
    "eq:grav_poly_energy"
    "eq:grav_test"
    "sec:applications"
    "tab:data"
    "tab:indices"
    "fig:temperatures"
    "fig:densities"
    "fig:pressures"
    "fig:emissivities"
    "fig:emissivities_pp"
    "fig:emissivities_cno"
    "fig:luminosities"
    "fig:mc_alpha_a"
    "fig:mc_alpha_b"
    "sec:alternative_methods"
    "sec:regime"
    "sec:improvs")
   (LaTeX-add-bibliographies
    "poly"))
 :latex)

