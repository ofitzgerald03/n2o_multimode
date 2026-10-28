"""
Appendix B specification: statement, tolerance, layer and test level of every
check. Entries are (label prefix, sentence, tolerance, layer), in log order.
"""

KG = r"\SI{0.02}{\kilo\gram}"
ARCH = {
    "single_elec": "single-electric",
    "multimode": "bipropellant multimode",
    "hybrid": "hybrid",
    "single_chem": "single-chemical",
    "mono_multimode": "monopropellant multimode",
}
ORDER5 = ("single_elec", "multimode", "hybrid", "single_chem", "mono_multimode")

SL = []
# LEO mission results
for k in ORDER5:
    SL.append(
        (
            f"{k} wet/dry/prop/tanks/eps",
            f"The {ARCH[k]} design of the LEO mission reproduces its frozen wet, dry, loaded-propellant, tank and power-system masses.",
            KG + " on each",
            "system",
        )
    )
    SL.append(
        (
            f"{k} converged",
            f"The {ARCH[k]} design of the LEO mission converges.",
            "exact",
            "system",
        )
    )
# differences between architectures
SL += [
    (
        "mm - se delta",
        "The bipropellant multimode minus single-electric wet-mass difference in the LEO mission reproduces its frozen value.",
        KG,
        "system",
    ),
    (
        "hy - mm delta",
        "The hybrid minus bipropellant multimode wet-mass difference in the LEO mission reproduces its frozen value.",
        KG,
        "system",
    ),
    (
        "mono_mm - mm delta",
        "The monopropellant minus bipropellant multimode wet-mass difference in the LEO mission reproduces its frozen value.",
        KG,
        "system",
    ),
]
for k in ORDER5:
    n = "no electric unit" if k == "single_chem" else "one electric unit"
    SL.append(
        (
            f"{k} n_ep ==",
            f"The {ARCH[k]} design of the LEO mission installs {n}.",
            "exact",
            "system",
        )
    )
# tank sizing, LEO
for k in ("multimode", "hybrid", "single_chem"):
    SL.append(
        (
            f"COPV {k} wet",
            f"With a composite \\ce{{N2O}} tank, the {ARCH[k]} design of the LEO mission reproduces its frozen wet mass.",
            KG,
            "system",
        )
    )
SL.append(
    (
        "COPV ordering",
        "With a composite \\ce{N2O} tank, the bipropellant multimode design of the LEO mission stays lighter than the hybrid.",
        "exact (sign)",
        "system",
    )
)
# power of the two electric units
SL.append(
    (
        "T_req",
        "The thrust needed to fly the GEO transfer within its window, computed in closed form at a fixed spacecraft mass, reproduces its frozen value.",
        r"\SI{0.1}{\milli\newton}",
        "system",
    )
)
for u, name in (("hall", "Hall"), ("ecr", "electrodeless")):
    SL += [
        (
            f"{u} n_units ==",
            f"The closed-form transfer check installs the frozen number of {name} units.",
            "exact",
            "system",
        ),
        (
            f"{u} P_dc",
            f"The bus power of that {name} fleet reproduces its frozen value.",
            r"\SI{0.01}{\kilo\watt}",
            "system",
        ),
        (
            f"{u} EPS increment",
            f"The added power-system mass of that {name} fleet reproduces its frozen value.",
            r"\SI{0.1}{\kilo\gram}",
            "system",
        ),
        (
            f"{u} electronics",
            f"The PPU mass of that {name} fleet reproduces its frozen value.",
            r"\SI{0.1}{\kilo\gram}",
            "system",
        ),
    ]
SL += [
    (
        "power-system limit separates",
        "In the closed-form transfer check the Hall fleet lies within the power-system limit and the electrodeless fleet does not.",
        "exact",
        "system",
    ),
    (
        "chem transfer propellant",
        "The propellant of a chemical transfer in the closed-form check reproduces its frozen value.",
        r"\SI{0.5}{\kilo\gram}",
        "system",
    ),
]
# allocation limits
SL += [
    (
        "f_EP = 0 offset within",
        "With no $\\Delta V$ assigned to its electric mode, the bipropellant multimode design of the LEO mission exceeds the single-chemical design by its idle electrodeless hardware.",
        r"\SI{20}{\percent} of that hardware",
        "system",
    ),
    (
        "expected f0 offset",
        "The idle electrodeless hardware (head and rated PPU with the structure margin) reproduces its frozen mass.",
        r"\SI{0.01}{\kilo\gram}",
        "system",
    ),
    (
        "expected f1 offset",
        "The idle bipropellant hardware (head and empty propene tank with the structure margin) reproduces its frozen mass.",
        r"\SI{0.01}{\kilo\gram}",
        "system",
    ),
    (
        "f_EP = 1 gap decomposes",
        "With all $\\Delta V$ assigned to its electric mode, the difference between the bipropellant multimode and single-electric designs of the LEO mission equals the sum of its named terms.",
        r"residual within \SI{20}{\percent} of the idle bipropellant hardware",
        "system",
    ),
]
# allocation sweep
SL += [
    (
        "EP loss-inclusive dv non-decreasing",
        "Over the electric-share grid of the LEO mission, the electric $\\Delta V$ never falls and the chemical $\\Delta V$ never rises.",
        r"\num{1e-9} in $\Delta V$",
        "sweep",
    ),
    (
        "all f_EP grid closures converge",
        "Every point of the electric-share grid of the LEO mission converges.",
        "exact",
        "sweep",
    ),
]
# GEO transfer mission results
for k in ORDER5:
    SL += [
        (
            f"GEO {k} wet/dry/prop/tanks/eps",
            f"The {ARCH[k]} design of the GEO transfer mission reproduces its frozen wet, dry, loaded-propellant, tank and power-system masses.",
            KG + " on each",
            "system",
        ),
        (
            f"GEO {k} converged",
            f"The {ARCH[k]} design of the GEO transfer mission converges.",
            "exact",
            "system",
        ),
        (
            f"GEO {k} n_ep ==",
            f"The {ARCH[k]} design of the GEO transfer mission installs its frozen number of electric units.",
            "exact",
            "system",
        ),
    ]
SL.append(
    (
        "GEO mono_multimode closes",
        "The monopropellant multimode design of the GEO transfer mission converges and lies above the envelope limit.",
        "exact",
        "system",
    )
)
for k in ("single_elec", "multimode", "hybrid", "single_chem"):
    SL.append(
        (
            f"GEO {k} nominal excluded",
            f"The {ARCH[k]} design of the GEO transfer mission passes every schedule, duty and power check and lies above the envelope limit.",
            "exact",
            "system",
        )
    )
SL.append(
    (
        "GEO multimode all-EP transfer",
        "A fully electric transfer by the bipropellant multimode architecture has no fixed point: the closure is flagged as diverged and returns no mass.",
        "exact",
        "system",
    )
)
# collision-avoidance check
for h in (2, 6, 24, 72):
    SL.append(
        (
            f"chem CA feasible at t_warn = {h} h",
            f"A chemical avoidance burn completes within a warning time of \\SI{{{h}}}{{\\hour}}.",
            "exact",
            "system",
        )
    )
SL += [
    (
        "single_elec CA verdicts",
        "The single-electric avoidance verdict over warning times of \\SIlist{2;6;24;72}{\\hour} is infeasible at \\SI{2}{\\hour} and feasible at the other three.",
        "exact",
        "sweep",
    ),
    (
        "1 d centroid-charge factor",
        "At a warning time of one day the finite-burn correction of an electric avoidance burn reproduces its frozen factor.",
        r"\num{0.005}",
        "system",
    ),
]
# allocation sweep structure
SL += [
    (
        "all sweep closures converge",
        "In both missions every point of the electric-share sweep either converges or is flagged as diverged; none ends without a result.",
        "exact",
        "sweep",
    ),
    (
        "the divergences on the grid are exactly",
        "The points of the electric-share sweep that diverge are exactly the frozen set, all in the GEO transfer mission.",
        "exact",
        "sweep",
    ),
    (
        "GEO mono_multimode f_EP divergence set",
        "The monopropellant multimode architecture diverges in the GEO transfer mission at exactly the frozen electric shares.",
        "exact",
        "sweep",
    ),
    (
        "GEO allocation monotone",
        "Over the electric-share grid of the GEO transfer mission, the electric $\\Delta V$ never falls and the chemical $\\Delta V$ never rises.",
        r"\num{1e-9} in $\Delta V$",
        "sweep",
    ),
    (
        "GEO f_EP = 0 gap decomposes",
        "With no $\\Delta V$ assigned to its electric mode, the difference between the bipropellant multimode and single-chemical designs of the GEO transfer mission equals the sum of its named terms.",
        r"residual within \SI{20}{\percent} of the idle electrodeless hardware",
        "system",
    ),
    (
        "GEO f_EP = 1 gap decomposes",
        "With all $\\Delta V$ assigned to its electric mode, the difference between the bipropellant multimode and single-electric designs of the GEO transfer mission equals the sum of its named terms.",
        r"residual within \SI{20}{\percent} of the idle bipropellant hardware",
        "system",
    ),
]
# interface to the thruster-level model
SL.append(
    (
        "ecr_unit(0.61) == ecr_unit()",
        "The electrodeless unit rebuilt at the measured xenon utilisation returns the nominal specific impulse, thrust, mass flow rate, bus power and head mass.",
        "exact",
        "thruster",
    )
)
for eu in ("0.2", "0.85"):
    SL += [
        (
            f"interface Isp at eta_u = {eu}",
            f"At a utilisation of \\num{{{eu}}} the specific impulse of the electrodeless unit reproduces its independently computed value.",
            r"relative \num{1e-9}",
            "thruster",
        ),
        (
            f"interface mdot at eta_u = {eu}",
            f"At a utilisation of \\num{{{eu}}} the mass flow rate of the electrodeless unit reproduces its independently computed value.",
            r"relative \num{1e-9}",
            "thruster",
        ),
        (
            f"T and P_dc held at eta_u = {eu}",
            f"At a utilisation of \\num{{{eu}}} the thrust and bus power of the electrodeless unit are unchanged from their nominal values.",
            "exact",
            "thruster",
        ),
    ]
SL += [
    (
        "raw pipeline parity op_point",
        "The thrust, specific impulse, mass flow rate and absorbed power computed by the thruster-level model at parity agree with the input values of the system-level model.",
        r"relative \num{2e-3}",
        "thruster",
    ),
    (
        "system EP Isp == chemical Isp",
        "The utilisation at which the electric specific impulse equals the bipropellant specific impulse reproduces its frozen value, and the unit rebuilt there has that specific impulse.",
        r"\num{1e-5} in utilisation; \SI{1e-6}{\second}",
        "thruster",
    ),
]
# utilisation sweep
SL.append(
    (
        "LEO eta_u grid: 41/41",
        "Every point of the utilisation sweep of the LEO mission converges and none is excluded.",
        "exact",
        "sweep",
    )
)
for eu in ("0.05", "0.61", "0.85"):
    SL.append(
        (
            f"leo sweep wet at eta_u = {eu}",
            f"The bipropellant multimode wet mass of the LEO mission at a utilisation of \\num{{{eu}}} reproduces its frozen value.",
            KG,
            "sweep",
        )
    )
SL += [
    (
        "leo f_EP = 0 reference",
        "The own-chemical reference mass of the bipropellant multimode architecture in the LEO mission reproduces its frozen value.",
        KG,
        "sweep",
    ),
    (
        "leo minimum viable eta_u (own_chem)",
        "The utilisation threshold of the bipropellant multimode architecture against its own chemical mode in the LEO mission reproduces its frozen value.",
        r"\num{1e-3}",
        "sweep",
    ),
    (
        "leo minimum viable eta_u (single_chem)",
        "The utilisation threshold of the bipropellant multimode architecture against the single-chemical architecture in the LEO mission reproduces its frozen value.",
        r"\num{1e-3}",
        "sweep",
    ),
    (
        "leo min-viable bracket",
        "In the LEO mission the own-chemical threshold separates the two sides: just below it the multimode design is heavier than the reference, just above it no heavier.",
        r"$\pm \num{0.005}$ in utilisation",
        "sweep",
    ),
    (
        "GEO eta_u grid: 41/41",
        "Every point of the utilisation sweep of the GEO transfer mission converges, and the point at a utilisation of \\num{0.05} reproduces its frozen mass.",
        r"\SI{0.5}{\kilo\gram}",
        "sweep",
    ),
    (
        "GEO eta_u grid: every point excluded",
        "Every point of the utilisation sweep of the GEO transfer mission lies above the envelope limit.",
        "exact",
        "sweep",
    ),
]
for eu in ("0.07", "0.61", "0.85"):
    SL.append(
        (
            f"geo sweep wet at eta_u = {eu}",
            f"The bipropellant multimode wet mass of the GEO transfer mission at a utilisation of \\num{{{eu}}} reproduces its frozen value.",
            r"\SI{0.5}{\kilo\gram}",
            "sweep",
        )
    )
SL += [
    (
        "geo f_EP = 0 reference",
        "The own-chemical reference mass of the bipropellant multimode architecture in the GEO transfer mission reproduces its frozen value.",
        KG,
        "sweep",
    ),
    (
        "geo minimum viable eta_u (own_chem)",
        "The utilisation threshold of the bipropellant multimode architecture against its own chemical mode in the GEO transfer mission reproduces its frozen value.",
        r"\num{1e-3}",
        "sweep",
    ),
    (
        "geo minimum viable eta_u (single_chem)",
        "The utilisation threshold of the bipropellant multimode architecture against the single-chemical architecture in the GEO transfer mission reproduces its frozen value.",
        r"\num{1e-3}",
        "sweep",
    ),
    (
        "geo min-viable bracket",
        "In the GEO transfer mission the own-chemical threshold separates the two sides of the reference mass.",
        r"$\pm \num{0.005}$ in utilisation",
        "sweep",
    ),
    (
        "LEO mono_multimode grid: 41/41",
        "Every point of the monopropellant multimode utilisation sweep of the LEO mission converges and none is excluded.",
        "exact",
        "sweep",
    ),
    (
        "LEO mono_multimode own-chem",
        "The own-chemical reference of the monopropellant multimode architecture in the LEO mission, flown at the monopropellant specific impulse, reproduces its frozen mass.",
        KG,
        "sweep",
    ),
    (
        "LEO mono_multimode parity grid point",
        "At the measured xenon utilisation the monopropellant multimode sweep reproduces the nominal LEO design.",
        KG,
        "sweep",
    ),
    (
        "LEO mono_multimode minimum viable eta_u (own_chem)",
        "The monopropellant multimode threshold against its own chemical mode in the LEO mission reproduces its frozen value.",
        r"\num{2e-4}",
        "sweep",
    ),
    (
        "LEO mono_multimode minimum viable eta_u (single_chem)",
        "The monopropellant multimode threshold against the single-chemical architecture in the LEO mission reproduces its frozen value.",
        r"\num{2e-4}",
        "sweep",
    ),
    (
        "LEO crossover ordering",
        "In the LEO mission both monopropellant multimode thresholds lie below the corresponding bipropellant thresholds.",
        "exact (order)",
        "sweep",
    ),
    (
        "GEO mono_multimode eta_u grid",
        "In the GEO transfer mission the monopropellant multimode utilisation sweep diverges at exactly the frozen low-utilisation points.",
        "exact",
        "sweep",
    ),
    (
        "GEO mono_multimode min-viable readings",
        "In the GEO transfer mission the monopropellant multimode architecture has no threshold, because its own-chemical reference has no converged design.",
        "exact",
        "sweep",
    ),
    (
        "default eta_u path untouched",
        "The monopropellant sweep does not alter the reference set of the bipropellant sweep.",
        "exact",
        "sweep",
    ),
]
# bus-fraction comparison and tank type
SL += [
    (
        "LEO envelope never binds",
        "In the LEO mission no design lies above the envelope limit at limits of \\numlist{0.35;0.5;1.0}.",
        "exact",
        "sweep",
    ),
    (
        "GEO envelope verdicts",
        "In the GEO transfer mission the envelope readings of the four architectures at limits of \\numlist{0.35;0.5;1.0} reproduce the frozen table.",
        "exact",
        "sweep",
    ),
]
for k in ("single_elec", "multimode"):
    SL.append(
        (
            f"GEO COPV {k} wet",
            f"With a composite \\ce{{N2O}} tank, the {ARCH[k]} design of the GEO transfer mission reproduces its frozen wet mass.",
            KG,
            "system",
        )
    )
SL.append(
    (
        "GEO COPV multimode n_ep",
        "With a composite \\ce{N2O} tank, the bipropellant multimode design of the GEO transfer mission installs two electric units.",
        "exact",
        "system",
    )
)
for k in ("hybrid", "single_chem"):
    SL.append(
        (
            f"GEO COPV {k} wet",
            f"With a composite \\ce{{N2O}} tank, the {ARCH[k]} design of the GEO transfer mission reproduces its frozen wet mass.",
            KG,
            "system",
        )
    )
SL += [
    (
        "GEO COPV multimode - hybrid gap",
        "With a composite \\ce{N2O} tank, the bipropellant multimode minus hybrid difference in the GEO transfer mission reproduces its frozen value.",
        r"\SI{0.1}{\kilo\gram}",
        "system",
    ),
    (
        "GEO COPV ordering",
        "With a composite \\ce{N2O} tank, the bipropellant multimode design of the GEO transfer mission is lighter than the hybrid.",
        "exact (sign)",
        "system",
    ),
    (
        "toggle at metal",
        "Declaring every chemical tank metal reproduces the nominal designs of the GEO transfer mission exactly.",
        "exact",
        "system",
    ),
    (
        "toggle copv ==",
        "Declaring every chemical tank composite gives the same designs as declaring the \\ce{N2O} and propene tanks composite one by one.",
        "exact",
        "system",
    ),
]
# warning-time sweep
for k in ORDER5:
    SL.append(
        (
            f"t_warn verdicts {k}",
            f"The {ARCH[k]} feasibility verdicts at warning times of \\SIlist{{2;6;24;72}}{{\\hour}} reproduce the frozen set.",
            "exact",
            "sweep",
        )
    )
SL += [
    (
        "single_elec CA flip within",
        "The warning time below which the single-electric architecture is excluded lies within half an hour of its independently computed value.",
        r"\SI{0.5}{\hour}",
        "sweep",
    ),
    (
        "multimode 6 h spot wet",
        "The bipropellant multimode wet mass of the LEO mission at a warning time of \\SI{6}{\\hour} reproduces its frozen value.",
        KG,
        "sweep",
    ),
    (
        "CA-class propellant monotone",
        "The avoidance propellant rises as the warning time shortens, for the single-electric architecture wherever the required drift can be reached.",
        r"\SI{1e-12}{\kilo\gram}",
        "sweep",
    ),
    (
        "every t_warn grid closure converged",
        "Every point of the warning-time sweep converges.",
        "exact",
        "sweep",
    ),
    (
        "a wet-mass-optimal feasible contender",
        "A feasible design exists at every point of the warning-time sweep.",
        "exact",
        "sweep",
    ),
    (
        "optimal contender over",
        "The lightest feasible design at warning times of \\SIlist{2;6;24;72}{\\hour} reproduces the frozen set.",
        "exact",
        "sweep",
    ),
]
# transfer-window sweep
SL += [
    (
        "closed form hall at 3 mo",
        "In the closed-form transfer check at a window of three months, the Hall unit count, bus power, power-system increment, PPU mass and power verdict reproduce their frozen values.",
        r"\SI{0.01}{\kilo\watt}; \SI{0.1}{\kilo\gram}; counts and verdict exact",
        "sweep",
    ),
    (
        "closed form ecr at 3 mo",
        "The same holds for the electrodeless fleet at three months.",
        r"\SI{0.01}{\kilo\watt}; \SI{0.1}{\kilo\gram}; counts and verdict exact",
        "sweep",
    ),
    (
        "closed form hall at 24 mo",
        "The same holds for the Hall fleet at 24 months.",
        r"\SI{0.01}{\kilo\watt}; \SI{0.1}{\kilo\gram}; counts and verdict exact",
        "sweep",
    ),
    (
        "closed form ecr at 24 mo",
        "The same holds for the electrodeless fleet at 24 months.",
        r"\SI{0.01}{\kilo\watt}; \SI{0.1}{\kilo\gram}; counts and verdict exact",
        "sweep",
    ),
    (
        "transfer-window set",
        "The transfer-window sweep runs over one week, two weeks, one month and three, six, 12 and 24 months.",
        r"\num{0.01} months",
        "sweep",
    ),
    (
        "single_elec nominal n_Hall",
        "The single-electric Hall unit counts over the window set reproduce the frozen set, with no fixed point at the two shortest windows.",
        "exact",
        "sweep",
    ),
]
for k in ("single_elec", "multimode", "hybrid", "mono_multimode"):
    SL.append(
        (
            f"nominal T_max verdicts {k}",
            f"The {ARCH[k]} feasibility verdicts over the window set, under the time-optimal strategy, reproduce the frozen set.",
            "exact",
            "sweep",
        )
    )
SL += [
    (
        "multimode / hybrid nominal wets are window-invariant",
        "Under the time-optimal strategy the bipropellant multimode and hybrid designs do not depend on the transfer window.",
        KG,
        "sweep",
    ),
    (
        "multimode EP-transfer 24 mo split",
        "Under the maximum-electric strategy at a window of 24 months, the bipropellant multimode transfer share and wet mass reproduce their frozen values.",
        r"\num{5e-4} in share; \SI{0.05}{\kilo\gram}",
        "sweep",
    ),
    (
        "multimode EP-transfer excluded at every window",
        "Under the maximum-electric strategy the bipropellant multimode design lies above the envelope limit at every window.",
        "exact",
        "sweep",
    ),
    (
        "hybrid EP-transfer wets at",
        "Under the maximum-electric strategy the hybrid wet masses at three, six, 12 and 24 months reproduce their frozen values.",
        KG,
        "sweep",
    ),
    (
        "hybrid EP-transfer verdicts",
        "Under the maximum-electric strategy the hybrid feasibility verdicts over the window set reproduce the frozen set.",
        "exact",
        "sweep",
    ),
    (
        "mono_multimode nominal transfer converges",
        "Under the time-optimal strategy the monopropellant multimode design converges to the same mass and unit count at every window and lies above the envelope limit.",
        KG,
        "sweep",
    ),
    (
        "mono_multimode EP-transfer: a feasible split",
        "Under the maximum-electric strategy a power-feasible transfer share exists for the monopropellant multimode architecture at every window, and every resulting design lies above the envelope limit.",
        "exact",
        "sweep",
    ),
    (
        "mono_multimode EP-transfer 24 mo",
        "Under the maximum-electric strategy at 24 months, the monopropellant multimode transfer share and unit count reproduce their frozen values, the design converges within the power-system limit, and it lies above the envelope limit.",
        r"\num{5e-4} in share; count exact",
        "sweep",
    ),
    (
        "toggle off: ep_transfer reading empty",
        "With the maximum-electric reading switched off, the sweep returns no maximum-electric designs.",
        "exact",
        "sweep",
    ),
    (
        "toggle off: nominal records bit-identical",
        "Switching the maximum-electric reading on or off leaves the time-optimal designs unchanged.",
        "exact",
        "sweep",
    ),
    (
        "n_max under the power-system limit without degradation",
        "Without array degradation, the largest number of units that the power-system limit allows on the \\SI{300}{\\kilo\\gram} bus reproduces its frozen value for the Hall and electrodeless units.",
        "exact",
        "system",
    ),
    (
        "n_max under the power-system limit at the nominal derate",
        "At the nominal array degradation, applied in full as for a fully electric transfer, the same largest numbers of units reproduce their frozen values.",
        "exact",
        "system",
    ),
    (
        "hybrid 1 mo split share",
        "Under the maximum-electric strategy at a window of one month, the hybrid transfer share reproduces its frozen value.",
        r"\num{5e-4}",
        "sweep",
    ),
    (
        "hybrid 1 mo split wet",
        "That hybrid design reproduces its frozen wet mass.",
        KG,
        "sweep",
    ),
    (
        "hybrid 1 mo split n_hall == the unit ceiling at its own derate",
        "Under the maximum-electric strategy at a window of one month, the hybrid design installs as many Hall units as the power-system limit allows at its own degradation, the full degradation scaled by its electric share of the transfer.",
        "exact",
        "sweep",
    ),
    (
        "hybrid 1 mo split EPS",
        "The added power-system mass of that hybrid design reproduces its frozen value and lies within the power-system limit.",
        KG,
        "sweep",
    ),
    (
        "hybrid 1 mo split transfer EP propellant",
        "Its electric transfer propellant equals what the Hall fleet can expel within the window.",
        KG,
        "sweep",
    ),
    (
        "hybrid 1 mo split: every class verdict",
        "That hybrid design passes every schedule, duty and power check, and lies above the envelope limit.",
        "exact",
        "sweep",
    ),
    (
        "hybrid 24 mo: s* = 1.0",
        "At a window of 24 months the hybrid flies the whole transfer electrically and reproduces its frozen wet mass.",
        KG,
        "sweep",
    ),
    (
        "single_elec toggle allocation == nominal",
        "For the single-electric architecture the maximum-electric allocation is its nominal allocation.",
        "exact",
        "sweep",
    ),
]
# sensitivity method
SL += [
    (
        "pf_metal = 6350 m",
        "At a metal-tank performance factor of \\SI{6350}{\\metre} the single-chemical design of the LEO mission reproduces its frozen mass.",
        KG,
        "sweep",
    ),
    (
        "power_margin_150 = 50 W",
        "Lowering the LEO power margin to \\SI{50}{\\watt} raises the single-electric wet mass by its independently computed amount.",
        r"\SI{0.01}{\kilo\gram}",
        "sweep",
    ),
    (
        "m_head_ecr override propagates",
        "Raising the electrodeless head mass to \\SI{1.5}{\\kilo\\gram} raises the bipropellant multimode wet mass of the LEO mission by its independently computed amount, so a changed input reaches the unit it describes.",
        r"\SI{0.01}{\kilo\gram}",
        "sweep",
    ),
    (
        "every 'oat' register row appears",
        "Every input marked for the one-at-a-time set appears exactly once among the sensitivity cases.",
        "exact",
        "sweep",
    ),
    (
        "the three thruster band bars",
        "The three bands on the electrodeless unit appear among the sensitivity cases.",
        "exact",
        "sweep",
    ),
    (
        "p_heater_mono is an oat bar",
        "The catalyst heater power appears among the sensitivity cases.",
        "exact",
        "sweep",
    ),
    (
        "p_heater_mono bar zero-width",
        "Varying the catalyst heater power from \\SIrange{10}{50}{\\watt} does not change the monopropellant multimode wet mass of the LEO mission.",
        KG,
        "sweep",
    ),
    (
        "d(mono_mm - mm) tornado baseline == the LEO difference",
        "The baseline of the sensitivity cases reproduces the monopropellant minus bipropellant multimode difference of the LEO mission.",
        KG,
        "sweep",
    ),
    (
        "d(mono_mm - mm) GEO tornado baseline",
        "The baseline of the sensitivity cases reproduces the frozen monopropellant minus bipropellant multimode difference of the GEO transfer mission.",
        r"\SI{0.05}{\kilo\gram}",
        "sweep",
    ),
    (
        "parameter values restored after the tornado",
        "After the sensitivity cases, the nominal bipropellant multimode design of the LEO mission is reproduced, so no varied input remains in force.",
        KG,
        "sweep",
    ),
    (
        "corner grid levels read",
        "The corner check takes the three levels of the monopropellant head mass and specific impulse from their nominal values and bounds.",
        "exact",
        "sweep",
    ),
    (
        "nominal grid point bit-reproduces",
        "The nominal point of the corner check reproduces the five nominal LEO designs.",
        r"\SI{1e-9}{\kilo\gram}",
        "sweep",
    ),
    (
        "mono_multimode wet matches the reference grid",
        "At all nine corner points the monopropellant multimode wet mass reproduces its frozen value.",
        KG,
        "sweep",
    ),
    (
        "all 45 grid closures converged",
        "All 45 designs of the corner check converge.",
        "exact",
        "sweep",
    ),
    (
        "the four non-mono architectures are grid-invariant",
        "The other four architectures do not change across the corner check.",
        r"\SI{1e-9}{\kilo\gram}",
        "sweep",
    ),
    (
        "unfavourable corner (head hi, isp lo)",
        "At the least favourable corner, the monopropellant minus bipropellant multimode difference reproduces its frozen value.",
        KG,
        "sweep",
    ),
    (
        "interaction term == affine cross term",
        "The interaction term of the corner check reproduces the cross term of the affine closure.",
        r"\SI{2e-4}{\kilo\gram}",
        "sweep",
    ),
    (
        "unfavourable-corner lead stays negative",
        "At the least favourable corner the monopropellant multimode architecture remains the lighter of the two.",
        "exact (sign)",
        "sweep",
    ),
    (
        "mono_multimode lightest of the five",
        "The monopropellant multimode architecture is the lightest of the five at every corner point.",
        "exact",
        "sweep",
    ),
    (
        "parameter values restored after the corner",
        "After the corner check, the nominal bipropellant multimode design of the LEO mission is reproduced.",
        KG,
        "sweep",
    ),
]
# input consistency
SL += [
    (
        "GEO transfer class EP budget",
        "The low-thrust transfer budget equals the chemical budget divided by the planning efficiency.",
        r"\SI{1e-9}{\metre\per\second}",
        "system",
    ),
    (
        "a dv_transfer_chem override carries",
        "Changing the chemical transfer budget changes the low-thrust budget and the closed-form thrust requirement in the same proportion.",
        r"relative \num{1e-12}",
        "system",
    ),
]

# Thruster-model checks
CHK = []
for p, tex in (
    ("Xe", r"\ce{Xe}"),
    ("N2O", r"\ce{N2O}"),
    ("H2O", r"\ce{H2O}"),
    ("Ar", r"\ce{Ar}"),
):
    for t in (2, 5, 10, 22):
        CHK.append(
            (
                f"Ec {p} @ {t} eV",
                f"The ionisation cost of {tex} at an electron temperature of \\SI{{{t}}}{{\\electronvolt}} reproduces its frozen value.",
                r"relative \num{1e-6}",
                "thruster",
            )
        )
CHK += [
    (
        "Xe T [N]",
        "At the measured xenon utilisation, the thrust of the xenon anchor reproduces the measured operating point.",
        r"relative \num{1e-6}",
        "thruster",
    ),
    (
        "Xe m_dot",
        "The mass flow rate of the xenon anchor reproduces the measured operating point.",
        r"relative \num{1e-6}",
        "thruster",
    ),
    (
        "Xe P_in",
        "The absorbed power of the xenon anchor reproduces the measured operating point.",
        r"relative \num{1e-6}",
        "thruster",
    ),
    (
        "Xe Isp_gross",
        "The specific impulse of the xenon anchor unit equals the value implied by its measured thrust and mass flow rate.",
        r"relative \num{1e-6}",
        "thruster",
    ),
    (
        "Xe R_T",
        "The xenon thrust ratio at the anchor is one.",
        r"relative \num{1e-9}",
        "thruster",
    ),
    (
        "Xe R_mdot_ion",
        "The xenon ion mass-flow ratio at the anchor is one.",
        r"relative \num{1e-9}",
        "thruster",
    ),
    (
        "Xe R_ue",
        "The xenon exhaust-velocity ratio at the anchor is one.",
        r"relative \num{1e-9}",
        "thruster",
    ),
    (
        "Xe R_Pin",
        "The xenon absorbed-power ratio at the anchor is one.",
        r"relative \num{1e-9}",
        "thruster",
    ),
    (
        "N2O_B T [N]",
        "The thrust of the \\ce{N2O} unit at utilisation parity reproduces its frozen value.",
        r"relative \num{1e-6}",
        "thruster",
    ),
    (
        "N2O_B m_dot",
        "The mass flow rate of the \\ce{N2O} unit at utilisation parity reproduces its frozen value.",
        r"relative \num{1e-6}",
        "thruster",
    ),
    (
        "N2O_B P_in",
        "The absorbed power of the \\ce{N2O} unit at utilisation parity reproduces its frozen value.",
        r"relative \num{1e-6}",
        "thruster",
    ),
    (
        "N2O_B Isp_gross",
        "The specific impulse of the \\ce{N2O} unit at utilisation parity reproduces its frozen value.",
        r"relative \num{1e-6}",
        "thruster",
    ),
    (
        "N2O_B Isp_intensive",
        "The ion-beam specific impulse of the \\ce{N2O} unit reproduces its frozen value.",
        r"relative \num{1e-6}",
        "thruster",
    ),
    (
        "Eq. 17",
        "The two forms of the thrust efficiency of Petro and Sedwick, their Eqs. (17) and (18), agree for mirror ratios of 1, 10 and 100.",
        r"\num{1e-10}",
        "thruster",
    ),
    (
        "Analytic limit",
        "With a mirror ratio of one and no ionisation cost, the thrust efficiency equals its analytic limit.",
        r"\num{1e-12}",
        "thruster",
    ),
    (
        "Ar vs P&S @ 5 eV",
        "The argon ionisation cost at \\SI{5}{\\electronvolt} lies \\SI{1.85}{\\percent} below the curve fit of Petro and Sedwick.",
        r"\SI{0.15}{\percent}",
        "thruster",
    ),
    (
        "Ar vs P&S @ 10 eV",
        "The argon ionisation cost at \\SI{10}{\\electronvolt} lies \\SI{0.32}{\\percent} below the curve fit of Petro and Sedwick.",
        r"\SI{0.15}{\percent}",
        "thruster",
    ),
]

# Spacecraft check sections, in order; each is one group of the appendix
GROUPS = (
    "LEO mission results",
    "Differences between architectures",
    "Tank sizing, LEO",
    "Power of the two electric units",
    "Allocation limits",
    "Allocation sweep",
    "GEO transfer mission results",
    "Collision-avoidance check",
    "Allocation sweep structure",
    "Interface to the thruster-level model",
    "Utilisation sweep",
    "Bus-fraction comparison and tank type",
    "Warning-time sweep",
    "Transfer-window sweep",
    "Sensitivity method",
    "Input consistency",
)

# Test levels, in the order of the appendix tables
LEVELS = (
    (
        "unit",
        "Unit tests",
        "A single relation or input agrees with an independent hand calculation, a published value or its specification",
    ),
    (
        "integration",
        "Integration tests",
        "The chained components converge and agree with one another or with an independent value",
    ),
    (
        "regression",
        "Regression tests",
        "Results reproduce frozen values from an earlier run, so that no change goes unnoticed",
    ),
    (
        "analytical",
        "Analytical tests",
        "The model reproduces a closed-form solution or a limiting case",
    ),
    (
        "property",
        "Property tests",
        "Monotonicity, sign, ordering, threshold crossings, invariance and completeness hold",
    ),
)


def _ids(prefix, spec):
    """'S1.1-5, 7' style -> list of IDs; ranges inclusive."""
    out = []
    for part in spec.split(","):
        part = part.strip()
        if "-" in part:
            a, b = part.split("-")
            out += [f"{prefix}{n}" for n in range(int(a), int(b) + 1)]
        else:
            out.append(f"{prefix}{part}")
    return out


_LEVEL_SPEC = {
    # thruster suite: E_c and N2O values are regression targets; the xenon
    # anchor reproduces Cannat et al. (2015) Table II; Eqs. (17) and (18) and
    # the argon fit are published relations; the anchor ratios and the
    # R_m = 1, E_c = 0 limit are analytic
    "unit": [
        ("T.", "17-20, 30, 32-33"),
        ("S10.", "2-3, 5-6"),  # Isp ~ eta_u, mdot ~ 1/eta_u
        ("S14.", "5"),  # transfer-window set
        ("S15.", "11"),  # corner levels read from the register brackets
        ("S16.", "1-2"),
    ],  # derived budget and its override
    "integration": [
        ("S1.", "2, 4, 6, 8, 10"),
        ("S2.", "4-8"),
        ("S6.", "2"),
        ("S7.", "2, 5, 8, 11, 14"),
        ("S9.", "1"),
        ("S10.", "8"),
        ("S11.", "1, 18"),
        ("S13.", "9"),
        ("S14.", "25"),
        ("S15.", "2-3, 14"),
    ],  # S15.2-3: input -> closure -> mass
    "regression": [
        ("T.", "1-16, 25-29"),
        ("S1.", "1, 3, 5, 7, 9"),
        ("S2.", "1-3"),
        ("S3.", "1-3"),
        ("S4.", "1-9, 11"),
        ("S5.", "2-3"),
        ("S7.", "1, 3-4, 6-7, 9-10, 12-13, 15-21"),
        ("S8.", "5-6"),
        ("S9.", "2-3"),
        ("S11.", "2-7, 9, 11-16, 19-22, 24"),
        ("S12.", "2-8"),
        ("S13.", "1-5, 7, 11"),
        ("S14.", "1-4, 6-10, 12, 14-15, 18, 21-24, 26-29"),
        ("S15.", "1, 8-10, 13, 16-17, 20"),
    ],
    "analytical": [
        ("T.", "21-24, 31"),
        ("S5.", "1, 4"),
        ("S9.", "5-6"),
        ("S10.", "1, 9"),
        ("S13.", "6"),
    ],
    "property": [
        ("S3.", "4"),
        ("S4.", "10"),
        ("S6.", "1"),
        ("S8.", "1-4"),
        ("S9.", "4"),
        ("S10.", "4, 7"),
        ("S11.", "8, 10, 17, 23, 25-26"),
        ("S12.", "1, 9-11"),
        ("S13.", "8, 10"),
        ("S14.", "11, 13, 16-17, 19-20, 30"),
        ("S15.", "4-7, 12, 15, 18-19"),
    ],
}
LEVEL = {}
for _lev, _parts in _LEVEL_SPEC.items():
    for _pre, _spec in _parts:
        for _i in _ids(_pre, _spec):
            assert _i not in LEVEL, ("two levels", _i)
            LEVEL[_i] = _lev
