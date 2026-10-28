"""
Chemical mode: N2O/propene oxidiser fraction from the O/F ratio.
"""

from .parameters import value


def oxidiser_fraction():
    """N2O mass fraction of the biprop load, O/F / (O/F + 1)."""
    of = value("of_ratio")
    return of / (of + 1.0)


def split_biprop(m_chem):
    """Biprop propellant mass split into (m_N2O, m_propene)."""
    fox = oxidiser_fraction()
    return fox * m_chem, (1.0 - fox) * m_chem
