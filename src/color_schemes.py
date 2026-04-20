"""Centralized color schemes for figures.

Switch ACTIVE_COLOR_SCHEME_NAME to change all figure colors in one place.
"""

from matplotlib.colors import LinearSegmentedColormap


GROUP_ORDER = [
    ("replete", "10NaCl"),
    ("replete", "45NaCl"),
    ("deplete", "10NaCl"),
    ("deplete", "45NaCl"),
]


COLOR_SCHEMES = {
    "paper_default": {
        "group_colors": {
            ("replete", "10NaCl"): "#67AFD2",
            ("replete", "45NaCl"): "#016895",
            ("deplete", "10NaCl"): "#F4795B",
            ("deplete", "45NaCl"): "#C74632",
        },
        "diverging": {
            "negative": "#016895",
            "neutral": "white",
            "positive": "#C74632",
        },
    },
    "high_contrast": {
        "group_colors": {
            ("replete", "10NaCl"): "#4DA3C8",
            ("replete", "45NaCl"): "#004E73",
            ("deplete", "10NaCl"): "#F08A24",
            ("deplete", "45NaCl"): "#B7371F",
        },
        "diverging": {
            "negative": "#004E73",
            "neutral": "white",
            "positive": "#B7371F",
        },
    },
        "novel": {
        "group_colors": {
            ("replete", "10NaCl"): "#48BF84", ##61D095
            ("replete", "45NaCl"): "#439775",
            ("deplete", "10NaCl"): "#AD6A6C",
            ("deplete", "45NaCl"): "#5D2E46",
        },
        "diverging": {
            "negative": "#016895",
            "neutral": "white",
            "positive": "#C74632",
        },
    },
}


ACTIVE_COLOR_SCHEME_NAME = "novel"


def list_color_schemes():
    """Return available scheme names."""
    return sorted(COLOR_SCHEMES.keys())


def get_color_scheme(name=None):
    """Return a named color scheme or the active default scheme."""
    scheme_name = name or ACTIVE_COLOR_SCHEME_NAME
    if scheme_name not in COLOR_SCHEMES:
        available = ", ".join(list_color_schemes())
        raise ValueError(f"Unknown color scheme: {scheme_name}. Available: {available}")
    return COLOR_SCHEMES[scheme_name]


def make_diverging_colormap(name=None, cmap_name="custom_diverging", bins=256):
    """Create a diverging colormap from scheme negative/neutral/positive colors."""
    div = get_color_scheme(name)["diverging"]
    return LinearSegmentedColormap.from_list(
        cmap_name,
        [div["negative"], div["neutral"], div["positive"]],
        N=bins,
    )
