"""
plotting.py
-----------
Matplotlib plotting utilities shared across the NIRCam WFSS pipeline.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt


def get_corner_pos(
    ax: plt.Axes,
    loc: int = 1,
    edge: float = 0.02,
) -> tuple[float, float]:
    """
    Return the (x, y) data-coordinates of a corner of ``ax``.

    Parameters
    ----------
    ax:
        Matplotlib ``Axes`` object.
    loc:
        Corner location: 1=upper-right, 2=upper-left, 3=lower-left,
        4=lower-right.
    edge:
        Fractional inset from the axes boundary (0–1).

    Returns
    -------
    (x_text, y_text):
        Data-space coordinates of the inset corner.
    """
    if loc not in (1, 2, 3, 4):
        raise ValueError("loc must be 1, 2, 3, or 4")

    xlim   = ax.get_xlim()
    ylim   = ax.get_ylim()
    xscale = ax.get_xscale()
    yscale = ax.get_yscale()

    vec_edge = np.array([1 - edge, edge])
    vec_x = np.array(
        [
            [xlim[1], xlim[0], xlim[0], xlim[1]][loc - 1],
            [xlim[0], xlim[1], xlim[1], xlim[0]][loc - 1],
        ]
    )
    vec_y = np.array(
        [
            [ylim[1], ylim[1], ylim[0], ylim[0]][loc - 1],
            [ylim[0], ylim[0], ylim[1], ylim[1]][loc - 1],
        ]
    )

    if xscale == "linear":
        x_text = np.dot(vec_edge, vec_x)
    else:
        x_text = 10 ** np.dot(vec_edge, np.log10(vec_x))

    if yscale == "linear":
        y_text = np.dot(vec_edge, vec_y)
    else:
        y_text = 10 ** np.dot(vec_edge, np.log10(vec_y))

    return (x_text, y_text)


def corner_text(
    ax: plt.Axes,
    s: str = "",
    loc: int = 1,
    edge: float = 0.02,
    **kwargs,
) -> plt.Text:
    """
    Place a text label at a corner of ``ax``.

    Parameters
    ----------
    ax:
        Matplotlib ``Axes`` object.
    s:
        Text string to display.
    loc:
        Corner: 1=upper-right, 2=upper-left, 3=lower-left, 4=lower-right.
    edge:
        Fractional inset from the axes boundary.
    **kwargs:
        Additional keyword arguments forwarded to ``ax.text()``.

    Returns
    -------
    text:
        The ``matplotlib.text.Text`` object created.
    """
    if loc not in (1, 2, 3, 4):
        raise ValueError("loc must be 1, 2, 3, or 4")

    ha_map = ["right",  "left",  "left",  "right"]
    va_map = ["top",    "top",   "bottom", "bottom"]

    return ax.text(
        *get_corner_pos(ax, loc, edge),
        s=s,
        ha=ha_map[loc - 1],
        va=va_map[loc - 1],
        **kwargs,
    )
