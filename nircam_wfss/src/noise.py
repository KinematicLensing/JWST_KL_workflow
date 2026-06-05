"""
noise.py
--------
1/f noise (correlated read-noise stripe) subtraction and robust background
estimation for NIRCam WFSS grism images.
"""

from __future__ import annotations

import warnings

import numpy as np
from astropy.stats import sigma_clipped_stats


def my_1overf_subtraction(
    sci_data: np.ndarray,
    backg_mask: np.ndarray | None = None,
    amplifiers: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Remove 1/f correlated read-noise stripes by subtracting the per-row median.

    Parameters
    ----------
    sci_data:
        2-D science image (rows × cols).
    backg_mask:
        Boolean or integer mask array of the same shape as ``sci_data``.
        Pixels with value 1 are treated as background and used to compute
        the row median.  If ``None``, all finite pixels are used.
    amplifiers:
        Number of HAWAII-2RG amplifier channels to handle separately.
        Use ``1`` for a single median over the full row; use ``4`` for
        the standard four-channel layout (512 columns each).

    Returns
    -------
    corrected:
        ``sci_data`` with the 1/f model subtracted.
    model_img:
        The 2-D model image that was subtracted (row-constant stripes).
    """
    if backg_mask is None:
        backg_mask = np.ones_like(sci_data)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN rows are expected for masked/edge rows
        if amplifiers == 1:
            model_img = np.array(
                [
                    np.ones(sci_data.shape[1])
                    * np.nanmedian(sci_data[i][backg_mask[i] == 1])
                    for i in range(sci_data.shape[0])
                ]
            )
        elif amplifiers == 4:
            idx_row_blocks = [
                np.arange(512),
                np.arange(512) + 512,
                np.arange(512) + 512 * 2,
                np.arange(512) + 512 * 3,
            ]
            blocks = []
            for tmp_idx_row in idx_row_blocks:
                block = np.array(
                    [
                        np.ones_like(tmp_idx_row)
                        * np.nanmedian(
                            sci_data[i, tmp_idx_row][backg_mask[i, tmp_idx_row] == 1]
                        )
                        for i in range(sci_data.shape[0])
                    ]
                )
                blocks.append(block)
            model_img = np.hstack(blocks)
        else:
            raise ValueError("amplifiers must be 1 or 4")

    return sci_data - model_img, model_img


def robust_median_bkg(
    img: np.ndarray,
    return_mask: bool = False,
) -> float | tuple[float, np.ndarray]:
    """
    Estimate a robust global background by masking the brightest 25th percentile
    of rows and columns.

    Parameters
    ----------
    img:
        2-D image array.
    return_mask:
        If ``True``, return ``(background, mask)`` where ``mask`` is a 2-D
        binary array (0 = masked, 1 = used for background estimation).

    Returns
    -------
    background:
        Sigma-clipped median background value.
    mask (optional):
        Binary mask array, only returned when ``return_mask=True``.
    """
    arr_x_avg = np.nanpercentile(img, 75, axis=0)
    arr_y_avg = np.nanpercentile(img, 75, axis=1)
    img_mask = np.ones_like(img)
    img_mask[np.where(arr_y_avg < np.nanmedian(arr_y_avg))[0], :] = 0
    img_mask[:, np.where(arr_x_avg < np.nanmedian(arr_x_avg))[0]] = 0

    tmp_global_bkg = sigma_clipped_stats(
        img[img_mask == 0].flatten(), sigma_upper=2.0, maxiters=5
    )[1]

    if return_mask:
        return tmp_global_bkg, img_mask
    return tmp_global_bkg
