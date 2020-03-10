def filterredwv(x, sigma_r, sigma_w, zeropad=None):
    #+
    # NAME:
    #       FILTERREDWV
    # PURPOSE:
    #       Returns the 1/f component from a data vector x that is assumed
    #       to be described for parameters sigma_r and sigma_w as 1/f
    #       noise plus white noise as detailed by Carter & Winn (2009)
    #
    # CALLING SEQUENCE:
    #       result = filterredwv(x, sigma_r, sigma_w)
    #
    # INPUTS:
    #       X - Data vector, must be a power of two
    #       SIGMA_R - Red noise amplitude parameter (not equal to
    #                 1/f-component RMS)
    #       SIGMA_W - White noise amplitude parameter (approximately equal
    #                 to white-component RMS)
    #
    # OUTPUTS:
    #       RESULT - 1/f component with same number of element as for input X.
    #
    # NOTES:
    #       Refer to the paper by Carter & Winn (2009) for theory and
    #       details.  In the notation of that paper, gamma=1 for this
    #       algorithm.
    #
    # REVISION HISTORY:
    #       Written,    J. Carter               September, 2009
    #-

    gamma = 1
    duplicate = x
    els = len(duplicate)
    pow = np.ceil(np.log(els) / np.log(2.))
    if (2**pow != els and zeropad is not None):
        diff = 2.**pow - els
        left = np.floor(diff / 2.)
        right = diff - left

        if diff > 1:
            x = [np.zeros(left), duplicate, np.zeros(right)]
            else:
            x = [duplicate, np.zeros(right)]
    else:
        x = duplicate

    J = np.log(len(x)) / np.log(2)
    if np.abs(J - int(J)) != 0:
        'Data length must be a power of two'

    J = int(J)

    info = wv_fn_daubechies(2, wavelet, scaling, ioff, joff)
    wv = wv_dwt(x, wavelet, scaling, ioff, joff)

    sm2 = sigma_r**2 * (gamma == 1 ? 1.0 / (2.0 * np.log(2.0)): 2 - 2 ** gamma) + sigma_w ** 2

    wv[0] *= 1 - (sigma_w ** 2) / sm2

    k = 1.
    for i in np.arange(J):
        sm2 = sigma_r**2 * 2**(-gamma * i * (1.0)) + sigma_w ** 2
        # for m = 0., 2 ** (i - 1) - 1, 1
        for m in np.arange(-1, 2**(i - 1) - 1):
            wv[k] *= 1 - sigma_w ** 2 / sm2
            k += 1

    redcomp = wv_dwt(wv, wavelet, scaling, ioff, joff, inverse=True)

    if zeropad is not None:
        if len(left) != 0:
            if left != 0 and right != 0:
                x = duplicate

                # whitecomp = whitecomp[left:els+left]
                redcomp = redcomp[left: (els + left)]

    duplicate = 0

    return redcomp
