def solveredwv(data, s_r, s_w, silent=None, zeropad=None):
    #+
    # NAME:
    #       SOLVEREDWV
    # PURPOSE:
    #       Given a data vector, data, this algorithm attempts to find the best fit model
    #       describing the data as 1/f noise plus white noise
    # EXPLANATION:
    #       Uses the wavelet technique as described by Carter & Winn
    #       (2009) to determine parameters of noise formed as an additive
    #       combination of (Gaussian) noise with power spectral density proportional to 1/f and
    #       (Gaussian) white noise.  In particular, AMOEBA is used to maximize a
    #       likelihood that is a function of two parameters sigma_r and
    #       sigma_w which are related to the standard deviations of the 1/f
    #       and white components, respectively.
    #
    #       Additionally, a wavelet filter is applied to separate the two
    #       components.  These results may be returned to the caller as optional outputs.
    #
    # CALLING SEQUENCE:
    #       solveredwv,data, s_r, s_w,[redcomp=redcomp, whitecomp=whitecomp,
    #       alpha=alpha, sol=sol, /silent]
    #
    # INPUTS:
    #       data    - Data vector, length must be a power of two
    #       SIGMA_R - Red noise amplitude (not equal to 1/f-component
    #                 RMS): If 2-element vector, then this input gives the
    #                 approximate range for SIGMA_R in which the solution
    #                 lies.  If scalar, SIGMA_R is fixed to its input
    #                 value (i.e., it does not vary in AMOEBA).
    #       SIGMA_W - White noise amplitude parameter (approximately equal
    #                 to white-component RMS): If 2-element vector, then this input
    #                 input gives the approximate range for SIGMA_W in
    #                 which the solution lies.  If scalar, SIGMA_W is
    #                 fixed to its input value (i.e., it does not vary in
    #                 AMOEBA).
    # OPTIONAL INPUTS:
    #       SILENT - Suppresses output.
    #
    # OUTPUTS:
    #       Prints the best fit values for SIGMA_R and SIGMA_W.  Also
    #       reports the RMS ratio of the 1/f component to the white
    #       component [referred to as "alpha" by Carter & Winn (2009)]
    #
    # OPTIONAL OUTPUTS:
    #       REDCOMP - The best fit 1/f-component solution for data vector data (of the
    #                 same length as data).
    #       WHITECOMP - The best fit white-component solution for data
    #                 vector data (of the same length as data).
    #       SOL     - 2-element vector of the best fit SIGMA_R (SOL[0]) and
    #                 SIGMA_W (SOL[1])
    #       ALPHA   - Number giving the RMS ratio of the 1/f component to
    #                 the white component.
    #
    # EXAMPLES:
    #
    #       EXAMPLE 1: White noise.
    #
    #       IDL> x = randomn(seed,1024)
    #       IDL> solveredwv,x,[0.0,1.0],[0.5,1.5],sol=sol,alpha=alpha,$
    #            redcomp=redcomp,whitecomp=whitecomp
    #       Sigma_r = -0.032677
    #       Sigma_w = 0.994064
    #       Gamma = 1.000000
    #       RMS of white component is 0.993101
    #       RMS of red component is 0.000021
    #       Ratio of red to white RMS is 0.000021
    #       IDL>
    #
    #       EXAMPLE 2: Correlated noise (RMS = 0.40) + white noise (RMS = 1.0)
    #
    #       IDL> X = smooth(randomn(seed,1024),4)+randomn(seed,1024)
    #       IDL> solveredwv,x,[0.0,1.0],[0.5,1.5],sol=sol,alpha=alpha,$
    #            redcomp=redcomp,whitecomp=whitecomp
    #       Sigma_r = 6.238085
    #       Sigma_w = 1.000332
    #       Gamma = 1.000000
    #       RMS of white component is 0.951405
    #       RMS of red component is 0.215253
    #       Ratio of red to white RMS is 0.226248
    #
    # NOTES:
    #       Refer to the paper by Carter & Winn (2009) for theory and
    #       details.  In the notation of that paper, gamma=1 for this
    #       algorithm.
    #
    # REVISION HISTORY:
    #       Written,    J. Carter               September, 2009
    #       Made default zero-padding option    November, 2009
    #       Corrected non-short-circuited ifs!
    #               Thanks to E. Ford           January, 2010
    #-

    # common solvewv, x, var, sigma_r, sigma_w, gamma

    sigma_r = s_r
    sigma_w = s_w
    gamma = 1

    els = len(data)
    pow = np.ceil(alog(els) / np.alog(2.))
    if (2**pow != els and zeropad is not None) then begin
        diff = 2**pow - els
        left = np.floor(diff / 2.)
        right = diff - left
        x = [np.zeros(left), data, np.zeros(right)]
    else:
        x = data

    var = np.zeros(3)

    if (len(sigma_r) != 1):
        var[0] = 1
    if (len(sigma_w) != 1):
        var[1] = 1
    if (len(gamma) != 1):
        var[0] = 1

    k = 0
    par = np.arange(total(var))
    scale = np.arange(total(var))
    if (var[0] == 1) then begin
        par[k] = mean(sigma_r)
        # scale(k++) = sigma_r(1) - sigma_r(0)
        scale[k + 1] = sigma_r[1] - sigma_r[0]

    if var[1] == 1:
        par[k] = np.mean(sigma_w)
        # scale(k + +) = sigma_w(1) - sigma_w(0)
        scale[k + 1] = sigma_w[1] - sigma_w[0]

    if var[2] == 1:
        par[k] = np.mean(gamma)
        scale[k + 1] = gamma[1] - gamma[0]

    # result = amoeba(1.0e-5, function_name='minfuncwv', p0=par, scale=scale)
    # FINDME: this needs to be fixed
    result = scipy.optimize.minimize(minfuncwv, 1.0e-5, p0=par, scale=scale)
    redcomp = filterredwv(x, sigma_r, sigma_w)

    if zeropad is not None:
        if left != 0 and right != 0:
            x = data

            # whitecomp = whitecomp(left:(els+left))
            redcomp = redcomp[left: els + left]

    whitecomp = x - redcomp

    if silent is False:
        print(f"Sigma_r = {sigma_r}")
        print(f"Sigma_w = {sigma_w}")
        print(f"Gamma = {gamma}")

        print(f"RMS of white component is {np.std(whitecomp)}")
        print(f"RMS of red component is {np.std(redcomp)}")
        print("Ratio of red to white RMS is "
              f"{np.std(redcomp) / np.std(whitecomp)}")

    sol = [sigma_r, sigma_w]
    alpha = np.std(redcomp) / np.std(whitecomp)

    return redcomp, whitecomp, alpha, sol
