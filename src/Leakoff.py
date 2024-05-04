import pandas as pd
import math

class leakOff:
    def __init__(self, isOIlField=True):
        self.PI = 3.141592654
        self.isOIlField = isOIlField

    def C_c(self, dP_T, k_resfluid, por, ct, visc_resfluid):
        """
        Corrected value of conductance fluid-loss coefficient.

        Params
        -------
        dP_T: float
            Total pressure loss
        k_resfluid: float
            permeability of reservoir fluid
        por: float
            porosity
        ct: float
            Compressibility of the reservoir mobile phase
        visc_resfluid: float
            viscosity of reservoir fluid

        Returns
        -------
        C_c: float
            Corrected value of C
        """
        if self.isOIlField:
            C_c = 0.00118 * dP_T * math.sqrt(k_resfluid * por * ct / visc_resfluid)

        else:
            C_c = dP_T * math.sqrt(1 / self.PI * k_resfluid * por * ct / visc_resfluid)

        return C_c
    
    def C_v(self, dP_T, k_resfluid, k_rel, por_eff_filt, visc_filtrate):
        """
        Calculate viscosity leakoff control coefficient

        Params
        -------
        dP_T: float
            Total pressure loss
        k_resfluid: float
            Permeability of reservoir fluid
        k_rel: float
            Relative permeability
        por_eff_filt: float
            Effective filtrate porosity
        visc_flitrate: float
            Filtrate viscosity

        Returns
        -------
        C_v: float
            Viscosity leakoff control coefficient
        """
        if self.isOIlField:
            C_v = 0.00148 * math.sqrt(dP_T * k_resfluid * k_rel * por_eff_filt / visc_filtrate)
        else:
            C_v = math.sqrt(0.5 * dP_T * k_resfluid * k_rel * por_eff_filt / visc_filtrate)
        
        return C_v
    
    def C_vc(self, Cc, Cv):
        """
        Combination of C_v and C_c, restores contribution of dP on each zone invaded by filtrate (dP_v) and reservoir (dP_c)

        Params
        -------
        Cc: float
            Corrected value of conductance fluid loss coef.
        Cv: float
            Viscosity leakoff control coefficient

        Returns
        -------
        C_vc: float
            Combination of C_v and C_c
        """
        C_vc = 2 * Cc / (1 + math.sqrt(1 + 4 * (Cc/Cv) ** 2))

        return C_vc
    
    def spurtCorrelation(self, dP_T, k_resfluid, visc_filtrate, guar_conc, isXL, isCleanGuar):
        """
        Calculate approximate spurt loss. If result is < 2 gal/100ft2, then the result is negligible.

        Params
        -------
        dP_T: float
            Total pressure loss
        k_resfluid: float
            Permeability of reservoir fluid
        visc_filtrate: float
            Filtrate viscosity
        guar_conc: float
            Guar concentration
        isXL: bool
            Tests whether the injected fluid is crosslinked or linear
        isCleanGuar: bool
            Tests whether the injected fluid is guar or HPG (hydroxypropyl guar)

        Returns
        -------
        spurt_coor: float
            Estimated spurt loss
        """
        if isXL == True or isXL == 'Y':
            isXL = True
        else:
            isXL = False

        if isCleanGuar == True or isCleanGuar == 'Y':
            isCleanGuar = True
        else:
            isCleanGuar = False
        
        XL_factor = 3 if isXL else 1
        residueFactor = 0.25 if isCleanGuar else 1

        m = 0.75
        a = 660000000#
        b = 1.75
        c = 0.8
        e = 4

        spurt_corr = a * ((dP_T / 1000) ** b) * ((k_resfluid / 1000) ** c) / (residueFactor * visc_filtrate) / ((XL_factor ** m) * guar_conc) ** e
        
        return spurt_corr