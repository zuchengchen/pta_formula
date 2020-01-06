# Constrain GW polarization from NANOGrav 11 year data set

* One should use enterprise from
    https://github.com/bunnyandbear/enterprise

## Inflation model
* number of frequency bins is set to 50 (or 30)
* priors range for A_TT, A_ST, A_VL, A_SL are set in "models/models.py -> modelInflation" and "utils.py -> JumpProposal -> draw_from_gwb_log_uniform_distribution"


## SMBH model

* number of frequency bins is set to 50 (or 30)
* priors range for A_TT, A_ST, A_VL, A_SL are set in "models/models.py -> modelSMBH" and "analysis-SMBH.py -> JumpProposal -> draw_from_gwb_log_uniform_distribution"