# License
This data was generated using data from the Redistricting Data Hub.  Any use of this project shall also comply with restrictions on use of data and attribution requirements set forth in the Redistricting Data Hub terms and conditions found at: [https://redistrictingdatahub.org/terms-and-conditions/](https://redistrictingdatahub.org/terms-and-conditions/).

Use of this project is further governed by the terms of the [Creative Commons Attribution Noncommercial 4.0 International](https://creativecommons.org/licenses/by-nc/4.0/legalcode.en)

# Texas Json

This json and shapefile were created by Professor Ellen Veomett and her student Ananya Agarwal using the corresponding jupyter notebook.  As part of the cleaning process, precincts were nested within counties and small rook adjacencies (under 30.5 m) were changed to queen adjacencies.

# **Sources**

Obtain the following data from Restricting Data Hub

[Population data](https://redistrictingdatahub.org/dataset/texas-block-pl-94171-2020-by-table/): based on the decennial census at the Census Block level on 2020 Census Redistricting Data

[Congressional District data](https://redistrictingdatahub.org/dataset/2021-texas-congressional-districts-adopted-plan/): 2021 Texas Congressional Districts plan enacted on 11/11/21

[State House District data](https://redistrictingdatahub.org/dataset/2021-texas-state-house-adopted-plan/): 2021 State House Approved Plan

[State Senate District data](https://redistrictingdatahub.org/dataset/2021-texas-state-senate-adopted-plan/): 2021 State Senate Approved Plan

[2020 election data](https://redistrictingdatahub.org/dataset/vest-2020-texas-precinct-boundaries-and-election-results/)**:**  VEST 2020 Texas precinct and election results

[2018 election data](https://redistrictingdatahub.org/dataset/vest-2018-texas-precinct-boundaries-and-election-results/)**:**  VEST 2018 Texas precinct and election results

[2016 election data](https://redistrictingdatahub.org/dataset/vest-2016-texas-precinct-and-election-results/)**:**  VEST 2016 Texas precinct and election results

[2020 County data](https://redistrictingdatahub.org/dataset/texas-county-pl-94171-2020/):  from 2020 Census Redistricting Data (P.L. 94-171) Shapefiles

# **Processing**

Data were cleaned and aggregated in the corresponding jupyter notebook using MGGG’s python library [maup](https://github.com/mggg/maup). 

## Metadata
- `CNTY18`: County number
- `COLOR18`: Column from Vest 2020 election data (has range of integers from 1 to 7)
- `PREC18`: Precinct number
- `PCTKEY18`: Precinct FIPS code from 2018
- `CNTYKEY18`: County code
- `G20VR18`: Voter registration
- `G20SSVR18`: Spanish surname voter registration
- `CD`: Congressional district for 2021 enacted plan
- `SEND`: State Senate district for 2021 State Senate Adopted Plan
- `HDIST`: State House district for 2021 State House of Representatives Districts Plan
- `TOTPOP`: Total population in 2020 Census
- `NH_WHITE`: White, non-hispanic, population in 2020 Census
- `NH_BLACK`: Black, non-hispanic, population in 2020 Census
- `NH_AMIN`: American Indian and Alaska Native, non-hispanic, population in 2020 Census
- `NH_ASIAN`: Asian, non-hispanic, population in 2020 Census
- `NH_NHPI`: Native Hawaiian and Pacific Islander, non-hispanic, population in 2020 Census
- `NH_OTHER`: Other race, non-hispanic, population in 2020 Census
- `NH_2MORE`: Two or more races, non-hispanic, population in 2020 Census
- `HISP`: Hispanic population in 2020 Census
- `H_WHITE`: White, hispanic, population in 2020 Census
- `H_BLACK`: Black, hispanic, population in 2020 Census
- `H_AMIN`: American Indian and Alaska Native, hispanic, population in 2020 Census
- `H_ASIAN`: Asian, hispanic, population in 2020 Census
- `H_NHPI`: Native Hawaiian and Pacific Islander, hispanic, population in 2020 Census
- `H_OTHER`: Other race, hispanic, population in 2020 Census
- `H_2MORE`: Two or more races, hispanic, population in 2020 Census
- `VAP`: Total voting age population in 2020 Census
- `HVAP`: Hispanic voting age population in 2020 Census
- `WVAP`: White, non-hispanic, voting age population in 2020 Census
- `BVAP`: Black, non-hispanic, voting age population in 2020 Census
- `AMINVAP`: American Indian and Alaska Native, non-hispanic, voting age population in 2020 Census
- `ASIANVAP`: Asian, non-hispanic, voting age population in 2020 Census
- `NHPIVAP`: Native Hawaiian and Pacific Islander, non-hispanic, voting age population in 2020 Census
- `OTHERVAP`: Other race, non-hispanic, voting age population in 2020 Census
- `2MOREVAP`: Two or more races, non-hispanic, voting age population in 2020 Census
- `AGR18D`: Number of votes for 2018 Democratic agriculture commissioner candidate
- `AGR18R`: Number of votes for 2018 Republican agriculture commissioner candidate
- `AGR18O`: Number of votes for 2018 other party's agriculture commissioner candidate
- `ATG18D`: Number of votes for 2018 Democratic attorney general candidate
- `ATG18R`: Number of votes for 2018 Republican attorney general candidate
- `ATG18O`: Number of votes for 2018 other party's attorney general candidate
- `COM18D`: Number of votes for 2018 Democratic comptroller candidate
- `COM18R`: Number of votes for 2018 Republican comptroller candidate
- `COM18O`: Number of votes for 2018 other party's comptroller candidate
- `GOV18D`: Number of votes for 2018 Democratic gubernatorial candidate
- `GOV18R`: Number of votes for 2018 Republican gubernatorial candidate
- `GOV18O`: Number of votes for 2018 other party's gubernatorial candidate
- `LAN18D`: Number of votes for 2018 Democratic Commissioner of General Land Office candidate
- `LAN18R`: Number of votes for 2018 Republican Commissioner of General Land Office candidate
- `LAN18O`: Number of votes for 2018 other party's Commissioner of General Land Office candidate
- `LTG18D`: Number of votes for 2018 Democratic Lieutenant Governor candidate
- `LTG18R`: Number of votes for 2018 Republican Lieutenant Governor candidate
- `LTG18O`: Number of votes for 2018 other party's Lieutenant Governor candidate
- `PRE16D`: Number of votes for 2016 Democratic President
- `PRE16R`: Number of votes for 2016 Republican President
- `PRE16O`: Number of votes for 2016 other party's President
- `PRE20D`: Number of votes for 2020 Democratic President
- `PRE20R`: Number of votes for 2020 Republican President
- `PRE20O`: Number of votes for 2020 other party's President
- `RRC16D`: Number of votes for 2016 Democratic railroad commissioner candidate
- `RRC16R`: Number of votes for 2016 Republican railroad commissioner candidate
- `RRC16O`: Number of votes for 2016 other party's railroad commissioner candidate
- `RRC18D`: Number of votes for 2018 Democratic railroad commissioner candidate
- `RRC18R`: Number of votes for 2018 Republican railroad commissioner candidate
- `RRC18O`: Number of votes for 2018 other party's railroad commissioner candidate
- `RRC20D`: Number of votes for 2020 Democratic railroad commissioner candidate
- `RRC20R`: Number of votes for 2020 Republican railroad commissioner candidate
- `RRC20O`: Number of votes for 2020 other party's railroad commissioner candidate
- `SCC16D`: Number of votes for 2016 Democratic State Court of Criminal Appeals candidate
- `SCC16R`: Number of votes for 2016 Republican State Court of Criminal Appeals candidate
- `SCC16O`: Number of votes for 2016 other party's State Court of Criminal Appeals candidate
- `SCC18D`: Number of votes for 2018 Democratic State Court of Criminal Appeals candidate
- `SCC18R`: Number of votes for 2018 Republican State Court of Criminal Appeals candidate
- `SCC18O`: Number of votes for 2018 other party's State Court of Criminal Appeals candidate
- `SCC20D`: Number of votes for 2020 Democratic State Court of Criminal Appeals candidate
- `SCC20R`: Number of votes for 2020 Republican State Court of Criminal Appeals candidate
- `SSC16D`: Number of votes for 2016 Democratic State Supreme Court candidate
- `SSC16R`: Number of votes for 2016 Republican State Supreme Court candidate
- `SSC16O`: Number of votes for 2016 other party's State Supreme Court candidate
- `SSC18D`: Number of votes for 2018 Democratic State Supreme Court candidate
- `SSC18R`: Number of votes for 2018 Republican State Supreme Court candidate
- `SSC20D`: Number of votes for 2020 Democratic State Supreme Court candidate
- `SSC20R`: Number of votes for 2020 Republican State Supreme Court candidate
- `SSC20O`: Number of votes for 2020 other party's State Supreme Court candidate
- `USS18D`: Number of votes for 2018 Democratic senate candidate
- `USS18R`: Number of votes for 2018 Republican senate candidate
- `USS18O`: Number of votes for 2018 other party's senate candidate
- `USS20D`: Number of votes for 2020 Democratic senate candidate
- `USS20R`: Number of votes for 2020 Republican senate candidate
- `USS20O`: Number of votes for 2020 other party's senate candidate

